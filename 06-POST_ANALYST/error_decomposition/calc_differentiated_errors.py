"""
calc_differentiated_errors.py — v15.0 Fehlerzerlegung nach Disk / Avg / Int / Ref

Thesis-Notation (sauberer Cut, keine Legacy-Brücken):

    rho_Total = rho_Disk · rho_Avg · rho_Int · rho_Ref

Fehlerklassen:
    Disk : Diskretisierung des FE-Spannungsfelds (FE-Linearisierung)
    Avg  : Nodale Mittelung (=1 bei STRESS=RAW)
    Int  : Integrationsverfahren (=1 bei INT=G26)
    Ref  : Wahl der Referenzspannung (=1 bei Pf wegen smax-Kürzung)
    Total: Gesamt = Q_method / Q_ana

Metriken pro Klasse:
    rho_X   = Q_x / Q_(x-1)        multiplikativer Faktor
    delta_X = rho_X - 1.0          relative Abweichung
    Delta_X = absolute Differenz   (siehe Definition pro Klasse)

Spalten-Schema in _extended.csv:
    Q_rho_Disk   Q_delta_Disk   Q_Delta_Disk
    Q_rho_Avg    Q_delta_Avg    Q_Delta_Avg
    Q_rho_Int    Q_delta_Int    Q_Delta_Int
    Q_rho_Ref    Q_delta_Ref    Q_Delta_Ref
    Q_rho_Total  Q_delta_Total  Q_Delta_Total

mit Q ∈ {Veff, Aeff, Pf}.

Public API:
    compute_v15_errors(df) -> DataFrame mit zusätzlichen Spalten
    run_error_analysis()   -> Standalone-Workflow (CLI)
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

# parse_case_id stammt aus method_naming.py (Single Source of Truth, helpers/).
# Frueher kam es aus 08-PLOTS/case_id_legend.py; dieses Modul (matplotlib-Labels/Farben
# fuer die Plot-Skripte) ist im Public-Release entfallen — parse_case_id ist identisch
# in helpers/method_naming.py vorhanden.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'helpers'))
import method_naming as mn

# === KONFIGURATION (nur für Standalone-Modus) ===
REFERENCE_CSV = os.path.join("06-POST_ANALYST", "REFERENCES", "Reference_Gauss26.csv")
ANALYSIS_RUN = r"05-RUNS/2026-04-29_23-19"
OUTPUT_SUFFIX = "_extended"

FILTER_MESH = []
FILTER_GEOM = []
FILTER_METHODS = []     # z.B. ["RAW-G5-NDX", "AVG-G5-GPX"]


# ----------------------------------------------------------------------
# v15-Konstanten
# ----------------------------------------------------------------------

ERROR_COMPONENTS = ["Disk", "Avg", "Int", "Ref", "Total"]
ERROR_METRICS = ["rho", "delta", "Delta"]


# ----------------------------------------------------------------------
# Hauptberechnungs-Funktion
# ----------------------------------------------------------------------

def compute_v15_errors(df, stress_token, int_token):
    """v15.0: Berechnet Disk/Avg/Int/Ref-Fehlerzerlegung pro Zielgröße.

    Erwartet im DataFrame folgende Spalten (je nach LOC):
        m, S_num, S_ana
        Veff_num, Veff_ana   (LOC=V)
        Aeff_num, Aeff_ana   (LOC=A)
        Pf_num, Pf_ana
        RAW_V_g26, RAW_A_g26, RAW_Pf_g26, RAW_S_g26   (immer benötigt)
        AVG_V_g26, AVG_A_g26, AVG_Pf_g26, AVG_S_g26   (nur bei stress_token=='AVG')

    Erzeugt pro Zielgröße Q (15 Spalten):
        Q_rho_Disk,   Q_delta_Disk,   Q_Delta_Disk
        Q_rho_Avg,    Q_delta_Avg,    Q_Delta_Avg
        Q_rho_Int,    Q_delta_Int,    Q_Delta_Int
        Q_rho_Ref,    Q_delta_Ref,    Q_Delta_Ref
        Q_rho_Total,  Q_delta_Total,  Q_Delta_Total

    Parameters
    ----------
    df : pd.DataFrame
        Gemergt mit RAW- und ggf. AVG-G26-Referenz.
    stress_token : str
        'RAW' oder 'AVG' (aus parse_case_id info['stress']).
    int_token : str
        'EM', 'G1'..'G9', 'G26' (aus parse_case_id info['int']).

    Returns
    -------
    pd.DataFrame
        Gleicher DataFrame mit zusätzlichen v15-Spalten.

    Raises
    ------
    RuntimeError
        Wenn benötigte G26-Referenz-Spalten fehlen oder komplett NaN sind.
    """
    s = stress_token.strip().upper()
    i = int_token.strip().upper()

    if s not in ("RAW", "AVG"):
        raise ValueError(f"Unbekanntes stress_token: {s!r} (erwartet RAW oder AVG)")

    # Validierung: RAW-G26-Referenz muss immer da sein
    if 'RAW_S_g26' not in df.columns or df['RAW_S_g26'].isna().all():
        raise RuntimeError(
            "v15-ERR_EXT: RAW-G26-Referenz fehlt. "
            "Erzeuge G26-Run mit STRESS=RAW + INT=G26 + REF=NDX und regeneriere "
            "Reference_Gauss26.csv via 06-generate_gauss26_reference.py."
        )

    # Validierung: AVG-Methode braucht zusätzlich AVG-G26-Referenz
    if s == "AVG":
        if 'AVG_S_g26' not in df.columns or df['AVG_S_g26'].isna().all():
            raise RuntimeError(
                "v15-ERR_EXT: STRESS=AVG-Methode benötigt AVG-G26-Referenz, aber "
                "keine AVG-G26-Eintrag in Reference_Gauss26.csv. "
                "Erzeuge G26-Run mit STRESS=AVG + INT=G26 + REF=NDX und regeneriere "
                "Reference_Gauss26.csv."
            )

    m_vals = df['m'].astype(float)
    S_a = df['S_ana'].astype(float)
    S_n = df['S_num'].astype(float)
    RAW_S = df['RAW_S_g26'].astype(float)
    AVG_S = df['AVG_S_g26'].astype(float) if s == "AVG" else None

    # Targets dynamisch (LOC=V → Veff, LOC=A → Aeff, immer Pf)
    target_specs = []
    if 'Veff_num' in df.columns and 'Veff_ana' in df.columns and 'RAW_V_g26' in df.columns:
        target_specs.append(('Veff', 'Veff_num', 'Veff_ana', 'RAW_V_g26', 'AVG_V_g26', False))
    if 'Aeff_num' in df.columns and 'Aeff_ana' in df.columns and 'RAW_A_g26' in df.columns:
        target_specs.append(('Aeff', 'Aeff_num', 'Aeff_ana', 'RAW_A_g26', 'AVG_A_g26', False))
    if 'Pf_num' in df.columns and 'Pf_ana' in df.columns and 'RAW_Pf_g26' in df.columns:
        target_specs.append(('Pf', 'Pf_num', 'Pf_ana', 'RAW_Pf_g26', 'AVG_Pf_g26', True))

    if not target_specs:
        raise RuntimeError(
            "v15-ERR_EXT: Keine vollständigen Spalten-Sets für Veff/Aeff/Pf gefunden."
        )

    is_g26 = (i == "G26")

    for target, num_col, ana_col, raw_g26_col, avg_g26_col, is_pf in target_specs:
        Q_method = df[num_col].astype(float)
        Q_ana = df[ana_col].astype(float)
        Q_raw_g26 = df[raw_g26_col].astype(float)
        Q_avg_g26 = df[avg_g26_col].astype(float) if (s == "AVG" and avg_g26_col in df.columns) else None

        # --- Korrigierte Werte (Ref-Korrektur) ---
        if is_pf:
            # Pf: smax kürzt sich strukturell in Weibull heraus -> keine Korrektur
            Q_method_corr = Q_method.copy()
            Q_raw_g26_corr = Q_raw_g26.copy()
            Q_avg_g26_corr = Q_avg_g26.copy() if Q_avg_g26 is not None else None
        else:
            # Veff/Aeff: Q_corr = Q · (S_used/S_ana)^m
            Q_method_corr = Q_method * np.power(S_n / S_a, m_vals)
            Q_raw_g26_corr = Q_raw_g26 * np.power(RAW_S / S_a, m_vals)
            Q_avg_g26_corr = (Q_avg_g26 * np.power(AVG_S / S_a, m_vals)
                              if Q_avg_g26 is not None else None)

        # --- Disk-Anteil ---
        rho_Disk = Q_raw_g26_corr / Q_ana
        Delta_Disk = Q_raw_g26_corr - Q_ana

        # --- Avg-Anteil ---
        if s == "RAW":
            rho_Avg = pd.Series(1.0, index=df.index)
            Delta_Avg = pd.Series(0.0, index=df.index)
        else:
            rho_Avg = Q_avg_g26_corr / Q_raw_g26_corr
            Delta_Avg = Q_avg_g26_corr - Q_raw_g26_corr

        # --- Int-Anteil ---
        Q_stress_g26_corr = Q_avg_g26_corr if s == "AVG" else Q_raw_g26_corr
        if is_g26:
            # Methode IST die G26-Referenz
            rho_Int = pd.Series(1.0, index=df.index)
            Delta_Int = pd.Series(0.0, index=df.index)
        else:
            rho_Int = Q_method_corr / Q_stress_g26_corr
            Delta_Int = Q_method_corr - Q_stress_g26_corr

        # --- Ref-Anteil ---
        if is_pf:
            # Pf: smax kürzt sich -> rho_Ref = 1, Delta_Ref = 0 (explizit)
            rho_Ref = pd.Series(1.0, index=df.index)
            Delta_Ref = pd.Series(0.0, index=df.index)
        else:
            rho_Ref = Q_method / Q_method_corr
            Delta_Ref = Q_method - Q_method_corr

        # --- Total ---
        rho_Total = Q_method / Q_ana
        Delta_Total = Q_method - Q_ana

        # --- delta = rho - 1 ---
        delta_Disk = rho_Disk - 1.0
        delta_Avg = rho_Avg - 1.0
        delta_Int = rho_Int - 1.0
        delta_Ref = rho_Ref - 1.0
        delta_Total = rho_Total - 1.0

        # --- Spalten schreiben ---
        df[f'{target}_rho_Disk'] = rho_Disk
        df[f'{target}_delta_Disk'] = delta_Disk
        df[f'{target}_Delta_Disk'] = Delta_Disk

        df[f'{target}_rho_Avg'] = rho_Avg
        df[f'{target}_delta_Avg'] = delta_Avg
        df[f'{target}_Delta_Avg'] = Delta_Avg

        df[f'{target}_rho_Int'] = rho_Int
        df[f'{target}_delta_Int'] = delta_Int
        df[f'{target}_Delta_Int'] = Delta_Int

        df[f'{target}_rho_Ref'] = rho_Ref
        df[f'{target}_delta_Ref'] = delta_Ref
        df[f'{target}_Delta_Ref'] = Delta_Ref

        df[f'{target}_rho_Total'] = rho_Total
        df[f'{target}_delta_Total'] = delta_Total
        df[f'{target}_Delta_Total'] = Delta_Total

    # --- Produktform-Validierung ---
    _validate_product_form(df, target_specs, tolerance_log=1e-6)

    return df


def _validate_product_form(df, target_specs, tolerance_log=1e-6):
    """Validiert numerisch: rho_Total ≈ rho_Disk · rho_Avg · rho_Int · rho_Ref.

    Vergleicht im log-space (robust gegen Vorzeichen). Pf-Validierung mit
    eigener Toleranz nur bei Pf_num > 1e-10 (vermeidet log-Instabilität).
    """
    for target, _num, _ana, _raw, _avg, is_pf in target_specs:
        rho_t = df[f'{target}_rho_Total'].astype(float)
        prod = (df[f'{target}_rho_Disk'].astype(float)
                * df[f'{target}_rho_Avg'].astype(float)
                * df[f'{target}_rho_Int'].astype(float)
                * df[f'{target}_rho_Ref'].astype(float))

        # Filter: nur wo beide positiv und nicht-trivial
        if is_pf:
            mask = (rho_t > 1e-10) & (prod > 1e-10)
        else:
            mask = (rho_t > 0) & (prod > 0)

        if not mask.any():
            print(f"   [ErrDecomp v15] {target}: kein gültiger Vergleich möglich (alle Werte 0/NaN)")
            continue

        log_diff = np.abs(np.log(rho_t[mask]) - np.log(prod[mask]))
        max_err = float(log_diff.max())

        if max_err < tolerance_log:
            print(f"   [ErrDecomp v15] {target}: Produktform OK (max log-Diff = {max_err:.2e})")
        else:
            print(f"   [ErrDecomp v15] WARNUNG: {target} Produktform-Bruch! "
                  f"max log-Diff = {max_err:.2e} (Toleranz {tolerance_log:.0e})")


# ----------------------------------------------------------------------
# Standalone-Workflow (CLI)
# ----------------------------------------------------------------------

def run_error_analysis():
    """Standalone-Workflow v15: Lädt Reference + ANALYSIS-CSVs, berechnet Fehler."""
    if not os.path.exists(REFERENCE_CSV):
        print(f"FEHLER: Reference-Datei nicht gefunden: {REFERENCE_CSV}")
        print("Bitte zuerst 06-generate_gauss26_reference.py ausführen!")
        return

    df_ref = pd.read_csv(REFERENCE_CSV, delimiter=';')
    print(f"Reference: {len(df_ref)} Zeilen aus {REFERENCE_CSV}")

    run_path = Path(ANALYSIS_RUN)
    if not run_path.exists():
        print(f"FEHLER: Run-Ordner nicht gefunden: {ANALYSIS_RUN}")
        return

    all_csvs = list(run_path.rglob("ANALYSIS_*.csv"))
    all_csvs = [p for p in all_csvs if OUTPUT_SUFFIX not in p.stem]
    print(f"Gefunden: {len(all_csvs)} ANALYSIS-Dateien in {ANALYSIS_RUN}")

    processed = 0
    skipped = 0
    failed = 0

    for csv_path in all_csvs:
        case_id_str = csv_path.stem.replace("ANALYSIS_", "")
        info = mn.parse_case_id(case_id_str)
        if "error" in info:
            print(f"  ! Konnte Case-ID nicht parsen: {csv_path.stem}")
            continue

        if FILTER_MESH and info["mesh"] not in FILTER_MESH:
            continue
        geom_key = f"{info['case']}-{info['geom']}"
        if FILTER_GEOM and geom_key not in FILTER_GEOM:
            continue
        if FILTER_METHODS and info["method_str"] not in FILTER_METHODS:
            continue

        try:
            df = pd.read_csv(csv_path, delimiter=';')
            # JSON-Metadata fuer LOAD_N/SIG_0_ABS/SIG_0 laden,
            # damit der G26-Merge nur passende Reference-Zeilen findet.
            full_name = csv_path.stem.replace("ANALYSIS_", "")
            json_path = csv_path.parent / f"{full_name}.json"
            if not json_path.exists():
                print(f"  ✗ {csv_path.stem}: JSON-Metadata fehlt ({json_path.name}); v20-Pflicht.")
                failed += 1
                continue
            with open(json_path, 'r', encoding='utf-8') as jf:
                meta = json.load(jf)
            curr_load_n    = float(meta['load_n'])
            curr_sig_0_abs = bool(meta['sig_0_abs'])
            curr_sig_0     = float(meta['sig_0'])
            curr_pia_fix   = bool(meta.get('pia_fix', False))
            df_merged = _merge_dual_g26_refs(df, df_ref, info,
                                             curr_load_n, curr_sig_0_abs, curr_sig_0,
                                             pia_fix=curr_pia_fix)
            df_result = compute_v15_errors(df_merged, info["stress"], info["int"])

            # Hilfsspalten entfernen
            drop_cols = ['Mesh', 'Loadcase_self', 'Loadcase_other']
            df_result = df_result.drop(columns=[c for c in drop_cols if c in df_result.columns])

            output_path = csv_path.with_name(csv_path.stem + OUTPUT_SUFFIX + ".csv")
            df_result.to_csv(output_path, sep=';', index=False)
            print(f"  ✓ {csv_path.stem} -> {output_path.name}")
            processed += 1

        except RuntimeError as e:
            print(f"  ✗ {csv_path.stem}: {e}")
            failed += 1
        except Exception as e:
            print(f"  ✗ {csv_path.stem}: Unerwarteter Fehler: {e}")
            failed += 1

    print(f"\nZusammenfassung: {processed} verarbeitet, {skipped} übersprungen, {failed} fehlgeschlagen")


def _merge_dual_g26_refs(df_analysis, df_ref, info, load_n, sig_0_abs, sig_0, pia_fix=False):
    """v15+v20.0+v20.2: Mergt ANALYSIS-DataFrame mit RAW- und (bei AVG) AVG-G26-Referenzen.

    Setzt Mesh, Loadcase_self, Loadcase_other-Hilfsspalten und führt zwei
    Left-Joins aus. Returns den gemergten DataFrame mit Spalten:
        RAW_V_g26, RAW_A_g26, RAW_Pf_g26, RAW_S_g26
        AVG_V_g26, AVG_A_g26, AVG_Pf_g26, AVG_S_g26  (nur bei STRESS=AVG)

    v20.0 (Bug #2): Pre-Filter df_ref auf passende LOAD_N/SIG_0_ABS/SIG_0 vor dem Merge.
    Verhindert dass Cases mit gleicher Geometrie aber unterschiedlicher Last/Sigma0
    falsch gegeneinander gemergt werden.

    v20.2: Hard-Fail wenn pia_fix=False. G26-Reference ist in v20.2 ausschliesslich
    fuer PIAFix-Cases gueltig (siehe 06-generate_gauss26_reference.py v20.2-Filter).
    Legacy-Cases mit ERR_EXT=1 muessen (a) auf PIA_FIX=1 umstellen oder (b) ERR_EXT=0
    setzen (kein Decomposition-Output, aber Pipeline laeuft durch).
    """
    if not pia_fix:
        raise RuntimeError(
            "v20.2: G26-Decomposition setzt PIA_FIX=1 voraus. Legacy-Pipeline "
            "(PIA_FIX=0) integriert PIA an Knoten + Interpolation — methodisch "
            "nicht aequivalent zur PIAFix-Variante (PIA-at-GP). Eine Decomposition "
            "gegen die PIAFix-G26-Reference waere ein methodischer Hybrid.\n"
            "Aktion: (1) PIA_FIX=1 in CSV setzen + Re-Run, oder "
            "(2) ERR_EXT=0 setzen (Basis-CSV ohne Decomposition)."
        )
    own_stress = info["stress"]
    other_stress = "AVG" if own_stress == "RAW" else "RAW"

    key_self = f"{info['case']}-{info['geom']}-{info['crit']}-{info['loc']}-{own_stress}"
    key_other = f"{info['case']}-{info['geom']}-{info['crit']}-{info['loc']}-{other_stress}"

    # df_ref auf LOAD_N/SIG_0_ABS/SIG_0 filtern, bevor wir mergen.
    required_cols = {'LOAD_N', 'SIG_0_ABS', 'SIG_0'}
    missing_cols = required_cols - set(df_ref.columns)
    if missing_cols:
        raise RuntimeError(
            f"v20.0 G26-Reference ({REFERENCE_CSV}) hat altes Schema (v18.x). "
            f"Fehlende Spalten: {missing_cols}. "
            f"Aktion: alte Reference loeschen + neu generieren mit "
            f"06-generate_gauss26_reference.py."
        )
    df_ref_filtered = df_ref[
        (df_ref['LOAD_N'].astype(float) == float(load_n)) &
        (df_ref['SIG_0_ABS'].astype(bool) == bool(sig_0_abs)) &
        (df_ref['SIG_0'].astype(float) == float(sig_0))
    ].copy()
    if df_ref_filtered.empty:
        print(f"  [Merge] WARN: kein Reference-Eintrag fuer LOAD_N={load_n} "
              f"SIG_0_ABS={sig_0_abs} SIG_0={sig_0} — alle G26-Spalten werden NaN.")

    df = df_analysis.copy()
    df['Mesh'] = info['mesh']
    df['Loadcase_self'] = key_self

    # Erst-Merge: own_stress G26-Referenz
    rename_self = {
        'V_mesh': f'{own_stress}_V_g26',
        'A_mesh': f'{own_stress}_A_g26',
        'Pf_mesh': f'{own_stress}_Pf_g26',
        'S_mesh': f'{own_stress}_S_g26',
    }
    df_self = df_ref_filtered.rename(columns={'Loadcase': 'Loadcase_self'})[
        ['Mesh', 'Loadcase_self', 'm'] + [c for c in ['V_mesh', 'A_mesh', 'Pf_mesh', 'S_mesh']
                                           if c in df_ref.columns]
    ]
    df_self = df_self.rename(columns=rename_self)
    df_merged = df.merge(df_self, on=['Mesh', 'Loadcase_self', 'm'], how='left')

    # Zweit-Merge: other_stress (bei STRESS=AVG → RAW-Ref für Avg-Anteil; bei STRESS=RAW → AVG-Ref optional)
    df_merged['Loadcase_other'] = key_other
    rename_other = {
        'V_mesh': f'{other_stress}_V_g26',
        'A_mesh': f'{other_stress}_A_g26',
        'Pf_mesh': f'{other_stress}_Pf_g26',
        'S_mesh': f'{other_stress}_S_g26',
    }
    df_other = df_ref_filtered.rename(columns={'Loadcase': 'Loadcase_other'})[
        ['Mesh', 'Loadcase_other', 'm'] + [c for c in ['V_mesh', 'A_mesh', 'Pf_mesh', 'S_mesh']
                                            if c in df_ref.columns]
    ]
    df_other = df_other.rename(columns=rename_other)
    df_merged = df_merged.merge(df_other, on=['Mesh', 'Loadcase_other', 'm'], how='left')

    return df_merged


if __name__ == "__main__":
    run_error_analysis()
