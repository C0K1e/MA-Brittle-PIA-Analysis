"""
Skript A: Reference Aggregator (v8.0)

Liest ANALYSIS-CSVs aus einem G26-Run-Ordner und aggregiert sie
zu einer zentralen Reference_Gauss26.csv.

Workflow:
    1. User fuehrt G26-Lauf manuell durch (00-cases.csv mit INT=G26)
    2. ANALYSIS_*.csv liegen in 05-RUNS/{timestamp}/{case}/tables/
    3. Dieses Skript liest alle ANALYSIS-CSVs, extrahiert Veff_num + Pf_num + S_num
    4. Umbenennung zu V_mesh / Pf_mesh / S_mesh (+ Legacy-Aliases)
    5. Aggregation zu Reference_Gauss26.csv

Output-Spalten (v8.0): Mesh;Loadcase;m;V_mesh;Pf_mesh;S_mesh;VF_Gauss26;PF_Gauss26
Output-Spalten (Legacy, ohne S_num in Quelle): Mesh;Loadcase;m;V_mesh;Pf_mesh;VF_Gauss26;PF_Gauss26
"""

import os
import sys
import json
import pandas as pd
from pathlib import Path

# parse_case_id stammt aus method_naming.py (Single Source of Truth, helpers/).
# Frueher kam es aus 08-PLOTS/case_id_legend.py; dieses Modul (matplotlib-Labels/Farben
# fuer die Plot-Skripte) ist im Public-Release entfallen — parse_case_id ist identisch
# in helpers/method_naming.py vorhanden.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'helpers'))
import method_naming as mn

# === KONFIGURATION ===
SOURCE_RUN = r"05-RUNS/2026-05-16_00-09"  # Ordner mit G26-Ergebnissen (anpassen!)
OUTPUT_CSV = os.path.join("06-POST_ANALYST", "REFERENCES", "Reference_Gauss26.csv")  # Ausgabe-Datei

# --- FILTER (optional, leer = alle) ---
FILTER_MESH = []    # z.B. ["4x6", "8x12"] — nur bestimmte Meshes
FILTER_GEOM = []    # z.B. ["PWH-60x0x1"] — nur bestimmte Geometrien


def generate_reference():
    """Aggregiert G26-ANALYSIS-CSVs zu Reference_Gauss26.csv."""
    run_path = Path(SOURCE_RUN)
    if not run_path.exists():
        print(f"FEHLER: Run-Ordner nicht gefunden: {SOURCE_RUN}")
        return

    all_csvs = list(run_path.rglob("ANALYSIS_*.csv"))
    # Bereits erweiterte CSVs ausschliessen
    all_csvs = [p for p in all_csvs if '_extended' not in p.stem]
    print(f"Gefunden: {len(all_csvs)} ANALYSIS-Dateien in {SOURCE_RUN}")

    records = []
    has_s_num = None  # Wird beim ersten CSV gesetzt

    n_skipped_int = 0
    n_skipped_ref = 0
    n_skipped_stress = 0
    n_skipped_pia_fix = 0  # Legacy-Cases ohne PIA_FIX werden nicht als Reference akzeptiert

    for csv_path in all_csvs:
        # Case-ID parsen (ANALYSIS_ Prefix entfernen)
        case_id_str = csv_path.stem.replace("ANALYSIS_", "")
        info = mn.parse_case_id(case_id_str)
        if "error" in info:
            print(f"  ! Konnte Case-ID nicht parsen: {csv_path.stem}")
            continue

        # Strikte Validierung — nur echte G26-Reference-Cases akzeptieren.
        # Frueher konnten faelschlich INT=G5/G9-Cases im Source-Run mitaggregiert werden.
        if info["int"] != "G26":
            n_skipped_int += 1
            continue
        if info["ref"] != "NDX":
            n_skipped_ref += 1
            continue
        if info["stress"] not in ("RAW", "AVG"):
            n_skipped_stress += 1
            continue

        # Filter anwenden
        if FILTER_MESH and info["mesh"] not in FILTER_MESH:
            continue
        geom_key = f"{info['case']}-{info['geom']}"
        if FILTER_GEOM and geom_key not in FILTER_GEOM:
            continue

        # CSV laden
        df = pd.read_csv(csv_path, delimiter=';')

        # JSON-Metadata neben ANALYSIS-CSV laden, um LOAD_N/SIG_0_ABS/SIG_0/sigma0
        # in den Reference-Key aufzunehmen. Verhindert dass zwei Cases mit gleicher Geometrie/
        # Mesh aber unterschiedlicher Last/Sigma0 dieselbe Reference-Zeile ueberschreiben.
        # JSON-Pfad: tables/{full_name}.json (vom report_generator.py geschrieben)
        full_name = csv_path.stem.replace("ANALYSIS_", "")
        json_path = csv_path.parent / f"{full_name}.json"
        if not json_path.exists():
            print(f"!!! [G26-Gen] FATAL: JSON-Metadata fehlt fuer {full_name}")
            print(f"!!! [G26-Gen]        Erwartet unter: {json_path}")
            print(f"!!! [G26-Gen]        Aktion: Case neu rechnen (ANALYSIS-CSV ohne JSON ist v18-Stand).")
            sys.exit(1)
        with open(json_path, 'r', encoding='utf-8') as jf:
            meta = json.load(jf)
        required_meta = {'load_n', 'sig_0_abs', 'sig_0', 'sigma0'}
        missing_meta = required_meta - set(meta.keys())
        if missing_meta:
            print(f"!!! [G26-Gen] FATAL: JSON-Metadata fuer {full_name} fehlt {missing_meta}")
            print(f"!!! [G26-Gen]        Aktion: Case neu rechnen, damit JSON v20-aktuell ist.")
            sys.exit(1)
        load_n_val    = float(meta['load_n'])
        sig_0_abs_val = bool(meta['sig_0_abs'])
        sig_0_val     = float(meta['sig_0'])
        sigma0_val    = float(meta['sigma0'])

        # nur Cases mit PIA_FIX=1 als G26-Reference akzeptieren.
        # Legacy-Cases (PIA_FIX=0) integrieren PIA an Knoten + interpolation —
        # methodisch nicht aequivalent zur PIAFix-Variante. Eine gemischte Reference
        # waere falsch (verschiedene Integrale gegen unterschiedliche Cases gemerged).
        pia_fix_val = bool(meta.get('pia_fix', False))
        if not pia_fix_val:
            n_skipped_pia_fix += 1
            continue

        # Loadcase = CASE-GEOM-CRIT-LOC-STRESS (STRESS={RAW,AVG} differenziert Spannungsfelder)
        loadcase = f"{info['case']}-{info['geom']}-{info['crit']}-{info['loc']}-{info['stress']}"

        # Relevante Spalten extrahieren + umbenennen — domain-aware (V vs A)
        cols_to_extract = ['m']
        rename_map = {}

        # Volumen-Pfad (LOC=V) — Standard
        if 'Veff_num' in df.columns:
            cols_to_extract.append('Veff_num')
            rename_map['Veff_num'] = 'V_mesh'
        # Surface-Pfad (LOC=A, v12.0)
        if 'Aeff_num' in df.columns:
            cols_to_extract.append('Aeff_num')
            rename_map['Aeff_num'] = 'A_mesh'
        # Pf (LOC-unabhaengig)
        if 'Pf_num' in df.columns:
            cols_to_extract.append('Pf_num')
            rename_map['Pf_num'] = 'Pf_mesh'

        # S_num extrahieren wenn vorhanden
        if 'S_num' in df.columns:
            cols_to_extract.append('S_num')
            rename_map['S_num'] = 'S_mesh'
            if has_s_num is None:
                has_s_num = True
                print("  [INFO] S_num Spalte gefunden -> S_mesh wird extrahiert")
        elif has_s_num is None:
            has_s_num = False
            print("  [INFO] Kein S_num in CSV -> S_mesh nicht verfuegbar (Legacy-Modus)")

        ref = df[cols_to_extract].copy()
        ref = ref.rename(columns=rename_map)

        # Legacy-Aliases fuer Rueckwaertskompatibilitaet (nur Volumen)
        if 'V_mesh' in ref.columns:
            ref['VF_Gauss26'] = ref['V_mesh']
        if 'Pf_mesh' in ref.columns:
            ref['PF_Gauss26'] = ref['Pf_mesh']

        ref['Mesh'] = info['mesh']
        ref['Loadcase'] = loadcase
        # Last/Sigma0-Disambiguierung im Reference-Key
        ref['LOAD_N']    = load_n_val
        ref['SIG_0_ABS'] = sig_0_abs_val
        ref['SIG_0']     = sig_0_val
        ref['sigma0']    = sigma0_val
        records.append(ref)

        print(f"  + {info['mesh']:>10s} | {loadcase:<30s} | LOAD_N={load_n_val:g} sig0={sigma0_val:g} | {len(df)} m-Werte")

    # Validierungs-Statistik ausgeben
    if n_skipped_int + n_skipped_ref + n_skipped_stress + n_skipped_pia_fix > 0:
        print(f"  [Validation] Skipped Cases: INT≠G26: {n_skipped_int}, "
              f"REF≠NDX: {n_skipped_ref}, STRESS∉(RAW,AVG): {n_skipped_stress}, "
              f"PIA_FIX=0: {n_skipped_pia_fix}")
        print(f"  [Validation] Erwartet werden nur INT=G26 + REF=NDX + STRESS∈{{RAW,AVG}} + PIA_FIX=1-Cases.")
        print(f"  [Validation] v20.2: Legacy-PIA-Knoten-Aggregation (PIA_FIX=0) ist methodisch")
        print(f"  [Validation]        nicht aequivalent zu PIAFix (PIA-at-GP) — kein Reference-Merge.")

    if not records:
        print("FEHLER: Keine passenden ANALYSIS-Dateien gefunden!")
        if n_skipped_int + n_skipped_ref + n_skipped_stress > 0:
            print("       Hinweis: Alle gefundenen Cases haben unzulaessige INT/REF/STRESS-Werte.")
        return

    # Aggregieren (nur neue Daten aus diesem Lauf)
    df_new = pd.concat(records, ignore_index=True)

    # Upsert: Bestehende Reference-Datei laden und zusammenfuehren
    # Key erweitert um LOAD_N/SIG_0_ABS/SIG_0 — verhindert dass
    # Cases mit gleicher Geometrie aber unterschiedlicher Last sich ueberschreiben.
    key_cols = ['Mesh', 'Loadcase', 'LOAD_N', 'SIG_0_ABS', 'SIG_0', 'm']
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    if os.path.exists(OUTPUT_CSV):
        df_existing = pd.read_csv(OUTPUT_CSV, delimiter=';')
        # Schema-Validierung: alte Refs (ohne neue Key-Spalten) hart abbrechen
        required_cols = {'LOAD_N', 'SIG_0_ABS', 'SIG_0', 'sigma0'}
        missing_cols = required_cols - set(df_existing.columns)
        if missing_cols:
            print(f"!!! [G26-Gen] FATAL: Bestehende Reference_Gauss26.csv hat altes Schema (v18.x).")
            print(f"!!! [G26-Gen]        Fehlende Spalten: {missing_cols}")
            print(f"!!! [G26-Gen]        Aktion: alte Reference loeschen + neu generieren mit v20.0:")
            print(f"!!! [G26-Gen]                  rm {OUTPUT_CSV}")
            print(f"!!! [G26-Gen]                  python 06-POST_ANALYST/error_decomposition/06-generate_gauss26_reference.py")
            sys.exit(1)
        n_existing = len(df_existing)
        # Neue Daten ueberschreiben Zeilen mit gleichem Key (keep='last')
        df_combined = pd.concat([df_existing, df_new], ignore_index=True)
        df_combined = df_combined.drop_duplicates(subset=key_cols, keep='last')
        df_combined = df_combined.sort_values(key_cols).reset_index(drop=True)
        n_added = len(df_combined) - n_existing
        n_updated = len(df_new) - max(0, n_added)
        print(f"  [Upsert] Bestehend: {n_existing} Zeilen | Neu: {len(df_new)} "
              f"-> +{n_added} hinzugefuegt, ~{n_updated} aktualisiert")
        df_ref = df_combined
    else:
        df_ref = df_new

    # Spaltenreihenfolge: Key zuerst (Mesh, Loadcase, Last/Sigma0), dann m + Werte, dann Legacy
    # LOAD_N/SIG_0_ABS/SIG_0/sigma0 sind Teil des disambiguierenden Keys.
    output_cols = ['Mesh', 'Loadcase', 'LOAD_N', 'SIG_0_ABS', 'SIG_0', 'sigma0',
                   'm', 'V_mesh', 'A_mesh', 'Pf_mesh']
    if 'S_mesh' in df_ref.columns:
        output_cols.append('S_mesh')
    output_cols.extend(['VF_Gauss26', 'PF_Gauss26'])
    # Nur vorhandene Spalten selektieren (Schutz nach Upsert aus gemischten Quellen)
    output_cols = [c for c in output_cols if c in df_ref.columns]
    df_ref = df_ref[output_cols]

    df_ref.to_csv(OUTPUT_CSV, sep=';', index=False)

    s_mesh_status = "mit S_mesh" if 'S_mesh' in df_ref.columns else "ohne S_mesh (Legacy)"
    print(f"\n--> Reference gespeichert: {OUTPUT_CSV} ({s_mesh_status})")
    print(f"    {len(df_ref)} Zeilen, {df_ref['Mesh'].nunique()} Meshes, "
          f"{df_ref['Loadcase'].nunique()} Loadcases")


if __name__ == "__main__":
    generate_reference()
