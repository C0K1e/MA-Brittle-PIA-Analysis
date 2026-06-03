# report_generator.py
# JSON-Metadata Export und Markdown-Report Generator fuer Post-Analyst.
#
# Vollstaendig neu strukturierter Markdown-Report mit 10 Sektionen
# (Kurzueberblick, Lastfall, Geometrie, Methodik, Spannungsreferenzen,
# Analytische Referenz, Numerische Ergebnisse, Fehlerauswertung,
# Visualisierung, Dateien). LOC-aware (Veff vs. Aeff). Automatische
# erklaerende Saetze fuer STRESS/INT/REF/CRIT/LOC. Bild-Auto-Detection
# im plots/-Ordner.
#
# Zwei oeffentliche Funktionen:
#   save_json_metadata(metadata, filepath)
#   generate_markdown_report(metadata, results, image_paths, filepath)

import json
import os
import glob


# ==========================================
# Hilfs-Dictionaries fuer Report-Texte
# ==========================================

CASE_DESCRIPTIONS = {
    "PWH": (
        "Zugversuch an quadratischer Lochplatte (Plate With Hole).\n"
        "- Viertelmodell mit Symmetrie an X=0 und Y=0\n"
        "- Einachsige Zugbelastung (sigma_inf) in X-Richtung\n"
        "- Spannungskonzentrationsfaktor K_t = 3"
    ),
    "PWHR": (
        "Zugversuch an rechteckiger Lochplatte (Plate With Hole, Rectangular).\n"
        "- Viertelmodell mit Symmetrie an X=0 und Y=0\n"
        "- Einachsige Zugbelastung (sigma_inf) in X-Richtung\n"
        "- Spannungskonzentrationsfaktor K_t = 3"
    ),
    "BEAM1": (
        "Reiner Biegebalken mit Biegemoment.\n"
        "- Vollmodell\n"
        "- Biegemoment via SFGRAD (linearer Druckgradient)\n"
        "- Lineare Spannungsverteilung ueber Balkenhoehe"
    ),
    "BEAM2": (
        "Kragarm mit gleichmaessiger Flaechenlast.\n"
        "- Vollmodell\n"
        "- Flaechenlast auf Oberseite (Druck)\n"
        "- Target sigma_max = 1500 MPa"
    ),
    "3PB": (
        "Drei-Punkt-Biegung (3-Point Bending).\n"
        "- Viertelmodell mit Symmetrie an X=0 und Z=0\n"
        "- Einzelkraft in Balkenmitte, Auflager bei L_span/2\n"
        "- 10 mm Ueberhang ueber Auflager"
    ),
    "POR": (
        "Pressure on Ring / Ringdruckversuch.\n"
        "- Viertelmodell mit Symmetrie an X=0 und Y=0\n"
        "- Flache Keramikscheibe unter gleichmaessigem Druck\n"
        "- Auflagerring bei R_sup = 22.73 mm\n"
        "- Material: E = 607000 MPa, Nu = 0.22"
    ),
}

INT_TYPE_LABELS = {
    "EM":  "Elemental Mean (EM)",
    "G1":  "Gauss Integration 1. Ordnung (G1)",
    "G2":  "Gauss Integration 2. Ordnung (G2)",
    "G3":  "Gauss Integration 3. Ordnung (G3)",
    "G4":  "Gauss Integration 4. Ordnung (G4)",
    "G5":  "Gauss Integration 5. Ordnung (G5)",
    "G6":  "Gauss Integration 6. Ordnung (G6)",
    "G7":  "Gauss Integration 7. Ordnung (G7)",
    "G8":  "Gauss Integration 8. Ordnung (G8)",
    "G9":  "Gauss Integration 9. Ordnung (G9)",
    "G26": "Gauss Integration 26. Ordnung (G26, Reference)",
}

# erklaerende Saetze pro Methoden-Baukasten-Token
STRESS_DESCRIPTIONS = {
    "RAW": (
        "Die Auswertung basiert auf ungemittelten Elementspannungen. Dadurch "
        "bleiben lokale Elementwerte erhalten; nodale Glaettung durch ANSYS "
        "wird vermieden."
    ),
    "AVG": (
        "Die Auswertung basiert auf nodal gemittelten Spannungen. Die "
        "Elementwerte werden aus dem geglaetteten FEM-Spannungsfeld abgeleitet."
    ),
}

INT_DESCRIPTIONS_DEFAULT = (
    "Die effektive Groesse wird ueber Gauss-Quadratur berechnet. Das "
    "Spannungsfeld wird innerhalb der Elemente an Gausspunkten ausgewertet "
    "und numerisch integriert. Die Zahl nach G gibt die Quadraturordnung an."
)
INT_DESCRIPTIONS = {
    "EM": (
        "Die effektive Groesse wird ueber Elemental Mean berechnet. Pro "
        "Element beziehungsweise Surface-Face wird ein repraesentativer "
        "Mittelwert verwendet und mit dem zugehoerigen Volumen oder der "
        "Flaeche gewichtet."
    ),
}

REF_DESCRIPTIONS = {
    "EMX": (
        "Die Normierung erfolgt mit dem maximalen Element-Mean-Wert. Diese "
        "Referenz ist besonders passend fuer EM-Auswertungen."
    ),
    "GPX": (
        "Die Normierung erfolgt mit dem maximalen Gausspunktwert der "
        "ausgewerteten Domaene. Diese Referenz ist besonders passend fuer "
        "Gauss-Integrationen."
    ),
    "NDX": (
        "Die Normierung erfolgt mit dem maximalen nodalen Spannungswert aus "
        "der FEM-Auswertung."
    ),
}

CRIT_DESCRIPTIONS = {
    "PIA": (
        "Das Versagenskriterium verwendet den Principal-of-Independent-Action-"
        "Ansatz und beruecksichtigt die positiven Hauptspannungen gemeinsam."
    ),
    "S1": (
        "Das Versagenskriterium verwendet ausschliesslich die erste "
        "Hauptspannung."
    ),
}

LOC_DESCRIPTIONS = {
    "V": "Ausgewertet wird das effektive Volumen Veff.",
    "A": "Ausgewertet wird die effektive Flaeche Aeff.",
}

# Backward-Compat: Legacy-NORM-Labels (intern: NOD/GP/EM)
NORM_LABELS = {
    "NOD": "Nodal (NOD)",
    "GP":  "Gauss-Punkt (GP)",
    "EM":  "Elemental Mean (EM)",
}


# ==========================================
# JSON Metadata Export
# ==========================================

def save_json_metadata(metadata, filepath):
    """Speichert NUR Metadaten als JSON (kein results-Array).

    Args:
        metadata: Dictionary mit Metadaten (von ResultAnalyst._build_metadata())
        filepath: Zielpfad fuer JSON-Datei
    """
    if not metadata:
        return

    try:
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
        print(f"   [Report] JSON gespeichert: {os.path.basename(filepath)}")
    except Exception as e:
        print(f"   [Report] Fehler beim JSON-Speichern: {e}")


# ==========================================
# Markdown Report Generator (10 Sektionen)
# ==========================================

def generate_markdown_report(metadata, results, image_paths, filepath):
    """Erzeugt einen vollstaendigen Markdown-Report mit 10 Sektionen.

    Sektionen:
        0. Kurzueberblick
        1. Lastfall und Auswertungsdomaene
        2. Geometrie und Vernetzung
        3. Spannungsfeld und Auswertungsmethodik
        4. Spannungsreferenzen und Normierung
        5. Analytische Referenz
        6. Numerische Ergebnisse
        7. Fehlerauswertung
        8. Visualisierung
        9. Dateien und Reproduzierbarkeit

    Args:
        metadata: Dictionary mit Metadaten (von ResultAnalyst._build_metadata())
        results: Liste von Result-Dictionaries (m, Veff_ana/Aeff_ana, ...)
        image_paths: Dict mit Standard-Keys 's1', 's2' (relative Pfade);
                     wird durch Auto-Detection im plots/-Ordner ergaenzt.
        filepath: Zielpfad fuer Markdown-Datei
    """
    if not metadata or not results:
        return

    sim = metadata.get("simulation", {}) or {}
    geo = metadata.get("geometry", {}) or {}
    mesh = metadata.get("mesh", {}) or {}
    stress_ref = metadata.get("stress_reference", {}) or {}
    case_id = metadata.get("case_id", "UNKNOWN")
    case_type = metadata.get("case_type", "UNKNOWN")
    timestamp = metadata.get("timestamp", "N/A")

    # Top-Level v20.0+: load_n, sig_0_abs, sig_0, sigma0, pia_fix
    top_load_n   = metadata.get("load_n", metadata.get("load", {}).get("load_n"))
    top_sig0     = metadata.get("sigma0", stress_ref.get("sigma0", 0.0))
    top_sig0_abs = metadata.get("sig_0_abs", sim.get("sig_0_abs", True))
    top_sig0_val = metadata.get("sig_0", sim.get("sig_0_value", 0.0))
    top_pia_fix  = metadata.get("pia_fix", False)

    stress_token = (sim.get("stress") or "").upper()
    int_token    = (sim.get("int_type") or "").upper()
    ref_token    = (sim.get("ref") or "").upper()
    crit_token   = (sim.get("crit") or "").upper()
    loc_token    = (sim.get("loc") or "V").upper()

    # LOC-aware Effective-Quantity Keys
    if loc_token == "A":
        eff_ana_key, eff_num_key = "Aeff_ana", "Aeff_num"
        err_rel_key, err_fold_key = "Err_A_rel", "Err_A_Fold"
        eff_label = "Aeff"
    else:
        eff_ana_key, eff_num_key = "Veff_ana", "Veff_num"
        err_rel_key, err_fold_key = "Err_V_rel", "Err_V_Fold"
        eff_label = "Veff"

    lines = []

    # === Header ===
    lines.append(f"# Report: {case_id}")
    lines.append("")

    # === 0. Kurzueberblick ===
    lines.append("## 0. Kurzüberblick")
    lines.append("")
    lines.append("| Feld | Wert |")
    lines.append("|---|---|")
    lines.append(f"| Erstellungsdatum | {timestamp} |")
    lines.append(f"| Lastfall (case_type) | {case_type} |")
    lines.append(f"| Auswertungsdomäne (LOC) | {loc_token} ({eff_label}) |")
    lines.append(f"| Spannungsfeld (STRESS) | {stress_token or 'N/A'} |")
    lines.append(f"| Integration (INT) | {INT_TYPE_LABELS.get(int_token, int_token or 'N/A')} |")
    lines.append(f"| Referenzierung (REF) | {ref_token or 'N/A'} |")
    lines.append(f"| Bruchkriterium (CRIT) | {crit_token or 'N/A'} |")
    mesh_str_short = _build_mesh_str(mesh)
    lines.append(f"| Mesh | {mesh_str_short} |")
    lines.append(f"| Last (load_n) | {_fmt_float(top_load_n)} |")
    lines.append(f"| sigma0 | {_fmt_float(top_sig0)} |")
    lines.append(f"| PIA_FIX | {'1 (PIA-at-GP)' if top_pia_fix else '0 (Legacy)'} |")
    lines.append("")

    # === 1. Lastfall und Auswertungsdomaene ===
    lines.append("## 1. Lastfall und Auswertungsdomäne")
    lines.append("")
    lines.append(CASE_DESCRIPTIONS.get(case_type, f"Unbekannter Lastfall: `{case_type}`"))
    lines.append("")
    lines.append(f"**Auswertung:** {LOC_DESCRIPTIONS.get(loc_token, f'LOC={loc_token}')}")
    lines.append("")

    # === 2. Geometrie und Vernetzung ===
    lines.append("## 2. Geometrie und Vernetzung")
    lines.append("")
    lines.append("| Parameter | Wert |")
    lines.append("|---|---|")
    lines.append(f"| CASE_X / CASE_Y / CASE_Z | {geo.get('case_x', 'N/A')} / {geo.get('case_y', 'N/A')} / {geo.get('case_z', 'N/A')} |")
    lines.append(f"| Mesh (NX × NY × NZ) | {mesh_str_short} |")
    lines.append(f"| Elementanzahl | {mesh.get('num_elements', 'N/A')} |")
    lines.append(f"| Elementtyp | SOLID186 (20-Knoten Hex) |")
    lines.append(f"| Symmetriefaktor | {geo.get('symmetry_factor', 'N/A')} |")
    lines.append("")

    # Vtot/Atot Vergleichs-Tabelle
    vtot_ana = geo.get("vtot_exact")
    vtot_num = geo.get("vtot_num")
    atot_ana = geo.get("atot_exact")
    atot_num = geo.get("atot_num")
    have_atot = (atot_ana is not None) or (atot_num not in (None, 0, 0.0))
    lines.append("| Größe | Analytisch | Numerisch |")
    lines.append("|---|---:|---:|")
    if vtot_ana is not None or vtot_num is not None:
        lines.append(f"| Gesamtvolumen Vtot | {_fmt_float(vtot_ana)} | {_fmt_float(vtot_num)} |")
    if have_atot:
        lines.append(f"| Gesamtfläche Atot | {_fmt_float(atot_ana)} | {_fmt_float(atot_num)} |")
    lines.append("")

    # === 3. Spannungsfeld und Auswertungsmethodik ===
    lines.append("## 3. Spannungsfeld und Auswertungsmethodik")
    lines.append("")
    lines.append("| Token | Wert | Bedeutung |")
    lines.append("|---|---|---|")
    lines.append(f"| STRESS | {stress_token or 'N/A'} | Spannungsdarstellung |")
    lines.append(f"| INT | {int_token or 'N/A'} | Integrationsverfahren |")
    lines.append(f"| REF | {ref_token or 'N/A'} | Referenz/Normierungswert |")
    lines.append(f"| CRIT | {crit_token or 'N/A'} | Bruchkriterium |")
    lines.append(f"| LOC | {loc_token or 'N/A'} | Auswertungsdomäne |")
    lines.append("")

    # Automatische erklaerende Saetze
    if stress_token in STRESS_DESCRIPTIONS:
        lines.append(f"**STRESS={stress_token}**: {STRESS_DESCRIPTIONS[stress_token]}")
        lines.append("")
    if int_token:
        if int_token in INT_DESCRIPTIONS:
            lines.append(f"**INT={int_token}**: {INT_DESCRIPTIONS[int_token]}")
        elif int_token.startswith("G"):
            lines.append(f"**INT={int_token}**: {INT_DESCRIPTIONS_DEFAULT}")
        lines.append("")
    if ref_token in REF_DESCRIPTIONS:
        lines.append(f"**REF={ref_token}**: {REF_DESCRIPTIONS[ref_token]}")
        lines.append("")
    if crit_token in CRIT_DESCRIPTIONS:
        lines.append(f"**CRIT={crit_token}**: {CRIT_DESCRIPTIONS[crit_token]}")
        lines.append("")
    if loc_token in LOC_DESCRIPTIONS:
        lines.append(f"**LOC={loc_token}**: {LOC_DESCRIPTIONS[loc_token]}")
        lines.append("")

    # PIA_FIX-Hinweis
    if top_pia_fix:
        lines.append(
            "**PIA_FIX = 1**: PIA-Auswertung am Gauss-Punkt (lokale Auswertung, "
            "tensile cutoff am GP). Mathematisch konsistent mit Spec ab v19.0."
        )
    else:
        lines.append(
            "**PIA_FIX = 0**: Legacy-Pfad (PIA-Aggregation an Knoten + Interpolation). "
            "Bei multiaxialen Cases methodisch nicht aequivalent zu PIA-at-GP."
        )
    lines.append("")

    # === 4. Spannungsreferenzen und Normierung ===
    lines.append("## 4. Spannungsreferenzen und Normierung")
    lines.append("")
    smax_ref   = stress_ref.get("smax_ref", 0.0)
    smax_nodal = stress_ref.get("smax_nodal", 0.0)
    smax_norm  = stress_ref.get("smax_norm", 0.0)
    sigma0_val = stress_ref.get("sigma0", top_sig0)

    # Source-Hint pro REF-Token
    norm_source = {
        "EMX": "max(Element-Mean S1) ueber alle Elemente",
        "GPX": "max(S1) ueber alle Gauss-Punkte (SMAX_GAUSS_VOL/SURF)",
        "NDX": "max(S1) am FEM-Knoten (Header-Wert)",
    }.get(ref_token, "abhängig von REF")

    sig0_mode = "Absolut (TRUE)" if top_sig0_abs else "Relativ (FALSE → sig_0 × smax_nodal)"

    lines.append("| Wert | Bedeutung | Quelle | Betrag |")
    lines.append("|---|---|---|---:|")
    lines.append(f"| smax_ref / S_ana | analytische Maximalspannung | analytische Formel / LOADCASE_REGISTRY | {_fmt_float(smax_ref)} |")
    lines.append(f"| smax_nodal | FEM-Maximalspannung | Header oder RAW-PRESOL | {_fmt_float(smax_nodal)} |")
    lines.append(f"| smax_norm / S_num | Normierungswert (REF={ref_token}) | {norm_source} | {_fmt_float(smax_norm)} |")
    lines.append(f"| sigma0 | Weibull-Skalenparameter | CSV ({sig0_mode}) | {_fmt_float(sigma0_val)} |")
    lines.append("")

    # === 5. Analytische Referenz ===
    lines.append("## 5. Analytische Referenz")
    lines.append("")
    lines.append(_build_analytical_description(case_type, loc_token, crit_token))
    lines.append("")

    # === 6. Numerische Ergebnisse ===
    lines.append("## 6. Numerische Ergebnisse")
    lines.append("")

    showcase_m = [1, 5, 10, 15, 20, 25, 30, 40, 50]
    result_by_m = {r.get("m"): r for r in results if r.get("m") is not None}

    lines.append(
        f"| m | {eff_label}_ana | {eff_label}_num | rel. Fehler [%] | Fehlerfaktor (Fold) | Pf_ana | Pf_num | Pf-Fehler [%] |"
    )
    lines.append("|--:|--------:|--------:|--------------:|-----------:|-------:|-------:|---------------:|")

    for m_val in showcase_m:
        r = result_by_m.get(m_val)
        if not r:
            continue
        eff_ana = r.get(eff_ana_key)
        eff_num = r.get(eff_num_key)
        err_rel = r.get(err_rel_key)
        err_fold = r.get(err_fold_key)
        pf_ana = r.get("Pf_ana")
        pf_num = r.get("Pf_num")
        err_pf_rel = r.get("Err_Pf_rel")

        lines.append(
            f"| {m_val} "
            f"| {_fmt_float(eff_ana)} "
            f"| {_fmt_float(eff_num)} "
            f"| {_fmt_pct(err_rel)} "
            f"| {_fmt_float(err_fold)} "
            f"| {_fmt_float(pf_ana)} "
            f"| {_fmt_float(pf_num)} "
            f"| {_fmt_pct(err_pf_rel)} |"
        )
    lines.append("")
    lines.append(f"*Vollständige Daten: `ANALYSIS_{case_id}.csv`*")
    lines.append("")

    # === 7. Fehlerauswertung ===
    lines.append("## 7. Fehlerauswertung")
    lines.append("")
    extended_csv = _find_extended_csv(filepath, case_id)
    if extended_csv:
        lines.append(f"**Erweiterte Fehlerzerlegung verfügbar**: `{os.path.basename(extended_csv)}`")
        lines.append("")
        lines.append(
            "Die v15-Decomposition zerlegt den Gesamt-Faktor `rho_Total` multiplikativ in:"
        )
        lines.append(
            "- `rho_Disk` — Diskretisierungs-Fehler (FE-Netz)"
        )
        lines.append(
            "- `rho_Avg` — Glättung durch nodale Mittelung (= 1 bei STRESS=RAW)"
        )
        lines.append(
            "- `rho_Int` — Integrationsmethode (= 1 bei INT=G26-Reference)"
        )
        lines.append(
            "- `rho_Ref` — Referenz-/Normierungs-Effekt (= 1 für Pf, weil smax kürzt)"
        )
        lines.append("")
        lines.append(
            "Identität: `rho_Total = rho_Disk · rho_Avg · rho_Int · rho_Ref` "
            "(multiplikativ, alle Faktoren ≈ 1 bei perfekter Konvergenz)."
        )
        lines.append("")

        # rho-Decomposition-Tabellen pro m
        ext_rows_by_m = _read_extended_csv_rho(extended_csv, eff_label, showcase_m)
        if ext_rows_by_m:
            # Tabelle 1: Veff/Aeff-Decomposition
            lines.append(
                f"**Tabelle 7.1 — Auf {eff_label} propagierte Fehleranteile "
                f"(rho-Faktoren, inklusive Referenzierung)**"
            )
            lines.append("")
            lines.append(
                f"| m | rho_Disk | rho_Avg | rho_Int | rho_Ref | rho_Total |"
            )
            lines.append("|--:|--------:|--------:|--------:|--------:|---------:|")
            for m_val in showcase_m:
                row = ext_rows_by_m.get(m_val)
                if not row:
                    continue
                lines.append(
                    f"| {m_val} "
                    f"| {_fmt_rho(row.get(f'{eff_label}_rho_Disk'))} "
                    f"| {_fmt_rho(row.get(f'{eff_label}_rho_Avg'))} "
                    f"| {_fmt_rho(row.get(f'{eff_label}_rho_Int'))} "
                    f"| {_fmt_rho(row.get(f'{eff_label}_rho_Ref'))} "
                    f"| {_fmt_rho(row.get(f'{eff_label}_rho_Total'))} |"
                )
            lines.append("")

            # Tabelle 2: Pf-Decomposition (LOC-unabhaengig)
            lines.append(
                "**Tabelle 7.2 — Auf Pf propagierte Fehleranteile "
                "(rho-Faktoren, Referenzfaktor exakt 1)**"
            )
            lines.append("")
            lines.append(
                "| m | rho_Disk | rho_Avg | rho_Int | rho_Ref | rho_Total |"
            )
            lines.append("|--:|--------:|--------:|--------:|--------:|---------:|")
            for m_val in showcase_m:
                row = ext_rows_by_m.get(m_val)
                if not row:
                    continue
                lines.append(
                    f"| {m_val} "
                    f"| {_fmt_rho(row.get('Pf_rho_Disk'))} "
                    f"| {_fmt_rho(row.get('Pf_rho_Avg'))} "
                    f"| {_fmt_rho(row.get('Pf_rho_Int'))} "
                    f"| {_fmt_rho(row.get('Pf_rho_Ref'))} "
                    f"| {_fmt_rho(row.get('Pf_rho_Total'))} |"
                )
            lines.append("")
            lines.append(
                "*Lese-Hilfe: `rho = 1.000` → exakte Übereinstimmung. "
                "`rho > 1` → numerischer Wert überschätzt, `rho < 1` → unterschätzt. "
                "Pro Zeile gilt: `rho_Total = rho_Disk · rho_Avg · rho_Int · rho_Ref`.*"
            )
            lines.append("")
            lines.append(
                f"*Hinweis zur Propagation: Beide Tabellen beschreiben **dieselbe "
                f"Fehlerursache** (Diskretisierung + Quadratur des Spannungsfelds), "
                f"aber projiziert auf zwei unterschiedliche Zielgrößen. Die Werte "
                f"unterscheiden sich, weil {eff_label} den Normierungs-Hebel "
                f"`(σ_ana/σ_num)^m` über `rho_Ref` mitführt (kann kompensieren oder "
                f"verstärken), während Pf eine nichtlineare Weibull-Sättigungs-"
                f"Abbildung `1 - exp(-H)` anwendet und `rho_Ref` algebraisch zu 1 "
                f"kürzt.*"
            )
            lines.append("")
        lines.append("Detail-Plots in `../plots/` (sofern v15-Plot-Skripte ausgeführt).")
    else:
        lines.append("Keine erweiterte Fehlerzerlegung vorhanden.")
        lines.append("")
        lines.append(
            "Hinweis: für die v15-Decomposition `ERR_EXT=1` in `00-cases.csv` setzen + "
            "`Reference_Gauss26.csv` muss vorliegen (siehe v20.0/v20.2-Workflow). Bei "
            "PIA_FIX=0 ist die Decomposition seit v20.2 blockiert."
        )
    lines.append("")

    # === 8. Visualisierung ===
    lines.append("## 8. Visualisierung")
    lines.append("")
    image_blocks = _detect_images(filepath, case_id, image_paths)
    if image_blocks:
        for label, rel_path in image_blocks:
            lines.append(f"### {label}")
            lines.append("")
            lines.append(f"![{label}]({rel_path})")
            lines.append("")
    else:
        lines.append("Keine Visualisierungen vorhanden im plots/-Ordner.")
        lines.append("")

    # === 9. Dateien und Reproduzierbarkeit ===
    lines.append("## 9. Dateien und Reproduzierbarkeit")
    lines.append("")
    lines.append(f"- ANALYSIS-CSV: `tables/ANALYSIS_{case_id}.csv`")
    if extended_csv:
        lines.append(f"- Extended-CSV: `tables/{os.path.basename(extended_csv)}`")
    lines.append(f"- JSON-Metadata: `tables/{case_id}.json`")
    lines.append(f"- Markdown-Report (diese Datei): `tables/{case_id}_REPORT.md`")
    lines.append(f"- Plots: `plots/`")
    lines.append(f"- ANSYS-Rohdaten: `_ANSYS_Dateien/`")
    lines.append("")

    # === Footer ===
    lines.append("---")
    lines.append("*Generiert von report_generator.py v20.3*")

    try:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write("\n".join(lines))
        print(f"   [Report] Report gespeichert: {os.path.basename(filepath)}")
    except Exception as e:
        print(f"   [Report] Fehler beim Report-Speichern: {e}")


# ==========================================
# Hilfsfunktionen
# ==========================================

def _fmt_float(val):
    """Formatiert Zahlenwerte: wissenschaftlich fuer sehr kleine/grosse, sonst 6 signifikant.

    Robust gegen None und nicht-numerische Werte (gibt 'N/A' zurueck).
    """
    if val is None:
        return "N/A"
    try:
        v = float(val)
        if v == 0:
            return "0"
        if abs(v) < 0.001 or abs(v) > 1e6:
            return f"{v:.6e}"
        return f"{v:.6g}"
    except (ValueError, TypeError):
        return str(val)


def _fmt_pct(val):
    """Formatiert relativen Fehler als Prozent (zwei Nachkommastellen).

    Args:
        val: Fehler-Faktor (z.B. 0.01 = 1%) — wird mit 100 multipliziert.
    """
    if val is None:
        return "N/A"
    try:
        return f"{float(val) * 100:.2f}"
    except (ValueError, TypeError):
        return "N/A"


def _build_mesh_str(mesh):
    """Baut Mesh-String '8x12x6' oder '8x12' aus mesh-Dict."""
    nx = mesh.get("mesh_x", 0)
    ny = mesh.get("mesh_y", 0)
    nz = mesh.get("mesh_z", 1)
    if not nx or not ny:
        return "N/A"
    s = f"{nx}x{ny}"
    try:
        if int(nz) > 1:
            s += f"x{nz}"
    except (ValueError, TypeError):
        pass
    return s


def _build_analytical_description(case_type, loc, crit):
    """Erzeugt Kurzbeschreibung der analytischen Loesung pro Lastfall × LOC × CRIT."""
    if case_type == "POR":
        if loc == "V" and crit == "PIA":
            return (
                "POR + LOC=V + CRIT=PIA: `Veff = Aeff_full · t / (2(m+1))` mit voller "
                "Plattenformel `_compute_aeff_por` (sigma_rr + sigma_phiphi-Integral)."
            )
        if loc == "V" and crit == "S1":
            return (
                "POR + LOC=V + CRIT=S1: **blockiert** (`NotImplementedError`). "
                "Durch-die-Dicke-Integration mit Membran-Term braucht saubere Ableitung. "
                "Workaround: LOC=A nutzen oder CRIT=PIA."
            )
        if loc == "A" and crit == "PIA":
            return (
                "POR + LOC=A + CRIT=PIA: `_compute_aeff_por` (volle Formel über "
                "sigma_rr^m + sigma_phiphi^m)."
            )
        if loc == "A" and crit == "S1":
            return (
                "POR + LOC=A + CRIT=S1: `_compute_aeff_por_s1` (v20.0 NEU). "
                "Geschlossene Form: `Aeff_S1 = pi·R_s² / α · (1 − (1-α)^(m+1)) / (m+1)` "
                "mit α = (1+3ν)/((1-ν)·(R_s/R_d)² + 2(1+ν)). "
                "Plattenbiegungs-Approximation (Membran-Term vernachlässigt)."
            )
    if case_type in ("BEAM1", "3PB"):
        return (
            f"{case_type}: Veff = V_total / (2·(m+1)) (uniaxial, ν=0 → S1=PIA identisch). "
            "Aeff = A_total für BEAM1 (uniform σ); für 3PB nicht analytisch (Auflager-Ausschluss)."
        )
    if case_type == "BEAM2":
        return (
            "BEAM2 (Kragarm): Veff = V_total / (2(m+1)·(2m+1)). "
            "Aeff = A_total / (2m+1) (quadratisches Moment)."
        )
    if case_type in ("PWH", "PWHR"):
        return (
            f"{case_type}: Veff via Tabelle in `analytical_reference_data.py` "
            f"(diskrete CASE_X-Keys). Aeff aktuell blockiert (biaxiales Lochrand-Feld). "
            f"Siehe LOADCASE_REGISTRY in `06-POST_ANALYST/analytical/analytical_helper.py` für vollständigen Status."
        )
    return f"Lastfall {case_type}: siehe `06-POST_ANALYST/analytical/analytical_helper.py` LOADCASE_REGISTRY."


def _find_extended_csv(report_filepath, case_id):
    """Sucht ANALYSIS_{case_id}_extended.csv im selben tables-Ordner.

    Returns:
        Absoluter Pfad falls vorhanden, sonst None.
    """
    tables_dir = os.path.dirname(report_filepath)
    candidate = os.path.join(tables_dir, f"ANALYSIS_{case_id}_extended.csv")
    return candidate if os.path.exists(candidate) else None


def _fmt_rho(val):
    """v20.5: Formatiert rho-Faktor (typisch 0.9..1.1) mit 4 Nachkommastellen.

    Robust gegen None/leer/nan: liefert 'N/A' bei nicht-numerischen Werten.
    """
    if val is None or val == "":
        return "N/A"
    try:
        v = float(val)
        if v != v:  # NaN-Check
            return "N/A"
        return f"{v:.4f}"
    except (ValueError, TypeError):
        return "N/A"


def _read_extended_csv_rho(extended_csv_path, eff_label, target_m_values):
    """v20.5: Liest die rho-Spalten aus _extended.csv fuer die Showcase-m-Werte.

    Args:
        extended_csv_path: Absoluter Pfad zu ANALYSIS_*_extended.csv
        eff_label: 'Veff' oder 'Aeff' (LOC-aware)
        target_m_values: Liste der m-Werte fuer die wir Zeilen brauchen

    Returns:
        Dict {m: row_dict}. Nur Zeilen mit m in target_m_values. Robust gegen
        Lese-Fehler — gibt {} zurueck bei Problem.
    """
    import csv as _csv
    target_set = set(target_m_values)
    rho_keys = []
    for q in (eff_label, "Pf"):
        for x in ("Disk", "Avg", "Int", "Ref", "Total"):
            rho_keys.append(f"{q}_rho_{x}")

    rows_by_m = {}
    try:
        with open(extended_csv_path, 'r', encoding='utf-8') as f:
            reader = _csv.DictReader(f, delimiter=';')
            for raw_row in reader:
                try:
                    m_val = int(float(raw_row.get("m", "")))
                except (ValueError, TypeError):
                    continue
                if m_val not in target_set:
                    continue
                # Nur die rho-Spalten extrahieren (kompakter Dict)
                row_filtered = {k: raw_row.get(k) for k in rho_keys if k in raw_row}
                rows_by_m[m_val] = row_filtered
    except Exception as e:
        print(f"   [Report] WARN: extended.csv-Read fehlgeschlagen ({e}); rho-Tabelle wird uebersprungen.")
        return {}
    return rows_by_m


def _detect_images(report_filepath, case_id, image_paths):
    """Auto-Detection der Visualisierungs-Bilder im plots/-Ordner.

    Sucht (in Reihenfolge):
        FEM:        *_s1.png, *_s2.png  (über image_paths-Defaults)
        Analytik:   *_analytical_s1.png, *_analytical_s2.png
        Pfad:       *_path_stresses.png, *_pfad_spannungen.png
        Fehler:     *_error_decomposition_*.png

    Returns:
        Liste von (label, relativer_pfad)-Tupeln. Pfade sind relativ zum Report
        (typisch '../plots/...').
    """
    blocks = []

    # tables_dir / ../plots/
    tables_dir = os.path.dirname(report_filepath)
    case_dir = os.path.dirname(tables_dir)
    plots_dir = os.path.join(case_dir, "plots")

    # 1) FEM S1/S2 aus image_paths-Standard (wenn existieren)
    if image_paths:
        for key, label in [("s1", "FEM Erste Hauptspannung (S1)"),
                           ("s2", "FEM Zweite Hauptspannung (S2)")]:
            rel_path = image_paths.get(key, "")
            if rel_path:
                abs_path = os.path.join(tables_dir, rel_path)
                if os.path.exists(abs_path):
                    blocks.append((label, rel_path))

    # 2) Auto-Scan im plots/-Ordner
    if not os.path.isdir(plots_dir):
        return blocks

    # Pattern -> Label
    auto_patterns = [
        ("*_analytical_s1.png", "Analytische Erste Hauptspannung (S1, Analytik)"),
        ("*_analytical_s2.png", "Analytische Zweite Hauptspannung (S2, Analytik)"),
        ("*_path_stresses.png", "Pfadspannungen entlang Kerbrand"),
        ("*_pfad_spannungen.png", "Pfadspannungen entlang Kerbrand"),
        ("*_error_decomposition_veff.png", "Fehlerzerlegung Veff (5-Klassen)"),
        ("*_error_decomposition_aeff.png", "Fehlerzerlegung Aeff (5-Klassen)"),
        ("*_error_decomposition_pf.png", "Fehlerzerlegung Pf (5-Klassen)"),
    ]

    seen_paths = {p for (_, p) in blocks}  # gegen Duplikate
    for pattern, label in auto_patterns:
        matches = sorted(glob.glob(os.path.join(plots_dir, pattern)))
        for abs_path in matches:
            rel_path = os.path.join("..", "plots", os.path.basename(abs_path)).replace("\\", "/")
            if rel_path in seen_paths:
                continue
            seen_paths.add(rel_path)
            blocks.append((label, rel_path))

    return blocks