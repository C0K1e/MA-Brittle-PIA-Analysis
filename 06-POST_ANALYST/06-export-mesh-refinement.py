# 06-export-mesh-refinement.py
# Standalone-Entry-Point fuer Mesh-Refinement-Indikator-VTK-Export.
#
# Laeuft NACH einem abgeschlossenen Run und erzeugt im Run-Verzeichnis
# einen neuen vtk_mesh_refinement/-Ordner mit VTK-Dateien (eine pro m-Wert)
# und einer .vtk.series-Datei fuer ParaView.
#
# Output pro Element (3 Skalare + 2 Vektoren = 9 Werte):
#   - Skalare:    E (Element-Fehler), W (Hazard-Gewicht), I = W*E
#   - Vektoren:   r_e_local (Refinement-Faktor pro lokaler Achse)
#                 h_ideal_local (Ziel-Element-Groesse pro lokaler Achse)
#
# WICHTIG: Vector-Komponenten sind LOKALE Element-Achsen, NICHT globale
# x/y/z! Bei strukturierten Meshes mit gekruemmter Topologie (z.B. Lochrand
# bei PWH) sind die Element-Achsen gegenueber dem globalen Frame rotiert.
# Bei achsparallelen Volumen (POR-Innenwuerfel) sind sie identisch.
#
# In ParaView: VECTOR oeffnen, Component "X/Y/Z" auswaehlen — entspricht
# semantisch loc1/loc2/loc3 nach ANSYS-Mesher-Konvention.
#
# Lokale Achsen (HEX8):
#   loc1 = Face (Knoten 0,3,4,7) -> Face (Knoten 1,2,5,6)
#   loc2 = Face (Knoten 0,1,4,5) -> Face (Knoten 2,3,6,7)
#   loc3 = Bottom-Face (0,1,2,3) -> Top-Face (4,5,6,7)
#
# Konfiguration unten im Skript editieren — keine CLI-Argumente.

import os
import sys
import re
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR / "output"))
sys.path.insert(0, str(SCRIPT_DIR))

import vtk_refinement_exporter as vre  # noqa: E402

# =========================================================
# KONFIGURATION (User editiert hier)
# =========================================================

# Liste von Run-Wurzelpfaden (relativ zum Repo-Root oder absolut).
# Skript prozessiert ALLE case-Subordner pro Run.
RUN_DIRS = [
    "05-RUNS/3001-05-13_16-58-PWH-ALL",
    # Weitere Runs:
    # "05-RUNS/2026-05-04_20-33",
]

EPSILON = 0.10         # Globales Fehlerbudget (10 %)

M_VALUES = [1, 5, 10, 15, 20, 25, 30, 40, 50]

# Filter (leer = alle akzeptieren)
FILTER_CASES = ["PWH-400x400x1-RAW-EM-NDX-PIA-V-24x24"]      # z.B. ["POR-50.8x0x1.5-RAW-EM-EMX-PIA-V-8x12x6"]
FILTER_LOC = ["V", "A"]   # beide Domains unterstuetzt (Hex und Quad)

# =========================================================
# Helper: Case-ID-Parsing (vereinfachte Version)
# =========================================================

CASE_ID_PATTERN = re.compile(
    r"^(?P<case>[A-Z0-9]+)"
    r"-(?P<geom>[\d.x]+)"
    r"-(?P<stress>RAW|AVG|ABN|AIN)"
    r"-(?P<int>EM|G\d+)"
    r"-(?P<ref>NDX|EMX|GPX)"
    r"-(?P<crit>PIA|S1)"
    r"-(?P<loc>V|A)"
    r"-(?P<mesh>[\dx]+)$"
)


def parse_case_dir_name(name):
    """Parst case-Dir-Name in Komponenten. None bei Fehler."""
    m = CASE_ID_PATTERN.match(name)
    if not m:
        return None
    return m.groupdict()


def shorten_apdl_name(case_name, info):
    """ANSYS-Job-Name = max 32 Zeichen. Bei POR: nur D_out, dann -<methode>-<mesh>.

    Replikation der Logik aus 01-run_manager.py (_shorten_geom_for_apdl).
    """
    geom = info["geom"]
    case = info["case"]

    if case == "POR":
        # Erste Zahl extrahieren (D_out)
        parts = geom.split("x")
        if parts:
            try:
                d_out = float(parts[0])
                # Integer wenn moeglich (fuer kompakteres Naming)
                if d_out == int(d_out):
                    geom_short = str(int(d_out))
                else:
                    geom_short = str(d_out)
            except ValueError:
                geom_short = geom
        else:
            geom_short = geom
    else:
        geom_short = geom

    return (
        f"{case}-{geom_short}-{info['stress']}-{info['int']}-"
        f"{info['ref']}-{info['mesh']}"
    )


def discover_cases(run_dir):
    """Findet alle case-Subordner in einem Run-Wurzelpfad.

    Returns: list of (case_path, case_name, case_info_dict)
    """
    run_path = Path(run_dir)
    if not run_path.exists():
        print(f"   [Discover] WARN: Run-Verzeichnis nicht gefunden: {run_dir}")
        return []

    cases = []
    for entry in run_path.iterdir():
        if not entry.is_dir():
            continue
        info = parse_case_dir_name(entry.name)
        if info is None:
            continue
        cases.append((entry, entry.name, info))
    return cases


def case_has_required_files(case_path, apdl_name, loc):
    """Prueft ob die noetigen Input-Dateien fuer Refinement-Export da sind.

    Args:
        case_path: Path zum case-Verzeichnis
        apdl_name: gekuerzter ANSYS-Job-Name
        loc:       "V" oder "A"

    Returns: (ok: bool, reason: str)
    """
    tables = case_path / "tables"
    if not tables.exists():
        return False, "tables/ fehlt"

    # Knoten-Geometrie: VTK-Macro-Files ODER Gauss-VEFF-Files
    has_em_nodes = (tables / f"{apdl_name}_VTK_Nodes.out").exists()
    has_gauss_nodes = (tables / f"VEFF_{apdl_name}_Nodes.out").exists()
    if not (has_em_nodes or has_gauss_nodes):
        return False, "weder EM-VTK-Knoten noch Gauss-VEFF-Knoten"

    if loc.upper() == "V":
        # Element-Konnektivitaet
        has_em_elems = (tables / f"{apdl_name}_VTK_Elemente.out").exists()
        has_gauss_elems = (tables / f"VEFF_{apdl_name}_Elemente.out").exists()
        if not (has_em_elems or has_gauss_elems):
            return False, "weder EM-VTK-Elemente noch Gauss-VEFF-Elemente"

        # Volumen
        if not (tables / f"{apdl_name}-V.out").exists():
            return False, f"{apdl_name}-V.out fehlt"

        # Stress (RAW: CSV; AVG: aus VTK_Nodes)
        has_csv = (tables / f"{apdl_name}-S-clean.csv").exists()
        if not has_csv and not has_gauss_nodes:
            return False, "weder -S-clean.csv (RAW) noch Gauss-Stress-Quelle"

        return True, "OK"

    elif loc.upper() == "A":
        # Surface-Face-Daten
        if not (tables / f"{apdl_name}-A.out").exists():
            return False, f"{apdl_name}-A.out fehlt"

        # RAW braucht zusaetzlich -A-S-clean.csv; AVG nicht (S1 aus VTK_Nodes)
        # Hier nur Existenz pruefen — Loader kann mit Fallback umgehen
        return True, "OK"

    else:
        return False, f"unbekannter LOC: {loc}"


# =========================================================
# Main
# =========================================================

def process_case(case_path, case_name, info, epsilon, m_values):
    """Verarbeitet einen einzelnen Case und schreibt VTK-Refinement-Series."""
    apdl_name = shorten_apdl_name(case_name, info)
    tables_dir = str(case_path / "tables")

    ok, reason = case_has_required_files(case_path, apdl_name, info["loc"])
    if not ok:
        print(f"   [Skip] {case_name}: {reason}")
        return False

    out_dir = str(case_path / "vtk_mesh_refinement")
    print(f"\n--> Case: {case_name}")
    print(f"   apdl_name: {apdl_name}")
    print(f"   tables:    {tables_dir}")
    print(f"   output:    {out_dir}")

    try:
        vre.export_refinement_series(
            tables_dir=tables_dir,
            vtk_refinement_dir=out_dir,
            case_name=case_name,
            apdl_name=apdl_name,
            stress_token=info["stress"],
            loc=info["loc"],
            epsilon=epsilon,
            m_values=m_values,
        )
        return True
    except NotImplementedError as exc:
        print(f"   [Skip] {case_name}: {exc}")
        return False
    except Exception as exc:
        print(f"   [Error] {case_name}: {exc}")
        import traceback
        traceback.print_exc()
        return False


def main():
    print("=" * 60)
    print("Mesh Refinement Indicator VTK Export")
    print("=" * 60)
    print(f"EPSILON  = {EPSILON}")
    print(f"M_VALUES = {M_VALUES}")
    print(f"Runs:    {RUN_DIRS}")
    print()

    n_total = 0
    n_ok = 0
    n_skip = 0

    for run_dir in RUN_DIRS:
        print(f"\n{'=' * 60}")
        print(f"Run: {run_dir}")
        print(f"{'=' * 60}")

        cases = discover_cases(run_dir)
        print(f"   {len(cases)} Cases entdeckt")

        for case_path, case_name, info in cases:
            n_total += 1

            # Filter
            if FILTER_CASES and case_name not in FILTER_CASES:
                continue
            if FILTER_LOC and info["loc"] not in FILTER_LOC:
                print(f"   [Filter-LOC] {case_name}: LOC={info['loc']} ausgeschlossen")
                n_skip += 1
                continue

            success = process_case(case_path, case_name, info, EPSILON, M_VALUES)
            if success:
                n_ok += 1
            else:
                n_skip += 1

    print(f"\n{'=' * 60}")
    print(f"Zusammenfassung: {n_ok} OK, {n_skip} skipped, {n_total} total")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
