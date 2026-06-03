import os
import shutil
import subprocess
import csv
import sys
from datetime import datetime
import importlib.util

# Zentrales Naming-Modul (STRESS/INT/REF Tokens + Validierung + APDL-Bruecke)
# method_naming.py liegt jetzt in 06-POST_ANALYST/helpers/
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "06-POST_ANALYST", "helpers"))
import method_naming as mn

# === KONFIGURATION ===

# Maschinen-spezifische Pfade (ANSYS-EXE, Kernzahl) liegen in config_local.py
# (gitignored, pro Rechner). Vorlage: config_local.example.py -> kopieren zu
# config_local.py und eigene Pfade eintragen. Haelt private/lokale Pfade aus
# dem oeffentlichen Repo heraus.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import config_local as cfg
except ModuleNotFoundError:
    print("!!! [Manager] FATAL: config_local.py nicht gefunden.")
    print("!!!         Aktion: config_local.example.py -> config_local.py kopieren")
    print("!!!                 und die ANSYS-Pfade fuer diesen Rechner eintragen.")
    sys.exit(1)

KONFIGURATION_SANSYS2 = True

if KONFIGURATION_SANSYS2 is True:
    ANSYS_EXE_PATH = cfg.ANSYS_EXE_PATH_SERVER
    NUM_CORES = cfg.NUM_CORES_SERVER
    WORKING_DIR = os.getcwd()
    RUNS_DIR = os.path.join(WORKING_DIR, "05-RUNS")
    # CSV-Templates wurden nach csv-templates/ verschoben (Repo-Reorganisation)
    CSV_FILE = os.path.join("csv-templates", "00-cases-PWH-400-ALL-2.csv")
    #CSV_FILE = "00-cases.csv"
else:
    ANSYS_EXE_PATH = cfg.ANSYS_EXE_PATH_LAPTOP
    NUM_CORES = cfg.NUM_CORES_LAPTOP
    WORKING_DIR = os.getcwd()
    RUNS_DIR = os.path.join(WORKING_DIR, "05-RUNS")
    CSV_FILE = "00-cases.csv"


# Schalter für automatisches Post-Processing
RUN_POST_PROCESSING = True

# Post-Analyst Modul (wird später geladen)
ANALYST_SCRIPT = os.path.join("06-POST_ANALYST", "06-post_analyst.py")
post_module = None

if RUN_POST_PROCESSING and os.path.exists(ANALYST_SCRIPT):
    try:
        spec = importlib.util.spec_from_file_location("PostAnalyst", ANALYST_SCRIPT)
        post_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(post_module)
        print(f"   [Manager] Post-Analyst Modul '{ANALYST_SCRIPT}' erfolgreich geladen.")
    except Exception as e:
        print(f"!!! [Manager] WARN: Konnte Post-Analyst nicht laden: {e}")
        post_module = None
        RUN_POST_PROCESSING = False

# === QUELLDATEIEN DEFINIEREN ===
SOURCE_FILES = [
    os.path.join("02-MAIN", "02-main.inp"),               # Liegt in 02-MAIN/
    os.path.join("03-GEOMETRY", "geom_BEAM1.inp"),     # Biegebalken (Moment via SFGRAD)
    os.path.join("03-GEOMETRY", "geom_PWH.inp"),       # Plate With Hole (Lochplatte)
    os.path.join("03-GEOMETRY", "geom_POR.inp"),       # Pressure on Ring (Ringdruckversuch)

    os.path.join("04-POST", "post_GAUSS_unified.inp"),      # NOD + GaussNorm konsolidiert
    os.path.join("04-POST", "post_EM_NOD.inp"),
    os.path.join("04-POST", "post_EM_RAW.inp"),
    os.path.join("04-POST", "post_EM_NOD_surf.inp"),     # Surface-EM-NOD (LOC=A)
    os.path.join("04-POST", "post_EM_RAW_surf.inp"),     # Surface-EM-RAW (LOC=A, AVG=RAW)
    os.path.join("04-POST", "x_Effektives_Volumen_NOD.mac"),  # NOD-Variante (umbenannt fuer Symmetrie zu _RAW.mac)

    # === Unified Fortran (aktive EXE fuer alle Gauss-Pipelines) ===
    # Eine EXE fuer alle vier Pipelines (NOD/RAW × NOD/GP-Norm) + Gauss 1-9, 26.
    # Steuerung ueber AVGModeRAW + GaussPNorm in effVol_Parameter.out.
    os.path.join("04-POST", "FORTRAN", "Effektives_Volumen_Unified.exe"),
    os.path.join("04-POST", "FORTRAN", "effektivesVol_unified.f90"),
    os.path.join("04-POST", "FORTRAN", "overhead_unified.f90"),

    # === PIAFix-Variante (PIA-at-Gauss-Point Korrektur) ===
    # Parallele EXE neben Legacy. Wird verwendet wenn PIA_FIX=1 in 00-cases.csv.
    # Optional: Kopie nur falls Datei existiert (kein FATAL, weil v19.0-Patch).
    os.path.join("04-POST", "FORTRAN", "Effektives_Volumen_Unified_PIAFix.exe"),
    os.path.join("04-POST", "FORTRAN", "effektivesVol_unified_piafix.f90"),
    os.path.join("04-POST", "FORTRAN", "overhead_unified_piafix.f90"),

    # APDL RAWNodes
    os.path.join("04-POST", "post_GAUSS_RAW.inp"),
    os.path.join("04-POST", "x_Effektives_Volumen_RAW.mac"),

    # VTK Geometrie-Export (EM-Cases)
    os.path.join("04-POST", "x_Export_VTK_Geometry.mac"),
]

def _shorten_geom_for_apdl(case_id, geom_str):
    """Verkuerzt geom_str fuer APDL-Jobname (32-Zeichen-Limit).

    POR: Nur Durchmesser (case_x), da case_y=0 und case_z kurz.
    Andere Cases: Unveraendert (passen in 32 Zeichen).
    """
    if case_id == "POR":
        return geom_str.split("x")[0]  # "50.8x0x1.5" -> "50.8"
    return geom_str


def _cleanup_case_dir(case_dir):
    """Verschiebt ANSYS-Ausgabedateien in _ANSYS_Dateien/ Unterordner.

    Behaelt current_params.inp im Case-Root.
    Ueberspringt Unterordner (plots/, tables/, _ANSYS_Dateien/).
    """
    ansys_dir = os.path.join(case_dir, "_ANSYS_Dateien")
    os.makedirs(ansys_dir, exist_ok=True)

    keep_at_root = {"current_params.inp"}

    for entry in os.listdir(case_dir):
        entry_path = os.path.join(case_dir, entry)
        if not os.path.isfile(entry_path):
            continue
        if entry.lower() in keep_at_root:
            continue
        shutil.move(entry_path, os.path.join(ansys_dir, entry))


def _cleanup_batch_dir(batch_dir):
    """Verschiebt Berechnungscodes in _Berechnungscodes/ Unterordner.

    Verschiebt alle SOURCE_FILES-Basenames.
    Laesst NAME_MAPPING.csv, PLOT_*.png und Case-Unterordner am Root.
    """
    code_dir = os.path.join(batch_dir, "_Berechnungscodes")
    os.makedirs(code_dir, exist_ok=True)

    source_basenames = {os.path.basename(f) for f in SOURCE_FILES}

    moved = 0
    for entry in os.listdir(batch_dir):
        if entry in source_basenames:
            src = os.path.join(batch_dir, entry)
            if os.path.isfile(src):
                shutil.move(src, os.path.join(code_dir, entry))
                moved += 1

    print(f"   [Cleanup] {moved} Dateien nach _Berechnungscodes/ verschoben.")


def run_study():
    # Pre-Run-Validierung — kritische Resourcen pruefen
    # bevor ANSYS startet. Verhindert silent-fails durch fehlende Unified-EXE.
    UNIFIED_EXE = os.path.join("04-POST", "FORTRAN", "Effektives_Volumen_Unified.exe")
    if not os.path.exists(UNIFIED_EXE):
        print(f"!!! [Manager] FATAL: Unified-Fortran-EXE fehlt: {UNIFIED_EXE}")
        print(f"!!! [Manager]        Aktion: Re-Compile noetig.")
        print(f"!!! [Manager]        Compile (in 04-POST/FORTRAN/): ifx /O2 overhead_unified.f90 effektivesVol_unified.f90 /Fe:Effektives_Volumen_Unified.exe")
        sys.exit(1)
    if not os.path.exists(CSV_FILE):
        print(f"!!! [Manager] FATAL: CSV-Datei fehlt: {CSV_FILE}")
        print(f"!!! [Manager]        KONFIGURATION_SANSYS2 = {KONFIGURATION_SANSYS2}")
        print(f"!!! [Manager]        Aktion: CSV erzeugen oder KONFIGURATION_SANSYS2 anpassen.")
        sys.exit(1)
    if not os.path.exists(ANSYS_EXE_PATH):
        print(f"!!! [Manager] FATAL: ANSYS-EXE nicht gefunden: {ANSYS_EXE_PATH}")
        print(f"!!! [Manager]        KONFIGURATION_SANSYS2 = {KONFIGURATION_SANSYS2}")
        print(f"!!! [Manager]        Aktion: ANSYS_EXE_PATH oder KONFIGURATION_SANSYS2 anpassen.")
        sys.exit(1)

    # PIAFix-EXE-Check, falls irgendein Case PIA_FIX=1 setzt
    UNIFIED_PIAFIX_EXE = os.path.join("04-POST", "FORTRAN", "Effektives_Volumen_Unified_PIAFix.exe")
    needs_piafix = False
    try:
        with open(CSV_FILE, 'r', encoding='utf-8') as _f:
            _reader = csv.DictReader(_f, delimiter=';')
            for _row in _reader:
                if _row.get('PIA_FIX', '0').strip() == '1':
                    needs_piafix = True
                    break
    except Exception:
        pass  # CSV-Parser-Fehler werden weiter unten beim Haupt-Loop behandelt
    if needs_piafix and not os.path.exists(UNIFIED_PIAFIX_EXE):
        print(f"!!! [Manager] FATAL: PIAFix-Fortran-EXE fehlt: {UNIFIED_PIAFIX_EXE}")
        print(f"!!! [Manager]        Mindestens ein Case in {CSV_FILE} hat PIA_FIX=1 gesetzt.")
        print(f"!!! [Manager]        Aktion: Re-Compile noetig.")
        print(f"!!! [Manager]        Quick-Compile (ifx):")
        print(f"!!! [Manager]          ifx /O2 04-POST/FORTRAN/overhead_unified_piafix.f90 \\")
        print(f"!!! [Manager]                  04-POST/FORTRAN/effektivesVol_unified_piafix.f90 \\")
        print(f"!!! [Manager]                  /Fe:Effektives_Volumen_Unified_PIAFix.exe")
        sys.exit(1)
    if needs_piafix:
        print(f"   [Manager] Pre-Run-Checks OK: Unified-EXE + PIAFix-EXE + CSV + ANSYS-EXE vorhanden.")
    else:
        print(f"   [Manager] Pre-Run-Checks OK: Unified-EXE, CSV, ANSYS-EXE vorhanden.")

    # Pre-Run-Kollisions-Detektor — verhindert dass zwei Cases mit identischem
    # full_name (z.B. selbe STRESS/INT/REF/CRIT/LOC, aber unterschiedlichem PIA_FIX) im
    # gleichen Batch laufen und sich gegenseitig die Output-Files ueberschreiben.
    # Hintergrund: full_name = {CASE}-{GEOM}-{STRESS}-{INT}-{REF}-{CRIT}-{LOC}-{MESH}.
    # PIA_FIX ist bewusst NICHT im full_name (Numerik-Variante, kein Methoden-Token).
    # Workaround fuer dieses Setup: PIA_FIX=0 und PIA_FIX=1 in zwei separate Batches
    # splitten (Option B aus dem v19.0-Plan).
    seen_full_names = {}  # full_name -> (csv-Zeile, pia_fix-Flag)
    try:
        with open(CSV_FILE, 'r', encoding='utf-8') as _f:
            _reader = csv.DictReader(_f, delimiter=';')
            for _line_idx, _row in enumerate(_reader, start=2):  # Header = Zeile 1
                # Nachbau der full_name-Logik aus der Case-Iteration (kompakt, nur Tokens)
                _case_id  = _row.get('CASE', '').strip()
                _case_x   = _row.get('CASE_X', '').strip()
                _case_y   = _row.get('CASE_Y', '').strip()
                _case_z   = _row.get('CASE_Z', '').strip()
                _stress   = _row.get('STRESS', '').strip()
                _int_type = _row.get('INT', '').strip()
                _ref      = _row.get('REF', '').strip()
                _crit     = _row.get('CRIT', '').strip()
                _loc      = _row.get('LOC', '').strip()
                _mesh_x   = _row.get('MESH_X', '').strip()
                _mesh_y   = _row.get('MESH_Y', '').strip()
                _mesh_z   = _row.get('MESH_Z', '1').strip()
                _pia_fix  = _row.get('PIA_FIX', '0').strip() or '0'
                if not _case_id:
                    continue  # leere Zeile, ueberspringen
                _mesh_str = f"{_mesh_x}x{_mesh_y}"
                try:
                    if int(_mesh_z) > 1:
                        _mesh_str += f"x{_mesh_z}"
                except ValueError:
                    pass
                _geom_str = f"{_case_x}x{_case_y}x{_case_z}"
                _fname = mn.build_case_id(_case_id, _geom_str, _stress, _int_type, _ref, _crit, _loc, _mesh_str)
                if _fname in seen_full_names:
                    _prev_line, _prev_pia = seen_full_names[_fname]
                    print(f"!!! [Manager] FATAL: Doppelter full_name '{_fname}' in {CSV_FILE}.")
                    print(f"!!!         Cases in CSV-Zeile {_prev_line} (PIA_FIX={_prev_pia}) und")
                    print(f"!!!                       Zeile {_line_idx} (PIA_FIX={_pia_fix}) kollidieren.")
                    print(f"!!!         Beide wuerden in den gleichen case_dir schreiben")
                    print(f"!!!         und sich die Output-Files ueberschreiben.")
                    print(f"!!!         Aktion: PIA_FIX=0 und PIA_FIX=1 in zwei separate")
                    print(f"!!!                 Batches splitten")
                    sys.exit(1)
                seen_full_names[_fname] = (_line_idx, _pia_fix)
        print(f"   [Manager] Kollisions-Check OK: {len(seen_full_names)} unique full_names in CSV.")
    except FileNotFoundError:
        # CSV-Existenz wurde schon im Pre-Run-Check verifiziert; ignoriere hier
        pass

    # 1. Zeitstempel Ordner erstellen
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
    batch_dir = os.path.join(RUNS_DIR, timestamp)

    if not os.path.exists(batch_dir):
        os.makedirs(batch_dir)

    print(f"\n=== [Manager] Neuer Batch: {batch_dir} ===")

    # 2. Skripte kopieren (Flattening)
    print("\n   [Manager] Sammle und kopiere Skripte...")

    for rel_path in SOURCE_FILES:
        src = os.path.join(WORKING_DIR, rel_path)
        filename = os.path.basename(rel_path)
        dst = os.path.join(batch_dir, filename)

        if os.path.exists(src):
            shutil.copy(src, dst)
        else:
            print(f"!!! [Manager] WARN: Quelldatei nicht gefunden: {rel_path}")

    # 2b. Name-Mapping-Datei (apdl_name <-> full_name)
    mapping_path = os.path.join(batch_dir, "NAME_MAPPING.csv")
    mapping_file = open(mapping_path, "w")
    mapping_file.write("apdl_name;full_name\n")

    # 3. CSV Lesen und Runs starten
    with open(CSV_FILE, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter=';')

        for row in reader:
            # Token Parsing (Methoden-Baukasten Stress|Int|Ref)
            case_id   = row['CASE']
            case_x    = row['CASE_X']
            case_y    = row['CASE_Y']
            case_z    = row['CASE_Z']
            stress    = row['STRESS'].strip().upper()
            int_type  = row['INT'].strip().upper()
            ref       = row['REF'].strip().upper()
            crit      = row['CRIT'].strip().upper()
            loc       = row['LOC'].strip().upper()
            mesh_x    = row['MESH_X']
            mesh_y    = row['MESH_Y']
            mesh_z    = row['MESH_Z']
            load_n    = row['LOAD_N']
            sig_0_abs = row['SIG_0_ABS']
            sig_0     = row['SIG_0']
            vtk_flag  = row.get('VTK', '0').strip()
            err_ext_str = row.get('ERR_EXT', '').strip()  # per-Case Extended Error Decomp
            pia_fix_str = row.get('PIA_FIX', '0').strip()  # per-Case PIA-at-Gauss-Point Korrektur

            # Zentrale Validierung der Methodenkombination (R1, R2, LOC=A-Regeln)
            ok, reason = mn.validate_method_combination(stress, int_type, ref, loc, crit)
            if not ok:
                print(f"!!! SKIP {case_id}: {reason}")
                continue

            # APDL-Bruecke: STRESS/REF -> my_avg/my_norm fuer current_params.inp
            avg_mode = mn.legacy_avg_from_stress(stress)   # 'RAW' oder 'NOD'
            norm     = mn.legacy_norm_from_ref(ref)        # 'EM', 'GP' oder 'NOD'

            # Mesh String
            mesh_str = f"{mesh_x}x{mesh_y}"
            if int(mesh_z) > 1:
                mesh_str += f"x{mesh_z}"

            # Geometry String (vollstaendig fuer full_name, verkuerzt fuer apdl_name)
            geom_str = f"{case_x}x{case_y}x{case_z}"
            apdl_geom = _shorten_geom_for_apdl(case_id, geom_str)

            # Case-ID = {CASE}-{GEOM}-{STRESS}-{INT}-{REF}-{CRIT}-{LOC}-{MESH}
            full_name = mn.build_case_id(case_id, geom_str, stress, int_type, ref, crit, loc, mesh_str)

            # APDL-Name: Ohne CRIT/LOC + ggf. verkuerzter Geom-String (ANSYS Limit 32 Zeichen)
            apdl_name = f"{case_id}-{apdl_geom}-{stress}-{int_type}-{ref}-{mesh_str}"
            if len(apdl_name) > 32:
                print(f"!!! WARNUNG: APDL-Name '{apdl_name}' hat {len(apdl_name)} Zeichen (Limit: 32)!")

            mapping_file.write(f"{apdl_name};{full_name}\n")

            print(f"\n   [Manager] --> Starte Case: {full_name}")

            # 4. Unterordner erstellen
            case_dir = os.path.join(batch_dir, full_name)
            os.makedirs(case_dir, exist_ok=True)

            os.makedirs(os.path.join(case_dir, "plots"), exist_ok=True)
            os.makedirs(os.path.join(case_dir, "tables"), exist_ok=True)

            # 5. Params schreiben
            is_gauss = 1 if int_type.startswith("G") else 0
            is_em    = 1 if int_type == "EM" else 0
            gauss_ord = int(int_type[1:]) if is_gauss else 0

            param_file = os.path.join(case_dir, "current_params.inp")
            with open(param_file, "w") as pf:
                pf.write("! Auto-generated\n")
                pf.write(f"my_case   = '{case_id}'\n")
                pf.write(f"my_int    = '{int_type}'\n")
                pf.write(f"my_norm   = '{norm}'\n")
                pf.write(f"my_avg    = '{avg_mode}'\n")
                pf.write(f"my_crit   = '{crit}'\n")
                pf.write(f"my_loc    = '{loc}'\n")
                pf.write(f"nelx      = {mesh_x}\n")
                pf.write(f"nely      = {mesh_y}\n")
                pf.write(f"nelz      = {mesh_z}\n")
                pf.write(f"load_F    = {load_n}\n")
                pf.write(f"full_name = '{apdl_name}'\n")
                pf.write(f"is_gauss  = {is_gauss}\n")
                pf.write(f"is_em     = {is_em}\n")
                pf.write(f"gauss_ord = {gauss_ord}\n")
                pf.write(f"case_x    = {case_x}\n")
                pf.write(f"case_y    = {case_y}\n")
                pf.write(f"case_z    = {case_z}\n")
                pf.write(f"my_vtk    = {vtk_flag}\n")
                # PIA_FIX-Schalter fuer APDL *IF-Branch in post_GAUSS_unified.inp.
                # 0 = Legacy-EXE (Default), 1 = PIAFix-EXE.
                # Wirkt nur bei INT=G* (silent ignore bei EM, der EM-Pfad braucht keinen Gauss-Fix).
                pf.write(f"my_pia_fix = {pia_fix_str if pia_fix_str else '0'}\n")

            # 6. ANSYS STARTEN (apdl_name fuer Jobname, <= 32 Zeichen)
            out_file_rel = f"{apdl_name}.out"
            input_script = "..\\02-main.inp"

            cmd = [
                ANSYS_EXE_PATH,
                "-dir", case_dir,
                "-j", apdl_name,
                "-i", input_script,
                "-o", out_file_rel,
                "-np", str(NUM_CORES),
                "-b", "-g", "off"
            ]

            try:
                subprocess.run(
                    cmd,
                    shell=True,
                    cwd=case_dir,
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE
                )
            except subprocess.CalledProcessError as e:
                print(f"!!! FEHLER bei ANSYS Run {full_name} !!!")
                if e.stderr:
                    print(e.stderr.decode('cp1252', errors='ignore'))

            # 7. POST-PROCESSING (Der Analyst) — v14: STRESS/REF statt AVG/NORM
            if RUN_POST_PROCESSING and post_module:
                try:
                    post_module.run_analysis_for_case(
                        run_dir=case_dir,
                        full_name=full_name,
                        apdl_name=apdl_name,
                        int_type=int_type,
                        mesh_y=mesh_y,
                        stress=stress,
                        ref=ref,
                        sig_0_abs=(sig_0_abs.strip().upper() == 'TRUE'),
                        sig_0=float(sig_0),
                        case_x=float(case_x),
                        case_y=float(case_y),
                        case_z=float(case_z),
                        load_n=float(load_n),
                        ansys_exe_path=ANSYS_EXE_PATH,
                        num_cores=NUM_CORES,
                        mesh_x=int(mesh_x),
                        mesh_z=int(mesh_z),
                        crit=crit,
                        loc=loc,
                        export_vtk=(vtk_flag == '1'),
                        err_ext=(True if err_ext_str == '1' else (False if err_ext_str == '0' else None)),
                        pia_fix=(pia_fix_str == '1'),  # PIA-at-Gauss-Point Korrektur
                    )
                except Exception as e:
                    print(f"!!! FEHLER im Post-Analyst: {e}")

            # 8. PER-CASE CLEANUP: ANSYS-Dateien in Unterordner verschieben
            try:
                _cleanup_case_dir(case_dir)
            except Exception as e:
                print(f"   [Cleanup] WARNUNG: Case-Cleanup fehlgeschlagen: {e}")

    mapping_file.close()
    print(f"   [Manager] Name-Mapping gespeichert: {mapping_path}")

    # 9. BATCH-LEVEL CLEANUP: Berechnungscodes in Unterordner verschieben
    try:
        _cleanup_batch_dir(batch_dir)
    except Exception as e:
        print(f"   [Cleanup] WARNUNG: Batch-Cleanup fehlgeschlagen: {e}")

    print("\n=== Fertig. Alle Runs abgeschlossen. ===")

if __name__ == "__main__":
    run_study()