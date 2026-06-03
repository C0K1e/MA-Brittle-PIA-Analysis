# config_local.example.py — Vorlage fuer lokale Maschinen-Konfiguration
# ============================================================================
# SO BENUTZEN:
#   1. Diese Datei kopieren zu  ->  config_local.py
#   2. Die ANSYS-EXE-Pfade und Kernzahlen fuer DEINEN Rechner eintragen.
#
# config_local.py ist in .gitignore und wird nie committet/veroeffentlicht.
# Diese .example-Datei IST im Repo und dient nur als Vorlage/Dokumentation,
# welche Variablen 01-run_manager.py erwartet.
#
# Profil-Auswahl ueber KONFIGURATION_SANSYS2 in 01-run_manager.py:
#   True  -> Server-Profil (die *_SERVER-Werte unten)
#   False -> Laptop-Profil (die *_LAPTOP-Werte unten)
# ============================================================================

# --- Server-Profil (KONFIGURATION_SANSYS2 = True) ---
ANSYS_EXE_PATH_SERVER = r"C:\Pfad\zu\ANSYS Inc\v251\ansys\bin\winx64\ANSYS251.exe"
NUM_CORES_SERVER = 16

# --- Laptop-Profil (KONFIGURATION_SANSYS2 = False) ---
ANSYS_EXE_PATH_LAPTOP = r"C:\Pfad\zu\ANSYS Inc\ANSYS Student\v252\ansys\bin\winx64\ANSYS252.exe"
NUM_CORES_LAPTOP = 4
