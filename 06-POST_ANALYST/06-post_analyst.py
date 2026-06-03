import os
import sys
import math
import csv
import re

# sys.path Fix fuer dynamisches Laden via importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Sub-Folder-Struktur — zentraler sys.path-Setup, damit alle bestehenden
# Sibling-Imports (from analytical_helper import ..., from report_generator import ...,
# from parse_presol_output import ..., import method_naming as mn) unveraendert
# weiterarbeiten, obwohl die Module jetzt in Sub-Foldern liegen.
_HERE = os.path.dirname(os.path.abspath(__file__))
for _sub in ("analytical", "error_decomposition", "output", "stress_parsing", "helpers"):
    sys.path.insert(0, os.path.join(_HERE, _sub))

from analytical_helper import (get_analytical_veff, get_analytical_aeff,
                               get_analytical_smax_ref, LOADCASE_REGISTRY)
from report_generator import save_json_metadata, generate_markdown_report
from parse_presol_output import (parse_presol_to_csv, read_stress_csv, compute_element_means,
                                 read_element_file_with_ids, write_fortran_stress_file,
                                 max_s1_from_clean_csv)
import method_naming as mn  # zentrales Naming-Modul (STRESS/INT/REF + APDL-Bruecke)
from pf_helper import pf_from_hazard  # numerisch robuste Pf-Berechnung
import subprocess
import shutil

# Fehlerzerlegung
INTEGRATE_ERROR_DIFFERENTIATION = True
REFERENCE_CSV = os.path.join("06-POST_ANALYST", "REFERENCES", "Reference_Gauss26.csv")

"""
Post-Processing Analyst fuer EM/Gauss-Integrationen.

v14: Methoden-Baukasten Stress|Int|Ref (siehe method_naming.py).

STRESS-Parameter (Spannungsdarstellung, externe v14-API):
    - "RAW" : Ungemittelte Element-Spannungen (ESEL-Trick / PRESOL)
    - "AVG" : Nodal-gemittelte Spannungen (Standard ANSYS)

INT-Parameter (Integrationsverfahren):
    - "EM"            : Elemental Mean (direkt aus ANSYS-Schleife)
    - "G1"..."G9", "G26" : Gauss-Integration via Fortran

REF-Parameter (Referenzspannung, externe v14-API):
    - "EMX" : Element-Mean-Maximum   (kompatibel mit INT=EM)
    - "GPX" : Gauss-Punkt-Maximum    (kompatibel mit INT=G*)
    - "NDX" : Knoten-Maximum         (kompatibel mit allen INT)

Interne Tech-State (APDL-Bruecke, NICHT Teil der externen API):
    self.avg  ∈ {RAW, NOD}            <-> my_avg in current_params.inp
    self.norm ∈ {EM, GP, NOD}         <-> my_norm in current_params.inp
"""

class ResultAnalyst:
    def __init__(self, run_dir, case_name, mesh_y,
                 case_x=0, case_y=0, case_z=0,
                 stress="AVG", ref="EMX", sig_0_abs=True, sig_0=100.0, load_n=1.0,
                 ansys_exe_path="", num_cores=0, apdl_name=None,
                 int_type=None, mesh_x=0, mesh_z=1, crit="PIA", loc="V",
                 pia_fix=False):
        self.run_dir = run_dir
        self.case_name = case_name
        # APDL-Name (<=32 Zeichen) fuer ANSYS-Datei-Lookups
        self.apdl_name = apdl_name or case_name
        self.tables_dir = os.path.join(run_dir, "tables")
        self.int_type = (int_type or "UNKNOWN").upper()
        self.mesh_x = int(mesh_x) if mesh_x else 0
        self.mesh_y = int(mesh_y)
        self.mesh_z = int(mesh_z) if mesh_z else 1

        # Geometrie-Parameter aus CSV
        self.case_x = float(case_x) if case_x else 0.0
        self.case_y = float(case_y) if case_y else 0.0
        self.case_z = float(case_z) if case_z else 0.0
        self.load_n = float(load_n) if load_n else 1.0

        # Externe API — STRESS/REF (Methoden-Baukasten)
        self.stress = stress.strip().upper()
        self.ref    = ref.strip().upper()
        # Interne Tech-State (APDL-Bruecke) — wird ueberall im Code unten verwendet:
        #   self.avg  ∈ {RAW, NOD}    fuer ANSYS-Mittelungs-Modus
        #   self.norm ∈ {EM, GP, NOD} fuer smax_norm-Quelle (siehe _read_gauss_norm_smax etc.)
        self.avg  = mn.legacy_avg_from_stress(self.stress)
        self.norm = mn.legacy_norm_from_ref(self.ref)

        # Sigma0-Steuerung per Case (aus CSV)
        self.sig_0_abs = sig_0_abs
        self.sig_0_value = float(sig_0) if sig_0 else 100.0
        self.crit = crit.upper() if crit else "PIA"
        self.loc = (loc or "V").upper()  # Auswertungsdomaene V (Volumen) oder A (Flaeche)
        # PIA_FIX-Flag — bei INT=G* steuert die EXE-Wahl (Legacy vs. PIAFix).
        # Bei INT=EM ist der Flag wirkungslos (silent ignore — EM hat keine Gauss-Punkte).
        self.pia_fix = bool(pia_fix)

        self.m_start = 1
        self.m_max = 50

        # --- Spannungsreferenzen ---
        self.smax_nodal = 0.0   # Maximale nodale Spannung (unabhaengig von NORM)
        self.smax_norm = 0.0    # Normierungswert (abhaengig von NORM: NOD/GP/EM)
        self.sigma0 = 0.0       # Sigma0 — ein Wert (absolut aus CSV oder relativ zu smax_nodal)

        self.vtot_num = 0.0   # Numerisches Volumen
        self.atot_num = 0.0   # Numerische Gesamtflaeche (LOC=A)

        # Exakte Analytik
        self.vtot_exact = 0.0
        self.atot_exact = None  # Analytische Gesamtflaeche (Wenn None: case-spezifisch nicht verfuegbar)
        self.case_type = "UNKNOWN"

        # Symmetrie-Faktor initialisieren (Standard = 1.0 fuer Vollmodelle)
        self.symmetry_factor = 1.0

        self.elements = []
        self.element_ids = []  # Parallele Liste der ANSYS-Element-IDs (fuer VTK-Export)

        # Metadaten aus Run-Manager
        self.ansys_exe_path = ansys_exe_path
        self.num_cores = num_cores
        self.run_metadata = {}

        # Initialisierung
        self._setup_geometry_and_physics()
        self.results = []

    def _compute_sigma0(self, max_stress):
        """Berechnet Sigma0 basierend auf per-Case Einstellungen aus CSV."""
        if self.sig_0_abs:
            return self.sig_0_value  # Absoluter Wert aus CSV
        else:
            return self.sig_0_value * max_stress  # Relativ zu max_stress

    def _setup_geometry_and_physics(self):
        """Bestimmt analytisches Volumen basierend auf Case-Namen und Geometrie-Parametern.

        Case-Type wird via split('-')[0] ermittelt (exaktes Matching, keine
        Substring-Reihenfolge-Abhaengigkeit wie PWHR/PWH).
        """
        # Case-Type = erstes Token im Namen (z.B. "PWH-60x0x1-..." -> "PWH")
        self.case_type = self.case_name.split("-")[0]

        entry = LOADCASE_REGISTRY.get(self.case_type)
        if entry is None:
            print(f"   [Analyst] WARNUNG: Unbekannter Case '{self.case_name}'. "
                  f"Analytik nutzt Mesh-Daten.")
            self.case_type = "UNKNOWN"
            self.vtot_exact = None
            self.atot_exact = None
            return

        self.symmetry_factor = entry["symmetry_factor"]
        self.vtot_exact = entry["vtot_exact"](self.case_x, self.case_y, self.case_z)

        atot_func = entry.get("atot_exact")
        self.atot_exact = atot_func(self.case_x, self.case_y, self.case_z) if atot_func else None

        # LOC=A nur fuer Lastfaelle mit analytischer Aeff-Formel zulassen
        if self.loc == "A":
            aeff_func = entry.get("aeff")
            if aeff_func is None:
                raise NotImplementedError(
                    f"Aeff-Analytik fuer Lastfall '{self.case_type}' nicht verfuegbar. "
                    f"In v12.0 unterstuetzt: BEAM1, BEAM2, POR. "
                    f"Lastfall '{self.case_type}' hat aeff=None in LOADCASE_REGISTRY (TODO)."
                )

    def _get_analytic_veff(self, m):
        """Delegiert an analytical_helper Modul."""
        vol_base = self.vtot_exact if self.vtot_exact else self.vtot_num
        return get_analytical_veff(
            self.case_type, m, vol_base,
            dim_x=self.case_x, dim_y=self.case_y, dim_z=self.case_z,
            crit=self.crit
        )

    def _get_analytic_aeff(self, m):
        """v12.0: Analytische effektive Flaeche fuer LOC=A.

        Delegiert an analytical_helper.get_analytical_aeff(). Wirft NotImplementedError
        wenn der Lastfall keine Aeff-Formel hat (PWH, PWHR, 3PB).
        """
        area_base = self.atot_exact if self.atot_exact else self.atot_num
        val = get_analytical_aeff(
            self.case_type, m, area_base,
            dim_x=self.case_x, dim_y=self.case_y, dim_z=self.case_z,
            crit=self.crit
        )
        if val is None:
            raise NotImplementedError(
                f"Aeff-Analytik fuer '{self.case_type}' nicht verfuegbar (m={m})."
            )
        return val

    def _override_smax_nodal_from_raw_csv(self):
        """v14.2: Bei STRESS=RAW smax_nodal aus -S-clean.csv max(s1) ueberschreiben.

        Motivation: Der APDL-Header schreibt nodal-gemittelte S1-Werte
        (`*GET,smax,SORT,0,MAX` nach `NSORT,S,1`). Bei STRESS=RAW soll auch
        die Normierungs-Quelle RAW sein, damit alle Pipeline-Stufen domain-
        konsistent arbeiten — sonst entsteht ein Inkonsistenz-Artefakt zwischen
        nodal-gemitteltem smax_nodal und ungemittelten Element-Knoten-Spannungen
        in der Veff/Aeff-Integration.

        Greift auf:
          - LOC=V: {apdl_name}-S-clean.csv
          - LOC=A: {apdl_name}-A-S-clean.csv
        Aktiv nur bei self.avg == "RAW" (= STRESS=RAW). Bei AVG-Cases keine Aktion.

        Bei numerischer Identitaet (|raw_max - header_max| < 1e-6) erfolgt
        keine Korrektur und kein Print — das ist typisch fuer glatte
        Spannungsfelder am Maximum (z.B. POR equi-biaxial bei r=0).
        """
        if self.avg != "RAW":
            return  # nur bei STRESS=RAW relevant

        if self.loc == "A":
            clean_csv = os.path.join(self.tables_dir, f"{self.apdl_name}-A-S-clean.csv")
        else:
            clean_csv = os.path.join(self.tables_dir, f"{self.apdl_name}-S-clean.csv")

        if not os.path.exists(clean_csv):
            return

        raw_max = max_s1_from_clean_csv(clean_csv)
        if raw_max is None or raw_max <= 0:
            return

        # Numerische Identitaet -> keine Korrektur, kein Log-Spam
        if self.smax_nodal > 0 and abs(raw_max - self.smax_nodal) < 1e-6:
            return

        prev = self.smax_nodal
        self.smax_nodal = raw_max
        # smax_norm bei NDX-Pfad nachfuehren (smax_norm = smax_nodal bei NORM=NOD)
        norm_upper = self.norm.upper() if self.norm else "NOD"
        if norm_upper == "NOD":
            self.smax_norm = raw_max
        # sigma0 ggf. neu berechnen (relativer Modus)
        if not self.sig_0_abs:
            self.sigma0 = self._compute_sigma0(raw_max)
        print(f"   [Analyst] STRESS=RAW Konsistenz: smax_nodal Header={prev:.4e} "
              f"-> RAW-Knoten-Max={raw_max:.4e} (aus {os.path.basename(clean_csv)})")

    def _parse_header_smax(self, filepath):
        """Liest SMAX aus dem Header der Datei."""
        smax = None
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    m = re.search(r"(?:smax|SMAX)\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)", line, re.IGNORECASE)
                    if m:
                        smax = float(m.group(1))
                        break
        except FileNotFoundError:
            return None
        return smax

    def _read_gauss_norm_smax(self):
        """Liest SMAX_GAUSS_VOL/SURF aus GauszInfo-File (nur bei NORM=GP).

        v13.1: Domain-spezifisch ohne Legacy-Fallback.
            LOC=V → SMAX_GAUSS_VOL
            LOC=A → SMAX_GAUSS_SURF
        Datei liegt in tables/ (v5.0+) oder run_dir (legacy).
        """
        filepath = os.path.join(self.tables_dir, f"VEFF_{self.apdl_name}_GauszInfo.out")
        if not os.path.exists(filepath):
            filepath = os.path.join(self.run_dir, "effVol_GauszInfo.out")
        if not os.path.exists(filepath):
            return None

        loc_upper = (getattr(self, 'loc', 'V') or 'V').upper()
        preferred_field = 'SMAX_GAUSS_SURF' if loc_upper == 'A' else 'SMAX_GAUSS_VOL'

        fields = {}
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                m = re.match(
                    r"\s*(SMAX_GAUSS_(?:VOL|SURF|LINE))\s*=\s*"
                    r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)",
                    line)
                if m:
                    fields[m.group(1)] = float(m.group(2))

        if preferred_field in fields and fields[preferred_field] > 0:
            print(f"   [Analyst] GaussNorm: {preferred_field} = {fields[preferred_field]:.4e} (LOC={loc_upper})")
            return fields[preferred_field]
        return None

    def _load_mesh_properties(self):
        """Liest Element-Daten (EM-Pfad).

        Zwei-Datei-Schema (v8.3, EM-RAW):  {apdl_name}-V.out  +  {apdl_name}-S.out
        Fallback Single-File (EM-NOD / alt): {apdl_name}.out   (6 Spalten: Elem,Vol,S1,S2,S3,dS1)
        """
        vol_path    = os.path.join(self.tables_dir, f"{self.apdl_name}-V.out")
        stress_path = os.path.join(self.tables_dir, f"{self.apdl_name}-S.out")
        legacy_path = os.path.join(self.tables_dir, f"{self.apdl_name}.out")

        use_two_file = os.path.exists(vol_path)

        if use_two_file:
            # ===== Neues Zwei-Datei-Schema (EM-RAW ab v8.3) =====
            if not os.path.exists(stress_path):
                print(f"   [Analyst] CRITICAL: -V.out gefunden, aber -S.out fehlt! {stress_path}")
                return False

            # 1. SMAX aus Header der -V-Datei lesen
            if self.smax_nodal == 0.0:
                val = self._parse_header_smax(vol_path)
                if val:
                    self.smax_nodal = val
                    self.sigma0 = self._compute_sigma0(self.smax_nodal)

            # 2. Volumen-Dictionary aus -V.out
            vol_dict = {}
            n_vol_rows = 0
            sum_raw_vol = 0.0
            n_collisions = 0
            try:
                with open(vol_path, 'r', encoding='utf-8', errors='ignore') as f:
                    start_reading = False
                    for line in f:
                        if "Elem" in line and "Vol" in line:
                            start_reading = True
                            continue
                        if not start_reading:
                            continue
                        parts = line.split()
                        if len(parts) < 2:
                            continue
                        try:
                            elem_id = int(round(float(parts[0])))
                            vol     = float(parts[1])
                            n_vol_rows += 1
                            sum_raw_vol += vol
                            if elem_id in vol_dict:
                                n_collisions += 1
                            vol_dict[elem_id] = vol
                        except ValueError:
                            continue
            except Exception as e:
                print(f"   [Analyst] Fehler beim Lesen der Volumen-Datei: {e}")
                return False

            if not vol_dict:
                print(f"   [Analyst] CRITICAL: Keine Elementdaten in {vol_path}")
                return False

            if n_collisions > 0:
                raise ValueError(
                    f"Element-ID-Kollisionen in -V.out: {n_collisions} Duplikate bei "
                    f"{n_vol_rows} Zeilen — APDL schreibt Elementnummern im E-Format "
                    f"statt I-Format (post_EM_RAW.inp: I14 statt E12.5 verwenden)."
                )

            print(f"   [Analyst] V-out: {n_vol_rows} Zeilen, {len(vol_dict)} unique IDs, "
                  f"Summe={sum(vol_dict.values()):.6e}")

            # 3. Elementmittelwerte der Hauptspannungen aus -S.out (PRESOL)
            #    Bereinigte CSV on-demand erzeugen (nur beim ersten Lauf noetig)
            clean_csv_path = stress_path.replace('-S.out', '-S-clean.csv')
            if not os.path.exists(clean_csv_path):
                print(f"   [Analyst] Erzeuge bereinigtes Stress-CSV: {os.path.basename(clean_csv_path)}")
                parse_presol_to_csv(stress_path, clean_csv_path)

            raw_data    = read_stress_csv(clean_csv_path)
            stress_dict = compute_element_means(raw_data)

            if not stress_dict:
                print(f"   [Analyst] CRITICAL: Keine Spannungsdaten aus {clean_csv_path} gelesen")
                return False

            # STRESS=RAW Konsistenz — smax_nodal aus PRESOL max(s1) lesen
            self._override_smax_nodal_from_raw_csv()

            # 4. Zusammenfuehren: (vol, s1, s2, s3) pro Element
            self.elements = []
            self.element_ids = []
            n_missing = 0
            for elem_id, vol in vol_dict.items():
                if elem_id in stress_dict:
                    s1, s2, s3 = stress_dict[elem_id]
                    self.elements.append((vol, s1, s2, s3))
                    self.element_ids.append(elem_id)
                else:
                    n_missing += 1
            if n_missing > 0:
                print(f"   [Analyst] Info: {n_missing} Elemente ohne Spannungsdaten (ignoriert).")

        else:
            # ===== Legacy Single-File-Schema (EM-NOD oder alte EM-RAW-Runs) =====
            if not os.path.exists(legacy_path):
                print(f"   [Analyst] CRITICAL: Keine Element-Datei gefunden!")
                print(f"               Gesucht: {vol_path}")
                print(f"               Gesucht: {legacy_path}")
                return False

            print(f"   [Analyst] Lese Legacy-Format (Single-File): {os.path.basename(legacy_path)}")

            # SMAX aus Header
            if self.smax_nodal == 0.0:
                val = self._parse_header_smax(legacy_path)
                if val:
                    self.smax_nodal = val
                    self.sigma0 = self._compute_sigma0(self.smax_nodal)

            # Elemente lesen (Spalten: Elem, Vol, S1, S2, S3, ...)
            self.elements = []
            self.element_ids = []
            n_legacy_rows = 0
            n_legacy_collisions = 0
            seen_legacy_ids = set()
            try:
                with open(legacy_path, 'r', encoding='utf-8', errors='ignore') as f:
                    start_reading = False
                    for line in f:
                        if "Elem" in line and "Vol" in line:
                            start_reading = True
                            continue
                        if not start_reading:
                            continue
                        parts = line.split()
                        if len(parts) < 5:
                            continue
                        try:
                            eid = int(round(float(parts[0])))
                            vol = float(parts[1])
                            s1  = float(parts[2])
                            s2  = float(parts[3])
                            s3  = float(parts[4])
                            n_legacy_rows += 1
                            if eid in seen_legacy_ids:
                                n_legacy_collisions += 1
                            seen_legacy_ids.add(eid)
                            self.elements.append((vol, s1, s2, s3))
                            self.element_ids.append(eid)
                        except ValueError:
                            continue
            except Exception as e:
                print(f"   [Analyst] Fehler beim Lesen der Mesh-Daten: {e}")
                return False

            if n_legacy_collisions > 0:
                raise ValueError(
                    f"Element-ID-Kollisionen in {os.path.basename(legacy_path)}: "
                    f"{n_legacy_collisions} Duplikate bei {n_legacy_rows} Zeilen — "
                    f"APDL schreibt Elementnummern im E-Format statt I-Format."
                )

        if not self.elements:
            return False

        if self.vtot_num == 0.0:
            self.vtot_num = sum(e[0] for e in self.elements) * self.symmetry_factor

        # --- Bestimmung von smax_norm (NORM-Dispatch, unveraendert) ---
        max_s1_found = 0.0
        for _, s1, _, _ in self.elements:
            if s1 > max_s1_found:
                max_s1_found = s1

        norm_upper = self.norm.upper() if self.norm else "EM"

        if norm_upper == "NOD":
            if self.smax_nodal > 0:
                self.smax_norm = self.smax_nodal
            else:
                self.smax_norm = max_s1_found
                print(f"   [Analyst] WARNUNG: NORM=NOD gesetzt, aber smax_nodal fehlt im Header.")
                print(f"   [Analyst] Fallback auf Element-Maximum: {self.smax_norm:.4e}")
        elif norm_upper == "EM":
            self.smax_norm = max_s1_found
        elif norm_upper == "GP":
            self.smax_norm = max_s1_found
        else:
            self.smax_norm = max_s1_found
            print(f"   [Analyst] WARNUNG: Unbekannter NORM-Wert '{self.norm}'. Zulaessig: EM, NOD, GP")
            print(f"   [Analyst] Fallback auf NORM=EM (Element-Maximum)")

        if self.smax_nodal == 0.0:
            print("   [Analyst] Info: Smax im Header nicht gefunden, nutze Mesh-Max als Referenz.")
            self.smax_nodal = self.smax_norm
            self.sigma0 = self._compute_sigma0(self.smax_nodal)

        return True

    def _load_gauss_properties(self):
        """Liest Mesh-Infos aus Gauss-spezifischen Dateien (kein EM-File noetig).
        - Element-Anzahl: aus VEFF_*_Elemente.out (Zeilenanzahl)
        - smax_nodal: aus VEFF Header
        - vtot_num: aus VEFF m=0 Zeile (PIA-V_eff) * symmetry_factor
        """
        # 1. Element-Anzahl zaehlen
        elem_file = os.path.join(self.tables_dir, f"VEFF_{self.apdl_name}_Elemente.out")
        self.num_elements_gauss = 0
        if os.path.exists(elem_file):
            try:
                with open(elem_file, 'r', encoding='utf-8', errors='ignore') as f:
                    self.num_elements_gauss = sum(1 for line in f if line.strip())
            except Exception as e:
                print(f"   [Analyst] Fehler beim Zaehlen der Elemente: {e}")

        # 2. smax_nodal aus VEFF Header
        veff_file = os.path.join(self.tables_dir, f"VEFF_{self.apdl_name}.out")
        if self.smax_nodal == 0.0:
            val = self._parse_header_smax(veff_file)
            if val:
                self.smax_nodal = val
                self.sigma0 = self._compute_sigma0(self.smax_nodal)

        # 3. vtot_num aus m=0 Zeile (Viertelmodell) * symmetry_factor
        if os.path.exists(veff_file):
            try:
                with open(veff_file, 'r') as f:
                    for line in f:
                        parts = line.split()
                        if len(parts) > 2 and parts[0] == '0':
                            self.vtot_num = float(parts[2]) * self.symmetry_factor
                            break
            except Exception as e:
                print(f"   [Analyst] Fehler beim Lesen von vtot aus VEFF: {e}")

        return self.vtot_num > 0

    def _run_gauss_raw_fortran(self):
        """Orchestriert Gauss-RAW: PRESOL parsen -> Stress-Datei -> Fortran RAWNodes.

        Schritt-fuer-Schritt:
          1. PRESOL parsen -> -S-clean.csv
          2. Element-Mapping lesen (effVol_Elemente.out mit ElemIDs)
          3. Fortran-Stress-Datei schreiben (effVol_RawStress.out)
          4. Dateien in run_dir kopieren (Fortran erwartet sie dort)
          5. GaussPNorm-Flag in Parameter-Datei setzen
          6. Fortran RAWNodes EXE ausfuehren
          7. Output umbenennen nach tables/

        Returns: True bei Erfolg, False bei Fehler.
        """
        tables = self.tables_dir

        # --- Schritt 1: PRESOL parsen ---
        stress_out = os.path.join(tables, f"{self.apdl_name}-S.out")
        clean_csv  = os.path.join(tables, f"{self.apdl_name}-S-clean.csv")

        if not os.path.exists(stress_out):
            print(f"   [Analyst] CRITICAL: PRESOL-Datei fehlt: {stress_out}")
            return False

        if not os.path.exists(clean_csv):
            print(f"   [Analyst] Erzeuge bereinigtes Stress-CSV fuer Gauss-RAW...")
            parse_presol_to_csv(stress_out, clean_csv)

        presol_data = read_stress_csv(clean_csv)
        if not presol_data:
            print(f"   [Analyst] CRITICAL: Keine PRESOL-Daten aus {clean_csv}")
            return False

        # STRESS=RAW Konsistenz — smax_nodal aus PRESOL max(s1) lesen
        self._override_smax_nodal_from_raw_csv()

        # --- Schritt 2: Element-Mapping lesen ---
        elem_file = os.path.join(tables, f"VEFF_{self.apdl_name}_Elemente.out")
        if not os.path.exists(elem_file):
            print(f"   [Analyst] CRITICAL: Elemente-Datei fehlt: {elem_file}")
            return False

        try:
            elem_order, node_order = read_element_file_with_ids(elem_file)
        except (ValueError, FileNotFoundError) as e:
            print(f"   [Analyst] CRITICAL: Element-Mapping fehlgeschlagen: {e}")
            return False

        # --- Schritt 3: Fortran-Stress-Datei schreiben ---
        raw_stress_path = os.path.join(self.run_dir, "effVol_RawStress.out")
        try:
            write_fortran_stress_file(presol_data, elem_order, node_order, raw_stress_path)
        except ValueError as e:
            print(f"   [Analyst] CRITICAL: Stress-Datei-Erzeugung fehlgeschlagen: {e}")
            return False

        # --- Schritt 4: Geometrie-Dateien in run_dir kopieren ---
        copy_map = {
            f"VEFF_{self.apdl_name}_Elemente.out":   "effVol_Elemente.out",
            f"VEFF_{self.apdl_name}_Faces.out":      "effVol_Faces.out",
            f"VEFF_{self.apdl_name}_NodeCoords.out":  "effVol_NodeCoords.out",
            f"VEFF_{self.apdl_name}_Parameter.out":  "effVol_Parameter.out",
        }
        for src_name, dst_name in copy_map.items():
            src = os.path.join(tables, src_name)
            dst = os.path.join(self.run_dir, dst_name)
            if not os.path.exists(src):
                print(f"   [Analyst] CRITICAL: Datei fehlt: {src}")
                return False
            shutil.copy2(src, dst)

        # --- Schritt 5: AVGModeRAW, GaussPNorm und smax in Parameter-Datei setzen ---
        # APDL schreibt F14.0 (mit Dezimalpunkt), Python ueberschreibt hier
        # mit sauberem Integer-Format fuer zuverlaessiges Fortran (11X,I14) Lesen.
        # smax-Zeile zusaetzlich mit RAW-Knoten-Maximum ueberschreiben, sonst
        # normalisiert Fortran die raw_stress mit dem nodal-gemittelten Header-Wert
        # ([effektivesVol_unified.f90:421-452]). Ohne diesen Patch waeren Veff/Aeff
        # bei STRESS=RAW + Gauss-Pfad mit AVG-smax normalisiert — methodisch inkonsistent.
        param_path = os.path.join(self.run_dir, "effVol_Parameter.out")
        norm_upper = self.norm.upper() if self.norm else "NOD"
        gnorm_val = 1 if norm_upper == "GP" else 0

        # RAW-Konsistenter smax aus PRESOL-CSV (gleiche Quelle wie smax_nodal-Override)
        raw_smax = max_s1_from_clean_csv(clean_csv)

        try:
            with open(param_path, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()
            for i, line in enumerate(lines):
                if line.startswith('AVGModeRAW:'):
                    # RAW-Modus immer 1 (wird nur im RAW-Pfad aufgerufen)
                    lines[i] = f"AVGModeRAW:{'1':>14s}    \n"
                elif line.startswith('GaussPNorm:'):
                    lines[i] = f"GaussPNorm:{str(gnorm_val):>14s}    \n"
                elif line.startswith('smax:') and raw_smax is not None and raw_smax > 0:
                    # smax durch RAW-Knoten-Max ersetzen.
                    # Format: (A8, ES16.8) — 8 Zeichen Label, 16 Zeichen Wert.
                    lines[i] = f"smax:   {raw_smax:16.8E}\n"
            with open(param_path, 'w', encoding='utf-8') as f:
                f.writelines(lines)
            if raw_smax is not None and raw_smax > 0:
                print(f"   [Analyst] STRESS=RAW Fortran-Konsistenz: "
                      f"effVol_Parameter.out smax ueberschrieben mit RAW-Max={raw_smax:.4e}")
        except Exception as e:
            print(f"   [Analyst] WARNUNG: Parameter-Flags konnten nicht gesetzt werden: {e}")

        # --- Schritt 6: Fortran Unified EXE ausfuehren ---
        # Legacy-Fallback auf RAWNodes.exe entfernt. Unified ist Pflicht
        # seit v13.0 — Pre-Run-Check in 01-run_manager.py verifiziert
        # die Unified-EXE schon vor dem ANSYS-Start. Falls die EXE hier trotzdem
        # fehlt, ist es ein echter Pipeline-Fehler (FATAL).
        # PIA_FIX-Flag steuert die EXE-Wahl im RAW-Pfad (analog zum APDL
        # *IF-Branch im NOD-Pfad in post_GAUSS_unified.inp).
        if self.pia_fix:
            exe_name_unified = "Effektives_Volumen_Unified_PIAFix.exe"
        else:
            exe_name_unified = "Effektives_Volumen_Unified.exe"
        exe_path = os.path.join(self.run_dir, "..", exe_name_unified)
        if not os.path.exists(exe_path):
            print(f"   [Analyst] CRITICAL: Fortran-EXE nicht gefunden:")
            print(f"               Pfad: {os.path.join(self.run_dir, '..', exe_name_unified)}")
            if self.pia_fix:
                print(f"               PIA_FIX=1 angefordert. Aktion: PIAFix-EXE compilieren")
                print(f"               ifx /O2 overhead_unified_piafix.f90 effektivesVol_unified_piafix.f90 /Fe:Effektives_Volumen_Unified_PIAFix.exe")
            else:
                print(f"               Aktion: Re-Compile noetig (ifx /O2 overhead_unified.f90 effektivesVol_unified.f90 /Fe:Effektives_Volumen_Unified.exe)")
            return False

        print(f"   [Analyst] Starte Fortran ({exe_name_unified}, RAW-Modus)...")
        try:
            result = subprocess.run(
                [exe_path],
                cwd=self.run_dir,
                capture_output=True, text=True, timeout=10000
            )
            if result.returncode != 0:
                print(f"   [Analyst] Fortran {exe_name_unified} Fehler (RC={result.returncode}):")
                if result.stderr:
                    print(f"   [Analyst] stderr: {result.stderr[:500]}")
                if result.stdout:
                    print(f"   [Analyst] stdout: {result.stdout[:500]}")
                return False
        except subprocess.TimeoutExpired:
            print(f"   [Analyst] CRITICAL: Fortran {exe_name_unified} Timeout (>10000)")
            return False
        except Exception as e:
            print(f"   [Analyst] CRITICAL: Fortran {exe_name_unified} Aufruf fehlgeschlagen: {e}")
            return False

        # --- Schritt 7: Output umbenennen nach tables/ ---
        effvol_out = os.path.join(self.run_dir, "effVol.out")
        veff_target = os.path.join(tables, f"VEFF_{self.apdl_name}.out")
        if os.path.exists(effvol_out):
            shutil.move(effvol_out, veff_target)
            print(f"   [Analyst] Fortran RAWNodes Output: {os.path.basename(veff_target)}")
        else:
            print(f"   [Analyst] CRITICAL: effVol.out nicht erzeugt!")
            return False

        # GauszInfo umbenennen (nur bei GaussPNorm)
        gauszinfo_out = os.path.join(self.run_dir, "effVol_GauszInfo.out")
        if os.path.exists(gauszinfo_out):
            gauszinfo_target = os.path.join(tables, f"VEFF_{self.apdl_name}_GauszInfo.out")
            shutil.move(gauszinfo_out, gauszinfo_target)
            print(f"   [Analyst] GauszInfo: {os.path.basename(gauszinfo_target)}")

        # --- Schritt 8: Intermediate Cleanup ---
        # effVol_RawStress.out nach tables/ verschieben (Archiv)
        raw_stress_src = os.path.join(self.run_dir, "effVol_RawStress.out")
        if os.path.exists(raw_stress_src):
            raw_stress_dst = os.path.join(tables, f"VEFF_{self.apdl_name}_RawStress.out")
            shutil.move(raw_stress_src, raw_stress_dst)
            print(f"   [Analyst] RawStress archiviert: {os.path.basename(raw_stress_dst)}")

        # Intermediate Kopien loeschen (Originale liegen bereits als VEFF_* in tables/)
        for fname in ["effVol_Elemente.out", "effVol_Faces.out",
                       "effVol_NodeCoords.out", "effVol_Parameter.out"]:
            fpath = os.path.join(self.run_dir, fname)
            if os.path.exists(fpath):
                try:
                    os.remove(fpath)
                except OSError:
                    pass

        print(f"   [Analyst] Gauss-RAW Fortran RAWNodes erfolgreich abgeschlossen.")
        return True

    def analyze_gauss(self):
        """Gauss-Analyse."""
        # RAWNodes: Gauss-RAW — Python orchestriert Fortran vor der Auswertung
        if self.avg == "RAW":
            if not self._run_gauss_raw_fortran():
                print("   [Analyst] Gauss-RAW: Fortran-Pipeline fehlgeschlagen.")
                return False

        if not self._load_gauss_properties():
            print("   [Analyst] Gauss: Konnte Gauss-Properties nicht laden.")

        filename = f"VEFF_{self.apdl_name}.out"
        filepath = os.path.join(self.tables_dir, filename)

        if not os.path.exists(filepath):
            print(f"   [Analyst] Gauss-Datei fehlt: {filepath}")
            return False

        if self.smax_nodal == 0.0:
            val = self._parse_header_smax(filepath)
            if val:
                self.smax_nodal = val

        # sigma0 unabhaengig von der smax_nodal-Quelle sicherstellen.
        # Bug-Hintergrund: Bei STRESS=RAW + LOC=V setzt _run_gauss_raw_fortran() via
        # _override_smax_nodal_from_raw_csv() bereits self.smax_nodal != 0, wodurch der
        # if-Block oben uebersprungen wurde. In Kombination mit SIG_0_ABS=TRUE blieb
        # self.sigma0 dann auf Default 0.0 -> Pf_num/Pf_ana == 0.0 in _compute_metrics.
        # RAW-A war zufaellig OK (Override fand -A-S-clean.csv vorher nicht und kehrte
        # frueh zurueck), AVG-Pfade sind generell OK (Override greift nicht, Block oben
        # setzt sigma0 mit). Selbes Pattern wie [_load_surface_properties():965-966].
        if self.sigma0 == 0.0 and self.smax_nodal > 0:
            self.sigma0 = self._compute_sigma0(self.smax_nodal)

        # STRESS=RAW Konsistenz auch im Gauss-Pfad. _run_gauss_raw_fortran
        # ruft den Override schon auf, aber dieser Block hier ueberschreibt smax_nodal
        # potentiell aus dem VEFF-Header (nodal-gemittelt). Daher nochmal nachziehen,
        # damit der RAW-Wert aus PRESOL die finale Quelle ist (greift nur bei AVG=RAW).
        self._override_smax_nodal_from_raw_csv()

        # LOC- und CRIT-basierter Spaltenindex
        # Fortran-Output VEFF_*.out: m | S1-V_eff | PIA-V_eff | S1-S_eff | PIA-S_eff
        #                                  Sp.1     Sp.2        Sp.3        Sp.4
        COL_MAP = {("V", "S1"): 1, ("V", "PIA"): 2,
                   ("A", "S1"): 3, ("A", "PIA"): 4}
        eff_col_idx = COL_MAP.get((self.loc, self.crit), 2)

        m_data = {}
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.split()
                if not parts: continue
                if re.fullmatch(r"\d+", parts[0]):
                    m = int(parts[0])
                    if len(parts) > eff_col_idx:
                        m_data[m] = float(parts[eff_col_idx]) * self.symmetry_factor

        # m=0 Fallback je nach Domain (Volumen oder Flaeche)
        if 0 not in m_data:
            if self.loc == "A":
                # atot_num kann aus Fortran-Output bei m=0 abgeleitet werden, sonst Analytik
                fallback_a = self.atot_num if self.atot_num > 0 else (
                    (self.atot_exact * self.symmetry_factor) if self.atot_exact else 0.0
                )
                if fallback_a > 0:
                    m_data[0] = fallback_a
            elif self.vtot_num > 0:
                m_data[0] = self.vtot_num  # vtot_num enthaelt bereits symmetry_factor

        # atot_num aus Gauss-Output ableiten, falls noch nicht gesetzt (fuer Metadaten)
        if self.loc == "A" and self.atot_num == 0.0 and 0 in m_data:
            self.atot_num = m_data[0]

        # NORM-spezifische smax_norm Bestimmung
        norm_upper = self.norm.upper() if self.norm else "NOD"  # Default NOD für Gauss

        if norm_upper == "GP":
            # GaussNorm-spezifischer Wert überschreibt Placeholder
            smax_gp = self._read_gauss_norm_smax()
            if smax_gp:
                self.smax_norm = smax_gp
                print(f"   [Analyst] Gauss-Norm: SMAX_GP = {smax_gp:.4e}")
            else:
                # Hart fail statt silent NDX-Fallback.
                # Methodenname sagt REF=GPX, Berechnung mit smax_nodal waere
                # methodisch divergent — Methodenname und Daten wuerden nicht mehr
                # zusammenpassen. Lieber Run abbrechen als falsche Ergebnisse erzeugen.
                raise RuntimeError(
                    f"REF=GPX (Gauss-Punkt-Maximum) angefordert, aber GauszInfo-Datei "
                    f"fehlt fuer Case '{self.apdl_name}'. Erwartet wurde "
                    f"'tables/VEFF_{self.apdl_name}_GauszInfo.out' oder "
                    f"'run_dir/effVol_GauszInfo.out'.\n"
                    f"Pruefe: (1) GaussPNorm=1 in effVol_Parameter.out, (2) Fortran-EXE "
                    f"hat ohne Fehler beendet (siehe fortran.log), (3) /RENAME-Zeile "
                    f"in post_GAUSS_unified.inp wurde ausgefuehrt.\n"
                    f"Run abgebrochen — kein silent NDX-Fallback (waere methodisch "
                    f"inkonsistent zum Methodenlabel)."
                )

        elif norm_upper == "NOD":
            # Setze smax_norm explizit (Gauss ruft _load_mesh_properties() nicht auf)
            if self.smax_nodal > 0:
                self.smax_norm = self.smax_nodal
            else:
                print(f"   [Analyst] WARNUNG: smax_nodal nicht verfügbar für NORM=NOD!")
                self.smax_norm = 1.0  # Notfall-Fallback

        elif norm_upper == "EM":
            # WARNUNG: EM-Normierung nicht für Gauss-Integration vorgesehen
            print(f"   [Analyst] WARNUNG: NORM=EM bei Gauss-Integration nicht vorgesehen!")
            print(f"   [Analyst] Fallback auf NORM=NOD (smax_nodal).")
            if self.smax_nodal > 0:
                self.smax_norm = self.smax_nodal
            # Falls smax_nodal fehlt, bleibt Placeholder aus _load_mesh_properties()

        else:
            # Unbekannter NORM-Wert
            print(f"   [Analyst] WARNUNG: Unbekannter NORM '{self.norm}'. Fallback auf NORM=NOD.")
            if self.smax_nodal > 0:
                self.smax_norm = self.smax_nodal

        self._compute_metrics(m_data)
        return True

    def analyze_em(self):
        """EM-Analyse."""
        # EM-Pfad fuer LOC=A in Iteration 2 (Schicht 3 — Surface-Faces). Iteration 1 nur Gauss.
        if self.loc == "A":
            if not self._load_surface_properties():
                print("   [Analyst] EM-Surface-Daten fehlen. LOC=A im EM-Pfad benoetigt -A.out.")
                return False
            return self._analyze_em_surface()

        if not self.elements:
            if not self._load_mesh_properties(): return False

        if self.smax_norm <= 0:
            print("   [Analyst] Warnung: smax_norm <= 0. Nutze Nodal Smax.")
            self.smax_norm = self.smax_nodal

        m_data = {0: self.vtot_num}
        active_entries = []
        for vol, s1, s2, s3 in self.elements:
            ratios = []
            if self.smax_norm > 0:
                if self.crit == "S1":
                    # S1: Nur erste Hauptspannung
                    if s1 > 0: ratios.append(s1 / self.smax_norm)
                else:
                    # PIA: Alle positiven Hauptspannungen
                    if s1 > 0: ratios.append(s1 / self.smax_norm)
                    if s2 > 0: ratios.append(s2 / self.smax_norm)
                    if s3 > 0: ratios.append(s3 / self.smax_norm)
            if ratios:
                active_entries.append((vol, ratios))

        for m in range(self.m_start, self.m_max + 1):
            veff_sum = 0.0
            for vol, ratios in active_entries:
                term = sum(r**m for r in ratios)
                veff_sum += vol * term
            m_data[m] = veff_sum * self.symmetry_factor


        self._compute_metrics(m_data)
        return True

    def _load_surface_properties(self):
        """v12.1: Liest Surface-Face-Daten aus {apdl_name}-A.out (LOC=A, EM-Pfad).

        Datei-Schemata:
            v12.0 / NOD (Single-File, 6 Spalten):
                face_id  parent_elem  area  S1_face  S2_face  S3_face
            v12.1 / NOD (Single-File, 10 Spalten — mit Knoten-IDs fuer VTK):
                face_id  parent_elem  area  S1_face  S2_face  S3_face  n1 n2 n3 n4
            v12.1 / RAW (Zwei-Datei-Schema):
                {apdl_name}-A.out:   face_id  parent  area  S1=0  S2=0  S3=0  n1 n2 n3 n4
                {apdl_name}-A-S.out: PRESOL,S,PRIN gefiltert auf Boundary-Solids
                Stress-Berechnung: 4 Face-Eckknoten mit Parent-RAW-Werten gemittelt

        Befuellt:
            self.surface_faces = [(area, s1, s2, s3, [n1, n2, n3, n4]), ...]
            self.atot_num      = sum(area) * symmetry_factor
            self.smax_norm     = je nach NORM (NOD/EM/GP-Fallback)
        """
        a_path = os.path.join(self.tables_dir, f"{self.apdl_name}-A.out")
        a_s_path = os.path.join(self.tables_dir, f"{self.apdl_name}-A-S.out")

        if not os.path.exists(a_path):
            print(f"   [Analyst] CRITICAL: Surface-Datei fehlt: {a_path}")
            return False

        # smax_nodal aus -A.out-Header lesen (Surface-Pfad).
        # Bisher wurde nur {apdl_name}-V.out gesucht — diese Datei
        # existiert aber nur im Volumen-Pfad (LOC=V). Bei LOC=A blieb smax_nodal=0
        # und der Fallback in Schritt 4 unten setzte smax_nodal = max(s1_face),
        # was die NDX/EMX-Diskriminierung kollabieren liess (REF=NDX vs REF=EMX
        # lieferten bit-identische Aeff/Pf-Werte).
        # Fix: APDL-Macros post_EM_*_surf.inp schreiben jetzt einen smax-Header
        # in -A.out (analog zu -V.out). Wir lesen primaer von dort; -V.out
        # bleibt als Legacy-Fallback fuer alte Run-Verzeichnisse.
        if self.smax_nodal == 0.0:
            val = self._parse_header_smax(a_path)
            if val:
                self.smax_nodal = val
            else:
                # Legacy-Fallback fuer alte Runs vor v14.1 (ohne smax-Header in -A.out)
                v_path = os.path.join(self.tables_dir, f"{self.apdl_name}-V.out")
                if os.path.exists(v_path):
                    val = self._parse_header_smax(v_path)
                    if val:
                        self.smax_nodal = val

        # AVG=RAW + PRESOL-File vorhanden → Zwei-Datei-Schema
        # Hard-Fail bei STRESS=RAW + LOC=A wenn -A-S.out fehlt.
        # Vorher: silent fallback auf Legacy-Pfad mit 0.0-Platzhaltern aus
        # post_EM_RAW_surf.inp -> stille Datenkorruption, Aeff/Pf basierend auf Nullen.
        if self.avg == "RAW" and not os.path.exists(a_s_path):
            print(f"   [Analyst] CRITICAL: STRESS=RAW + LOC=A benoetigt {a_s_path}")
            print(f"   [Analyst]            APDL-Macro post_EM_RAW_surf.inp muss diese Datei erzeugen.")
            print(f"   [Analyst]            Abbruch (kein silent-Fallback auf 0.0-Platzhalter).")
            return False
        use_raw = (self.avg == "RAW")

        # 1) Geometrie aus -A.out lesen (Spalten: face_id, parent, area, S1, S2, S3, [n1..n4])
        face_geom = []  # [(parent, area, s1_embed, s2_embed, s3_embed, n_list_or_None), ...]
        try:
            with open(a_path, 'r', encoding='utf-8', errors='ignore') as f:
                start_reading = False
                for line in f:
                    if "face_id" in line and "area" in line:
                        start_reading = True
                        continue
                    if not start_reading:
                        continue
                    parts = line.split()
                    if len(parts) < 6:
                        continue
                    try:
                        parent = int(round(float(parts[1])))
                        area = float(parts[2])
                        s1 = float(parts[3])
                        s2 = float(parts[4])
                        s3 = float(parts[5])
                        # Knoten-IDs (10-Spalten-Format ab v12.1) — None wenn 6-Spalten Legacy
                        n_list = None
                        if len(parts) >= 10:
                            n_list = [int(round(float(parts[6]))), int(round(float(parts[7]))),
                                      int(round(float(parts[8]))), int(round(float(parts[9])))]
                        if area > 0:
                            face_geom.append((parent, area, s1, s2, s3, n_list))
                    except ValueError:
                        continue
        except Exception as e:
            print(f"   [Analyst] Fehler beim Lesen der Surface-Datei: {e}")
            return False

        if not face_geom:
            print(f"   [Analyst] CRITICAL: Keine Face-Daten in {a_path}")
            return False

        # 2) Stress-Mapping je nach AVG-Modus
        self.surface_faces = []

        if use_raw:
            # PRESOL parsen → {elem_id: [(node_id, s1, s2, s3), ...]}
            print(f"   [Analyst] EM-Surface RAW: parse PRESOL aus {os.path.basename(a_s_path)}")
            clean_csv = a_s_path.replace('-A-S.out', '-A-S-clean.csv')
            if not os.path.exists(clean_csv):
                parse_presol_to_csv(a_s_path, clean_csv)
            raw_data = read_stress_csv(clean_csv)

            # STRESS=RAW Konsistenz — smax_nodal aus PRESOL max(s1) lesen
            self._override_smax_nodal_from_raw_csv()

            n_skipped = 0
            for parent, area, _s1, _s2, _s3, n_list in face_geom:
                if parent not in raw_data or n_list is None:
                    n_skipped += 1
                    continue
                # Pro Face: 4 Eckknoten in Parent-Element-Knotenliste finden
                parent_nodes = raw_data[parent]  # list of (node_id, s1, s2, s3)
                node_to_stress = {nid: (s1, s2, s3) for nid, s1, s2, s3 in parent_nodes}
                stresses = [node_to_stress[n] for n in n_list if n in node_to_stress]
                if len(stresses) < 4:
                    n_skipped += 1
                    continue
                s1_face = sum(s[0] for s in stresses) / len(stresses)
                s2_face = sum(s[1] for s in stresses) / len(stresses)
                s3_face = sum(s[2] for s in stresses) / len(stresses)
                self.surface_faces.append((area, s1_face, s2_face, s3_face, n_list))

            if n_skipped > 0:
                # klarere Formulierung — User soll wissen, dass das normal ist
                # bei Surface-Pfaden: Knoten am Innenrand oder Symmetrie-Schnitt haben keine
                # PRESOL-Stresses, weil sie nicht auf der ausgewaehlten Surface liegen.
                print(f"   [Analyst] Info: {n_skipped} Surface-Faces ohne PRESOL-Stress-Daten "
                      f"uebersprungen (normal bei Innenrand-/Symmetrie-Knoten in RAW-Pfad).")
        else:
            # NOD-Modus oder Legacy: Stresses sind bereits in -A.out eingebettet
            for parent, area, s1, s2, s3, n_list in face_geom:
                self.surface_faces.append((area, s1, s2, s3, n_list))

        if not self.surface_faces:
            print(f"   [Analyst] CRITICAL: Keine Surface-Faces nach Stress-Mapping")
            return False

        # 3) Numerische Gesamtflaeche
        self.atot_num = sum(a for a, *_ in self.surface_faces) * self.symmetry_factor

        # 4) smax_norm — NORM-aware (NORM=EM via max(s1_face) entsperrt)
        max_s1 = max(s1 for _, s1, _, _, _ in self.surface_faces)
        norm_upper = self.norm.upper() if self.norm else "NOD"

        if norm_upper == "NOD" and self.smax_nodal > 0:
            self.smax_norm = self.smax_nodal
        elif norm_upper == "EM":
            # Element-Maximum-Aequivalent fuer Surface = max Face-Mittelwert (s1)
            self.smax_norm = max_s1
            print(f"   [Analyst] NORM=EM bei LOC=A: smax_norm = max(s1_face) = {max_s1:.4e}")
        else:
            # GP nicht unterstuetzt (im Run Manager blockiert), aber Fallback: Face-Maximum
            self.smax_norm = max_s1

        if self.smax_nodal == 0.0:
            # Sollte nicht mehr greifen, weil -A.out jetzt einen smax-Header hat.
            # Falls doch (korrupter Output / alter Run): warnen + auf max_s1 zurueckfallen.
            print(f"   [Analyst] WARNUNG: smax_nodal nicht aus -A.out lesbar — "
                  f"Fallback auf max(s1_face)={max_s1:.4e}. NDX-Normierung methodisch "
                  f"degradiert (gleich wie EMX).")
            self.smax_nodal = max_s1
        if self.sigma0 == 0.0:
            self.sigma0 = self._compute_sigma0(self.smax_nodal)

        return True

    def _analyze_em_surface(self):
        """v12.0+v12.1: EM-Analyse fuer LOC=A — summiert ueber Surface-Faces statt Volumen-Elemente.

        Aeff_S1  = sum_i A_i * (s1_i / smax_norm)^m
        Aeff_PIA = sum_i A_i * sum_j positive(s_j_i / smax_norm)^m

        v12.1: surface_faces enthaelt jetzt 5-Tupel (area, s1, s2, s3, n_list) — n_list
        wird hier ignoriert (nur fuer VTK-Surface-Export relevant).
        """
        if self.smax_norm <= 0:
            print("   [Analyst] Warnung: smax_norm <= 0. Nutze Nodal Smax.")
            self.smax_norm = self.smax_nodal

        m_data = {0: self.atot_num}
        active_entries = []
        for area, s1, s2, s3, _n_list in self.surface_faces:
            ratios = []
            if self.smax_norm > 0:
                if self.crit == "S1":
                    if s1 > 0: ratios.append(s1 / self.smax_norm)
                else:
                    # PIA: alle positiven Hauptspannungen
                    if s1 > 0: ratios.append(s1 / self.smax_norm)
                    if s2 > 0: ratios.append(s2 / self.smax_norm)
                    if s3 > 0: ratios.append(s3 / self.smax_norm)
            if ratios:
                active_entries.append((area, ratios))

        for m in range(self.m_start, self.m_max + 1):
            aeff_sum = 0.0
            for area, ratios in active_entries:
                term = sum(r**m for r in ratios)
                aeff_sum += area * term
            m_data[m] = aeff_sum * self.symmetry_factor

        self._compute_metrics(m_data)
        return True

    def _compute_metrics(self, m_data_num):
        """Berechnet Pf, Fehler und Fold Change fuer Veff und Pf."""
        self.results = []

        # Numerische Referenz: smax_norm (case-abhaengig: NOD/GP/EM)
        smax_calc = self.smax_norm
        sigma0_calc = self.sigma0

        # Analytische Referenz (smax_ref je Case, z.B. K_t=3 fuer PWH/PWHR)
        # Als Instanzvariable speichern fuer CSV-Export (S_ana)
        self.smax_ref = get_analytical_smax_ref(
            self.case_type, self.smax_nodal, self.load_n,
            dim_x=self.case_x, dim_y=self.case_y, dim_z=self.case_z
        )
        smax_ref = self.smax_ref
        sigma0_ref = self.sigma0

        # LOC-aware Dispatch — bei LOC=A wird Aeff statt Veff berechnet
        if self.loc == "A":
            ana_func = self._get_analytic_aeff
            eff_ana_key = "Aeff_ana"
            eff_num_key = "Aeff_num"
            err_rel_key = "Err_A_rel"
            err_fold_key = "Err_A_Fold"
        else:
            ana_func = self._get_analytic_veff
            eff_ana_key = "Veff_ana"
            eff_num_key = "Veff_num"
            err_rel_key = "Err_V_rel"
            err_fold_key = "Err_V_Fold"

        for m in range(self.m_start, self.m_max + 1):
            if m not in m_data_num: continue

            v_num = m_data_num[m]
            v_ana = ana_func(m)

            # --- Pf Berechnung (Formel strukturell identisch fuer V und A) ---
            # pf_from_hazard() rechnet im Log-Raum mit Clamping
            # gegen Auslöschung (kleines H) und Overflow (grosses H, Power-Term).
            # Mathematisch identisch zu 1 - exp(-Veff * (smax/sigma0)^m), aber
            # numerisch robust bei Extremwerten.
            pf_ana = pf_from_hazard(v_ana, smax_ref,  sigma0_ref,  m)
            pf_num = pf_from_hazard(v_num, smax_calc, sigma0_calc, m)

            # --- Fehler Effective Quantity (Veff oder Aeff je nach LOC) ---
            err_v_rel = (v_num - v_ana) / v_ana if v_ana != 0 else 0.0

            err_v_fold = 0.0
            if v_ana != 0 and v_num != 0:
                ratio = v_num / v_ana
                if ratio >= 1.0:
                    err_v_fold = ratio
                else:
                    err_v_fold = -1.0 / ratio
            elif v_ana == 0 and v_num == 0:
                err_v_fold = 1.0
            else:
                err_v_fold = 0.0

            # --- Fehler Pf (LOC-unabhaengig) ---
            err_pf_rel = (pf_num - pf_ana) / pf_ana if pf_ana != 0 else 0.0
            err_pf_abs = pf_num - pf_ana

            err_pf_fold = 0.0
            if pf_ana != 0 and pf_num != 0:
                ratio_pf = pf_num / pf_ana
                if ratio_pf >= 1.0:
                    err_pf_fold = ratio_pf
                else:
                    err_pf_fold = -1.0 / ratio_pf
            elif pf_ana == 0 and pf_num == 0:
                err_pf_fold = 1.0
            else:
                err_pf_fold = 0.0

            row = {
                "m": m,
                eff_ana_key: v_ana,
                eff_num_key: v_num,
                err_rel_key: err_v_rel,
                err_fold_key: err_v_fold,
                "Pf_ana": pf_ana,
                "Pf_num": pf_num,
                "Err_Pf_rel": err_pf_rel,
                "Err_Pf_Fold": err_pf_fold,
                "Err_Pf_abs": err_pf_abs,
                "S_num": smax_calc,
                "S_ana": smax_ref,
            }

            self.results.append(row)

    def _build_metadata(self):
        """Sammelt alle Metadaten fuer JSON/Report-Export."""
        from datetime import datetime

        self.run_metadata = {
            "case_id": self.case_name,
            "timestamp": datetime.now().isoformat(),
            "case_type": self.case_type,
            # Top-Level-Felder fuer G26-Reference-Disambiguierung.
            # Werden vom G26-Generator + calc_differentiated_errors._merge_dual_g26_refs
            # gelesen, um Reference-Zeilen pro (LOAD_N, SIG_0_ABS, SIG_0)-Tupel zu trennen.
            "load_n": self.load_n,
            "sig_0_abs": self.sig_0_abs,
            "sig_0": self.sig_0_value,
            "sigma0": self.sigma0,
            # pia_fix als Top-Level-Identitaetsmerkmal. G26-Generator akzeptiert
            # in v20.2 nur Cases mit pia_fix=True (PIA-at-GP korrekt). Legacy-Cases
            # (pia_fix=False) bekommen keine G26-Reference -> ERR_EXT-Decomposition
            # blockiert mit klarer Fehlermeldung.
            "pia_fix": self.pia_fix,
            "simulation": {
                "ansys_exe_path": self.ansys_exe_path,
                "num_cores": self.num_cores,
                "stress": self.stress,   # Methoden-Baukasten Stress|Int|Ref
                "int_type": self.int_type,
                "ref": self.ref,
                "crit": self.crit,
                "loc": self.loc,
                "sig_0_abs": self.sig_0_abs,
                "sig_0_value": self.sig_0_value,
            },
            "geometry": {
                "case_x": self.case_x,
                "case_y": self.case_y,
                "case_z": self.case_z,
                "vtot_exact": self.vtot_exact,
                "vtot_num": self.vtot_num,
                "atot_exact": getattr(self, "atot_exact", None),
                "atot_num": getattr(self, "atot_num", 0.0),
                "symmetry_factor": self.symmetry_factor,
            },
            "mesh": {
                "mesh_x": self.mesh_x,
                "mesh_y": self.mesh_y,
                "mesh_z": self.mesh_z,
                "num_elements": len(self.elements) if self.elements else getattr(self, 'num_elements_gauss', 0),
            },
            "load": {
                "load_n": self.load_n,
            },
            "stress_reference": {
                "smax_nodal": self.smax_nodal,
                "smax_norm": self.smax_norm,
                "smax_ref": getattr(self, 'smax_ref', 0.0),
                "ref_type": self.ref,   # REF-Token (EMX/GPX/NDX)
                "sigma0": self.sigma0,
            },
            "results_summary": {
                "m_range": [self.m_start, self.m_max],
                "num_m_values": len(self.results),
            }
        }

    def save_csv(self):
        out_file = os.path.join(self.tables_dir, f"ANALYSIS_{self.case_name}.csv")
        if not self.results: return

        # domain-spezifische Spaltennamen (Veff_* fuer LOC=V, Aeff_* fuer LOC=A)
        if self.loc == "A":
            cols = ["m", "Aeff_ana", "Aeff_num", "Err_A_rel", "Err_A_Fold",
                    "Pf_ana", "Pf_num", "Err_Pf_rel", "Err_Pf_Fold", "Err_Pf_abs",
                    "S_num", "S_ana"]
        else:
            cols = ["m", "Veff_ana", "Veff_num", "Err_V_rel", "Err_V_Fold",
                    "Pf_ana", "Pf_num", "Err_Pf_rel", "Err_Pf_Fold", "Err_Pf_abs",
                    "S_num", "S_ana"]

        try:
            with open(out_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=cols, delimiter=';')
                writer.writeheader()
                writer.writerows(self.results)
            print(f"   [Analyst] Gespeichert: {os.path.basename(out_file)}")
        except Exception as e:
            print(f"   [Analyst] Fehler beim Speichern: {e}")

def run_analysis_for_case(run_dir, full_name, int_type, mesh_y,
                          stress="AVG", ref="EMX", sig_0_abs=True, sig_0=100.0,
                          case_x=0, case_y=0, case_z=0, load_n=1.0,
                          ansys_exe_path="", num_cores=0, apdl_name=None,
                          mesh_x=0, mesh_z=1, crit="PIA", loc="V", export_vtk=False,
                          err_ext=None, pia_fix=False):
    """v14: Externe API mit STRESS/REF (Methoden-Baukasten Stress|Int|Ref).

    v19.0: pia_fix steuert PIA-at-Gauss-Point Korrektur (silent ignore bei INT=EM).
    """
    analyst = ResultAnalyst(
        run_dir, full_name, mesh_y,
        case_x=case_x, case_y=case_y, case_z=case_z,
        stress=stress, ref=ref, sig_0_abs=sig_0_abs, sig_0=sig_0, load_n=load_n,
        ansys_exe_path=ansys_exe_path, num_cores=num_cores,
        apdl_name=apdl_name,
        int_type=int_type, mesh_x=mesh_x, mesh_z=mesh_z,
        crit=crit, loc=loc, pia_fix=pia_fix
    )

    success = False
    if int_type.startswith("G"):
        success = analyst.analyze_gauss()
    elif int_type == "EM":
        success = analyst.analyze_em()

    if success:
        analyst._build_metadata()
        analyst.save_csv()

        # JSON-Metadata (nur Metadaten, ohne Results-Array)
        json_path = os.path.join(analyst.tables_dir, f"{analyst.case_name}.json")
        save_json_metadata(analyst.run_metadata, json_path)

        # Report-Generierung ans ENDE der Pipeline verschoben (siehe unten),
        # damit _extended.csv (v15-Decomposition) bereits existiert und vom
        # Auto-Detection-Block in generate_markdown_report() gefunden wird.

        # --- VTK Export (optional, v10.0 + v12.1 Surface) ---
        # Gauss-Pfad hat keine pro-Element/Face-Aufschluesselung der Stresses
        # (Fortran liefert nur integrierte Veff/Seff pro m-Wert in VEFF_*.out).
        # Daher ist VTK-RRI-Export fuer Gauss-Cases methodisch nicht sinnvoll
        # und wird komplett uebersprungen — sowohl LOC=V (Hex-RRI) als auch
        # LOC=A (Surface-Quad-RRI). VTK ist dem EM-Pfad vorbehalten.
        _vtk_skip_gauss = export_vtk and int_type and int_type.upper().startswith("G")
        if _vtk_skip_gauss:
            print(f"   [VTK] uebersprungen: VTK-Export fuer Gauss-Pfad nicht vorgesehen "
                  f"(keine pro-Element-Spannungen verfuegbar). VTK ist EM-only.")

        if export_vtk and not _vtk_skip_gauss:
            try:
                vtk_dir = os.path.join(analyst.run_dir, "vtk")

                # LOC=A → Surface-VTK (Quad-Cells) + Volume-Context (graue Hex-Cells)
                if analyst.loc == "A":
                    from vtk_exporter import (export_vtk_volume_context,
                                              export_vtk_surface_series)
                    # 1) Volumen-Kontext (m-unabhaengig, einmalig)
                    export_vtk_volume_context(
                        tables_dir=analyst.tables_dir,
                        vtk_dir=vtk_dir,
                        case_name=full_name,
                        apdl_name=analyst.apdl_name,
                        int_type=int_type,
                    )
                    # 2) Surface-RRI-Series (pro m-Wert, Quad-Cells)
                    if hasattr(analyst, 'surface_faces') and analyst.surface_faces:
                        export_vtk_surface_series(
                            tables_dir=analyst.tables_dir,
                            vtk_dir=vtk_dir,
                            case_name=full_name,
                            apdl_name=analyst.apdl_name,
                            int_type=int_type,
                            surface_faces=analyst.surface_faces,
                            smax_norm=analyst.smax_norm,
                            sigma0=analyst.sigma0,
                            crit=analyst.crit,
                        )
                    else:
                        print("   [VTK-Surface] WARNUNG: surface_faces leer, Surface-Export uebersprungen.")
                else:
                    # LOC=V (Standard): Volumen-RRI-Series (Hex-Cells)
                    from vtk_exporter import export_vtk_series
                    export_vtk_series(
                        tables_dir=analyst.tables_dir,
                        vtk_dir=vtk_dir,
                        case_name=full_name,
                        apdl_name=analyst.apdl_name,
                        int_type=int_type,
                        avg=analyst.avg,
                        crit=analyst.crit,
                        smax_norm=analyst.smax_norm,
                        sigma0=analyst.sigma0,
                        symmetry_factor=analyst.symmetry_factor,
                        elements_em=analyst.elements if int_type == "EM" else None,
                        element_ids_em=analyst.element_ids if int_type == "EM" else None,
                        case_type=analyst.case_type,
                        case_x=analyst.case_x,
                        case_y=analyst.case_y,
                        case_z=analyst.case_z,
                        load_n=analyst.load_n,
                    )
            except ImportError:
                print("   [Analyst] WARNUNG: vtk_exporter nicht verfuegbar.")
            except Exception as e:
                print(f"   [Analyst] VTK-Export fehlgeschlagen: {e}")

        # --- Fehlerzerlegung (v15: 5-Klassen Decomposition; v16: 3-Stufen ERR_EXT-Logik) ---
        # Drei-Stufen-Klassifikation der ERR_EXT-Behandlung
        #   ERR_EXT=0 (explizit): clean SKIP mit Info-Log
        #   ERR_EXT leer/Default: Decomposition optional, Warning bei fehlender Reference
        #   ERR_EXT=1 (explizit) + Reference fehlt: Case als DEGRADED markiert (.degraded-File + JSON-Flag)
        err_ext_explicit_zero = (err_ext is False)  # User hat explizit "0" gesetzt
        err_ext_explicit_one = (err_ext is True)    # User hat explizit "1" gesetzt
        _do_err_ext = err_ext if err_ext is not None else INTEGRATE_ERROR_DIFFERENTIATION

        if err_ext_explicit_zero:
            print(f"   [ErrDecomp] INFO: ERR_EXT=0 — Fehlerzerlegung uebersprungen (Basis-CSV genuegt fuer diesen Case).")

        case_degraded = False
        degraded_reason = None

        if _do_err_ext and REFERENCE_CSV:
            try:
                from calc_differentiated_errors import compute_v15_errors, _merge_dual_g26_refs
                import pandas as pd

                analysis_csv = os.path.join(analyst.tables_dir, f"ANALYSIS_{full_name}.csv")
                df_analysis = pd.read_csv(analysis_csv, delimiter=';')

                # Merge-Keys aus full_name via zentralen Parser
                info = mn.parse_case_id(full_name)
                if "error" in info:
                    raise ValueError(f"Konnte full_name nicht parsen: {full_name} ({info['error']})")

                if not os.path.exists(REFERENCE_CSV):
                    # 3-Stufen-Verhalten je nach ERR_EXT-Flag-Status
                    if err_ext_explicit_one:
                        # User hat ERR_EXT=1 explizit angefordert — Case wird als DEGRADED markiert
                        case_degraded = True
                        degraded_reason = (f"ERR_EXT=1 explizit angefordert, aber "
                                           f"Reference_Gauss26.csv fehlt unter {REFERENCE_CSV}")
                        print(f"!!! [ErrDecomp] DEGRADED: ERR_EXT=1 angefordert, aber "
                              f"Reference_Gauss26.csv fehlt — Case als unvollstaendig markiert.")
                        print(f"    [ErrDecomp]            Basis-CSV ist geschrieben, _extended.csv FEHLT.")
                        print(f"    [ErrDecomp]            Aktion: G26-Runs (RAW+AVG) + "
                              f"'06-generate_gauss26_reference.py' und Case neu rechnen.")
                        # Marker-Datei im Case-Verzeichnis
                        degraded_path = os.path.join(analyst.run_dir, "_DEGRADED.txt")
                        with open(degraded_path, 'w', encoding='utf-8') as df:
                            df.write(f"DEGRADED Case: {full_name}\n")
                            df.write(f"Reason: {degraded_reason}\n")
                            df.write(f"Status: Basis-CSV geschrieben, _extended.csv FEHLT.\n")
                            df.write(f"Aktion: G26-Reference erzeugen + Case re-runnen.\n")
                    else:
                        # ERR_EXT leer/Default — Warning + Basis-CSV reicht
                        print(f"   [ErrDecomp] WARN: Reference_Gauss26.csv fehlt — "
                              f"Fehlerzerlegung uebersprungen, Basis-CSV reicht (ERR_EXT war optional).")
                        print(f"   [ErrDecomp]       Erwartet: {REFERENCE_CSV}")
                else:
                    df_ref = pd.read_csv(REFERENCE_CSV, delimiter=';')

                    # Dual-Merge mit RAW + AVG-G26-Referenzen
                    # LOAD_N/SIG_0_ABS/SIG_0 aus analyst-Instance ueber-
                    # geben, damit G26-Merge nur passende Reference-Zeilen findet.
                    # pia_fix als Hard-Fail-Bedingung uebergeben.
                    df_merged = _merge_dual_g26_refs(
                        df_analysis, df_ref, info,
                        analyst.load_n, analyst.sig_0_abs, analyst.sig_0_value,
                        pia_fix=analyst.pia_fix
                    )

                    # Fehlerzerlegung berechnen (kann RuntimeError werfen bei fehlender G26-Ref)
                    df_result = compute_v15_errors(df_merged, info["stress"], info["int"])

                    # Hilfsspalten entfernen
                    drop_cols = ['Mesh', 'Loadcase_self', 'Loadcase_other']
                    df_result = df_result.drop(columns=[c for c in drop_cols if c in df_result.columns])

                    out_path = analysis_csv.replace('.csv', '_extended.csv')
                    df_result.to_csv(out_path, sep=';', index=False)
                    print(f"   [Analyst] v15-Fehlerzerlegung: {os.path.basename(out_path)}")
            except RuntimeError as e:
                # klare Fehlermeldung bei fehlender G26-Referenz
                # bei explizitem ERR_EXT=1 als DEGRADED markieren
                print(f"   [Analyst] ERR_EXT v15 abgebrochen: {e}")
                if err_ext_explicit_one:
                    case_degraded = True
                    degraded_reason = f"ERR_EXT=1 explizit, aber Fehlerzerlegung fehlgeschlagen: {e}"
                    degraded_path = os.path.join(analyst.run_dir, "_DEGRADED.txt")
                    with open(degraded_path, 'w', encoding='utf-8') as df:
                        df.write(f"DEGRADED Case: {full_name}\n")
                        df.write(f"Reason: {degraded_reason}\n")
                        df.write(f"Status: Basis-CSV geschrieben, _extended.csv FEHLT.\n")
            except ImportError as e:
                print(f"   [Analyst] WARNUNG: calc_differentiated_errors nicht verfuegbar: {e}")
            except Exception as e:
                print(f"   [Analyst] Fehlerzerlegung fehlgeschlagen: {e}")

        # Wenn Case als DEGRADED markiert wurde, JSON-Metadata nachtraeglich aktualisieren
        if case_degraded:
            try:
                import json
                json_path = os.path.join(analyst.tables_dir, f"{analyst.case_name}.json")
                if os.path.exists(json_path):
                    with open(json_path, 'r', encoding='utf-8') as f:
                        meta = json.load(f)
                    meta["status"] = "DEGRADED"
                    meta["degraded_reason"] = degraded_reason
                    meta["err_ext_requested"] = True
                    meta["extended_csv_written"] = False
                    with open(json_path, 'w', encoding='utf-8') as f:
                        json.dump(meta, f, indent=2, ensure_ascii=False)
                    print(f"   [Analyst] JSON-Metadata mit DEGRADED-Status aktualisiert.")
            except Exception as e:
                print(f"   [Analyst] WARN: JSON-Update fehlgeschlagen: {e}")

        # Markdown-Report ans ENDE der Pipeline verschoben.
        # Vorher: Report wurde direkt nach JSON-Metadata erzeugt — bevor
        # _extended.csv vom v15-Decomposition-Block geschrieben war. Folge:
        # _find_extended_csv() in report_generator fand die Datei nicht und der
        # Fehlerauswertungs-Abschnitt zeigte "Keine erweiterte Fehlerzerlegung
        # vorhanden", obwohl sie 2 Schritte spaeter erzeugt wurde.
        # Jetzt: Report wird nach VTK-Export + ERR_EXT + DEGRADED-Update generiert.
        image_paths = {
            's1': f"../plots/{analyst.apdl_name}-s1.png",
            's2': f"../plots/{analyst.apdl_name}-s2.png",
        }
        report_path = os.path.join(analyst.tables_dir, f"{analyst.case_name}_REPORT.md")
        generate_markdown_report(analyst.run_metadata, analyst.results, image_paths, report_path)

    else:
        print("   [Analyst] Analyse fehlgeschlagen.")
