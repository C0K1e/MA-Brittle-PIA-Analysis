"""
parse_presol_output.py
======================
Standalone-Modul zum Parsen von ANSYS MAPDL PRESOL-Ausgabedateien (-S.out).

Aufgabe:
  - Rohe MAPDL-PRESOL-Datei einlesen und bereinigen
  - Elementbloecke mit Knotenspannungen S1, S2, S3 extrahieren
  - Bereinigte Daten als flache CSV ausgeben: element, node, s1, s2, s3

Oeffentliche API:
  parse_presol(filepath)          -> {elem_id: [(node_id, s1, s2, s3), ...]}
  write_stress_csv(data, outpath) -> None
  parse_presol_to_csv(in, out)    -> {elem_id: [(node_id, s1, s2, s3), ...]}
  read_stress_csv(csv_path)       -> {elem_id: [(node_id, s1, s2, s3), ...]}
  compute_element_means(data)     -> {elem_id: (s1_mean, s2_mean, s3_mean)}

  V9 Gauss-RAW Hilfsfunktionen:
  read_element_file_with_ids(fp)  -> (elem_order, node_order)
  write_fortran_stress_file(...)  -> None  (schreibt effVol_RawStress.out)

CLI-Aufruf:
  python parse_presol_output.py path/to/{name}-S.out
  -> erzeugt: path/to/{name}-S-clean.csv

Robustheit (v8.5 — 4-Schichten-Filter):
  Schicht 1 — SKIP_KEYWORDS:  bekannte Banner-/Header-Zeilen sofort verwerfen
  Schicht 2 — ALPHA_RE:       Zeilen mit Buchstaben (ausser E/e) im Datenbereich
                               sofort verwerfen; fangt Modellnamen (z.B. PWH-200x200x1-...)
                               und unbekannte Header ab
  Schicht 3 — Node-Validierung: node muss 1 <= node <= MAX_NODE_ID sein;
                               verhindert negative IDs aus Modellnamen (z.B. "-200")
  Schicht 4 — Mindestfelder:  mindestens 4 numerische Felder (NODE S1 S2 S3)
"""

import os
import re
import csv
import sys


# ---------------------------------------------------------------------------
# Konstanten fuer den Parser
# ---------------------------------------------------------------------------

# Regex: findet einzelne Float-Zahlen (inkl. konkatenierter negativer Werte)
FLOAT_RE = re.compile(r'[+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?')

# Regex: erkennt "ELEMENT = 42" oder "ELEMENT=42" (case-insensitive)
ELEM_RE = re.compile(r'ELEMENT\s*=\s*(\d+)', re.IGNORECASE)

# Regex: findet alphabetische Zeichen AUSSER E/e (wissenschaftliche Notation).
# Echte PRESOL-Datenzeilen enthalten ausschliesslich Ziffern, Vorzeichen, Punkte
# und E/e in Exponenten.  Jeder andere Buchstabe zeigt eine Header-/Titelzeile an.
#   Muster: [a-d f-z] (alle Kleinbuchstaben ausser 'e')
#         + [A-D F-Z] (alle Grossbuchstaben ausser 'E')
ALPHA_RE = re.compile(r'[a-df-zA-DF-Z]')

# Schluesselwoerter in MAPDL-Bannern und Seitenkoepfen -> Zeile sofort ueberspringen.
# "***" markiert zusaetzlich den Start einer neuen MAPDL-Seite und setzt in_data zurueck.
SKIP_KEYWORDS = (
    '***',            # MAPDL-Seitentrenner  (setzt in_data zurueck)
    'POST1 ELEMENT',
    'LOAD STEP',
    'DEGREE OF FREEDOM',
    'THIS PRINTOUT',
    'THE FOLLOWING',
    'LOAD CASE',
    'ANSYS MAPDL',
    'MAXIMUM ABSOLUTE',
    'NOTE',
    'WARNING',
    'ANSYS ACADEMIC',
    'PRINT ELEMENT',
    'TIME=',
    'SUBSTEP=',
    'VERSION=',       # z.B. "01055371  VERSION=WINDOWS x64 ..." (Lizenzzeile)
    'Ansys',          # z.B. "Ansys Mechanical Enterprise Academic Student"
    'RELEASE',        # z.B. "RELEASE 2025 R2" in MAPDL-Banner
    'COPYRIGHT',      # Lizenz-/Copyright-Bannerzeilen
    'ELEMENT NODAL',  # "***** POST1 ELEMENT NODAL STRESS LISTING *****"
    'SOLUTION PER',   # "PRINT S PRIN ELEMENT SOLUTION PER ELEMENT"
)

# Obergrenze fuer plausible Knoten-IDs.
# Selbst sehr grosse Netze (z.B. 200x200x1 mit SOLID186, 20 Knoten/Elem)
# erzeugen keine Node-IDs jenseits dieser Grenze.
MAX_NODE_ID = 10_000_000

# CSV-Ausgabe-Spalten
CSV_COLUMNS = ('element', 'node', 's1', 's2', 's3')


def _format_stress(v):
    """v20.0 (Bug #5): Formatiert Spannungswert mit 12 signifikanten Stellen in
    Wissenschafts-Notation. Vorher (v14.2): 5 sig figs + min. 3 Nachkommastellen.

    Motivation: ANSYS exportiert RAW-Spannungen mit `/FORMAT, 9, G, 20, 8`
    (9 sig figs). Die alte 5-sig-Rundung in dieser Datei verwarf 4 sig figs
    Praezision, was bei hohem Weibull-Modul m sichtbare Drift in Veff/Pf
    erzeugen konnte. Die bereinigte CSV ist Quelle fuer:
      - RAW-smax (max(s1) fuer smax_nodal-Override v14.2)
      - EM-Mittelwerte (compute_element_means)
      - effVol_RawStress.out (Fortran-Eingabe im RAW-Modus)

    Mit 12 sig figs in Wissenschafts-Notation ist die CSV stabil ueber alle
    Wertebereiche (sehr klein bis sehr gross), eindeutig parsebar und
    behaelt mehr Praezision als ANSYS aktuell liefert.

    Beispiele:
        911.50518   -> "9.115051800000E+02"   (vorher: "911.505")
        0.038695    -> "3.869500000000E-02"   (vorher: "0.038695")
        0.0         -> "0.000000000000E+00"
        1.5e-08     -> "1.500000000000E-08"

    Parameters
    ----------
    v : float
        Spannungswert.

    Returns
    -------
    str
        Formatierter String mit 12 sig figs in Wissenschafts-Notation.
    """
    return f"{v:.12E}"


# ---------------------------------------------------------------------------
# Kern-Parser
# ---------------------------------------------------------------------------

def parse_presol(filepath):
    """
    Liest eine rohe MAPDL-PRESOL-Datei (-S.out) und extrahiert
    elementweise Knotenspannungen.

    PRESOL-Spaltenreihenfolge:  NODE  S1  S2  S3  SINT  SEQV
    -> SINT und SEQV werden ignoriert.

    Parameters
    ----------
    filepath : str
        Pfad zur rohen PRESOL-Datei (-S.out).

    Returns
    -------
    dict
        {elem_id (int): [(node_id (int), s1 (float), s2 (float), s3 (float)), ...]}
        Leer bei Fehler oder leerer Datei.

    Notes
    -----
    4-Schichten-Filter (v8.5):
      1. SKIP_KEYWORDS   — bekannte Banner-/Kopfzeilen sofort verwerfen;
                           "***" setzt zusaetzlich in_data=False (Seitengrenze)
      2. ALPHA_RE        — Zeilen mit Nicht-E/e-Buchstaben verwerfen (Modellnamen,
                           unbekannte Header); greift nur im aktiven Datenbereich
      3. Node-Validierung — node muss 1 <= node <= MAX_NODE_ID (kein negatives Artefakt
                           aus Modellnamen, kein unplausibler Headerausdruck)
      4. Mindestfelder   — mindestens 4 Zahlen (NODE S1 S2 S3) erwartet
    """
    elem_data  = {}   # {elem_id: [(node_id, s1, s2, s3), ...]}
    current_id = None
    in_data    = False

    # Diagnosezaehler
    n_lines_total    = 0
    n_skipped_kw     = 0   # Schicht 1: durch SKIP_KEYWORDS verworfen
    n_rejected_alpha = 0   # Schicht 2: Buchstaben ausser E/e gefunden
    n_rejected_node  = 0   # Schicht 3: Node-ID ausserhalb 1..MAX_NODE_ID
    n_rejected_count = 0   # Schicht 4: weniger als 4 numerische Felder
    n_rejected_parse = 0   # Parse-Fehler (ValueError / IndexError)

    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                n_lines_total += 1
                stripped = line.strip()

                # Leerzeilen ueberspringen
                if not stripped:
                    continue

                # Schicht 1: MAPDL-Banner und Seitenkopfzeilen ueberspringen.
                # "***" markiert einen neuen MAPDL-Seitenblock -> in_data zuruecksetzen,
                # damit eventuelle Folgezeilen des Banners nicht als Daten landen.
                if any(kw in stripped for kw in SKIP_KEYWORDS):
                    n_skipped_kw += 1
                    if '***' in stripped:
                        in_data = False
                    continue

                # ELEMENT= Zeile -> neuer Elementblock beginnt
                m_elem = ELEM_RE.search(stripped)
                if m_elem:
                    current_id = int(m_elem.group(1))
                    in_data    = False
                    if current_id not in elem_data:
                        elem_data[current_id] = []
                    continue

                # Spalten-Header-Zeile erkennen: enthaelt sowohl "NODE" als auch "S1"
                if current_id is not None and 'S1' in stripped and 'NODE' in stripped:
                    in_data = True
                    continue

                # Datenzeilen: nur wenn wir im aktiven Datenbereich eines Elements sind
                if current_id is not None and in_data:

                    # Schicht 2: Alphabetische Zeichen ausser E/e -> kein Datensatz.
                    # Echte PRESOL-Zeilen enthalten nur Ziffern, Vorzeichen, Punkte
                    # und E/e (Exponent). Alles andere ist Header oder Titelzeile.
                    # Beispiele, die hier gefangen werden:
                    #   "PWH-200x200x1-EM-EM-RAW-4x4"  (Modellname, 'P','W','H','x',...)
                    #   "01055371  VERSION=WINDOWS x64" (falls SKIP_KW fehlt: 'V','R',...)
                    if ALPHA_RE.search(stripped):
                        n_rejected_alpha += 1
                        continue

                    # Schicht 4: Zahlen extrahieren, Mindestanzahl pruefen
                    nums = FLOAT_RE.findall(stripped)
                    if len(nums) < 4:
                        n_rejected_count += 1
                        continue

                    # Parse-Versuch: nums[0]=NODE, nums[1]=S1, nums[2]=S2, nums[3]=S3
                    try:
                        node_id = int(float(nums[0]))
                        s1      = float(nums[1])
                        s2      = float(nums[2])
                        s3      = float(nums[3])
                    except (ValueError, IndexError):
                        n_rejected_parse += 1
                        continue

                    # Schicht 3: Knoten-ID-Validierung.
                    # Negative IDs entstehen z.B. aus "PWH-200x200x1-..." (node=-200).
                    # Sehr grosse IDs entstehen z.B. aus "01055371 VERSION=..." (1055371
                    # waere noch plausibel; der SKIP_KW-Check greift vorher).
                    if not (1 <= node_id <= MAX_NODE_ID):
                        n_rejected_node += 1
                        continue

                    elem_data[current_id].append((node_id, s1, s2, s3))

    except FileNotFoundError:
        print(f"[parse_presol] FEHLER: Datei nicht gefunden: {filepath}")
        return {}
    except Exception as e:
        print(f"[parse_presol] FEHLER beim Lesen: {e}")
        return {}

    n_elem  = len(elem_data)
    n_nodes = sum(len(v) for v in elem_data.values())

    if n_elem == 0:
        print(f"[parse_presol] WARNUNG: Keine Elementbloecke gefunden in "
              f"{os.path.basename(filepath)}")
    else:
        print(f"[parse_presol] {n_elem} Elemente, {n_nodes} Knotenzeilen akzeptiert "
              f"aus {os.path.basename(filepath)}")

    # Diagnose-Ausgabe (immer, damit Fehlinterpretationen schnell auffallen)
    total_rejected = n_rejected_alpha + n_rejected_node + n_rejected_count + n_rejected_parse
    print(f"[parse_presol] Diagnose | Zeilen gesamt: {n_lines_total} | "
          f"KW-Skip: {n_skipped_kw} | "
          f"Alpha-Verwurf: {n_rejected_alpha} | "
          f"Node-Verwurf: {n_rejected_node} | "
          f"Feld-Verwurf: {n_rejected_count} | "
          f"Parse-Fehler: {n_rejected_parse} | "
          f"Verwurfe gesamt: {total_rejected}")
    if n_rejected_alpha > 0 or n_rejected_node > 0:
        print(f"[parse_presol] HINWEIS: {n_rejected_alpha} Zeilen durch Alpha-Filter, "
              f"{n_rejected_node} durch Node-Validierung verworfen — "
              f"pruefe ob die .out-Datei unerwartete Headerzeilen enthaelt.")

    return elem_data


# ---------------------------------------------------------------------------
# CSV-Ausgabe
# ---------------------------------------------------------------------------

def write_stress_csv(data, outpath):
    """
    Schreibt die PRESOL-Daten als flache, saubere CSV-Datei.

    Spalten: element, node, s1, s2, s3
    Sortiert: aufsteigend nach element, dann nach node.

    Parameters
    ----------
    data : dict
        {elem_id: [(node_id, s1, s2, s3), ...]}
    outpath : str
        Zielpfad fuer die CSV-Datei (z.B. .../tables/{name}-S-clean.csv).
    """
    try:
        with open(outpath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(CSV_COLUMNS)

            for elem_id in sorted(data.keys()):
                nodes = sorted(data[elem_id], key=lambda t: t[0])  # sort by node_id
                for node_id, s1, s2, s3 in nodes:
                    # einheitliches Stress-Format (5 sig figs, min. 3 Nachkomma)
                    writer.writerow([elem_id, node_id,
                                     _format_stress(s1),
                                     _format_stress(s2),
                                     _format_stress(s3)])

        n_rows = sum(len(v) for v in data.values())
        print(f"[parse_presol] CSV geschrieben: {os.path.basename(outpath)} ({n_rows} Zeilen)")

    except Exception as e:
        print(f"[parse_presol] FEHLER beim Schreiben der CSV: {e}")


# ---------------------------------------------------------------------------
# Convenience-Funktion
# ---------------------------------------------------------------------------

def parse_presol_to_csv(in_path, out_path):
    """
    Liest rohe PRESOL-Datei und schreibt bereinigte CSV in einem Schritt.

    Parameters
    ----------
    in_path : str
        Pfad zur rohen -S.out Datei.
    out_path : str
        Zielpfad fuer die bereinigte -S-clean.csv.

    Returns
    -------
    dict
        {elem_id: [(node_id, s1, s2, s3), ...]} — wie parse_presol().
    """
    data = parse_presol(in_path)
    if data:
        write_stress_csv(data, out_path)
    return data


# ---------------------------------------------------------------------------
# CSV-Einlesen (fuer Wiederverwendung ohne erneutes PRESOL-Parsen)
# ---------------------------------------------------------------------------

def read_stress_csv(csv_path):
    """
    Liest eine bereinigte -S-clean.csv zurueck in die interne dict-Struktur.

    Parameters
    ----------
    csv_path : str
        Pfad zur bereinigten CSV (Spalten: element, node, s1, s2, s3).

    Returns
    -------
    dict
        {elem_id (int): [(node_id (int), s1 (float), s2 (float), s3 (float)), ...]}
        Leer bei Fehler.
    """
    data = {}
    try:
        with open(csv_path, 'r', encoding='utf-8', errors='ignore') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    elem_id = int(row['element'])
                    node_id = int(row['node'])
                    s1      = float(row['s1'])
                    s2      = float(row['s2'])
                    s3      = float(row['s3'])
                    if elem_id not in data:
                        data[elem_id] = []
                    data[elem_id].append((node_id, s1, s2, s3))
                except (ValueError, KeyError):
                    continue
    except FileNotFoundError:
        print(f"[parse_presol] FEHLER: CSV nicht gefunden: {csv_path}")
        return {}
    except Exception as e:
        print(f"[parse_presol] FEHLER beim Lesen der CSV: {e}")
        return {}

    if not data:
        print(f"[parse_presol] WARNUNG: Keine Daten in CSV: {os.path.basename(csv_path)}")
    return data


# ---------------------------------------------------------------------------
# RAW-konsistenter smax_nodal aus -S-clean.csv
# ---------------------------------------------------------------------------

def max_s1_from_clean_csv(csv_path):
    """v14.2: Liest max(s1) aus einer bereinigten PRESOL-CSV.

    Bei STRESS=RAW soll auch der Normierungs-Maxwert (smax_nodal) aus dem
    RAW-Spannungsfeld kommen und nicht aus dem nodal-gemittelten APDL-Header.
    PRESOL liefert pro Element pro Eckknoten einen Spannungswert — dieselbe
    Knoten-ID kann in verschiedenen Elementen unterschiedliche s1-Werte haben
    (Sprung an Element-Grenzen). max(s1) ueber alle PRESOL-Eintraege ist der
    RAW-konsistente smax_nodal.

    Parameters
    ----------
    csv_path : str
        Pfad zur bereinigten -S-clean.csv (oder -A-S-clean.csv).

    Returns
    -------
    float | None
        Maximum von s1 ueber alle Zeilen, oder None bei Fehler/leerer Datei.
    """
    max_s1 = None
    try:
        with open(csv_path, 'r', encoding='utf-8', errors='ignore') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    s1 = float(row['s1'])
                except (ValueError, KeyError):
                    continue
                if max_s1 is None or s1 > max_s1:
                    max_s1 = s1
    except FileNotFoundError:
        return None
    except Exception:
        return None
    return max_s1


# ---------------------------------------------------------------------------
# Elementmittelwerte
# ---------------------------------------------------------------------------

def compute_element_means(data):
    """
    Berechnet elementweise Mittelwerte der Hauptspannungen aus Knotendaten.

    Parameters
    ----------
    data : dict
        {elem_id: [(node_id, s1, s2, s3), ...]}

    Returns
    -------
    dict
        {elem_id: (s1_mean, s2_mean, s3_mean)}
        Nur Elemente mit mindestens einem Knoteneintrag.
    """
    result = {}
    for elem_id, nodes in data.items():
        if not nodes:
            continue
        n = len(nodes)
        result[elem_id] = (
            sum(t[1] for t in nodes) / n,   # s1_mean
            sum(t[2] for t in nodes) / n,   # s2_mean
            sum(t[3] for t in nodes) / n,   # s3_mean
        )

    if result:
        print(f"[parse_presol] {len(result)} Elemente, Mittelwerte berechnet.")
    else:
        print("[parse_presol] WARNUNG: Keine Elementmittelwerte berechnet (leere Eingabe).")
    return result


# ---------------------------------------------------------------------------
# Gauss-RAW Hilfsfunktionen (V9)
# ---------------------------------------------------------------------------

def read_element_file_with_ids(filepath):
    """
    Liest effVol_Elemente.out im RAW-Format (9 Spalten fuer Hex).

    Spalte 1 = globale ANSYS Element-ID
    Spalten 2-9 = 8 Eckknoten-IDs (ANSYS-Konnektivitaetsreihenfolge)

    Parameters
    ----------
    filepath : str
        Pfad zur Elemente-Datei (z.B. tables/VEFF_{name}_Elemente.out).

    Returns
    -------
    tuple (elem_order, node_order)
        elem_order : list[int]
            Globale Element-IDs in sequentieller Reihenfolge (Zeile 1 -> Index 0).
        node_order : dict[int, list[int]]
            {glob_elem_id: [n1, n2, ..., n8]} — Eckknoten pro Element.

    Raises
    ------
    FileNotFoundError
        Wenn die Datei nicht existiert.
    ValueError
        Wenn eine Zeile nicht das erwartete Format hat.
    """
    elem_order = []
    node_order = {}

    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line_no, line in enumerate(f, 1):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            try:
                values = [int(float(p)) for p in parts]
            except ValueError as e:
                raise ValueError(
                    f"[read_element_file_with_ids] Zeile {line_no}: "
                    f"Konnte nicht parsen: '{stripped}' — {e}"
                )

            elem_id = values[0]
            node_ids = values[1:]
            elem_order.append(elem_id)
            node_order[elem_id] = node_ids

    n_elem = len(elem_order)
    if n_elem == 0:
        raise ValueError(
            f"[read_element_file_with_ids] Keine Elemente gefunden in {filepath}"
        )

    eckn = len(node_order[elem_order[0]])
    print(f"[parse_presol] Elemente-Datei gelesen: {n_elem} Elemente, "
          f"{eckn} Eckknoten/Element aus {os.path.basename(filepath)}")

    return elem_order, node_order


def write_fortran_stress_file(presol_data, elem_order, node_order, outpath):
    """
    Schreibt effVol_RawStress.out — elementlokale Spannungen, Fortran-freundlich.

    Fuer jedes Element (in elem_order-Reihenfolge) und fuer jeden seiner
    Eckknoten (in Konnektivitaetsreihenfolge aus node_order) werden
    die Hauptspannungen s1, s2, s3 aus presol_data nachgeschlagen und
    zeilenweise geschrieben.

    Format: nelem * eckn Zeilen, je 3 Spalten (s1, s2, s3), ES16.8-kompatibel.

    Parameters
    ----------
    presol_data : dict
        {elem_id: [(node_id, s1, s2, s3), ...]} — aus parse_presol() oder
        read_stress_csv().
    elem_order : list[int]
        Globale Element-IDs in sequentieller Reihenfolge (aus
        read_element_file_with_ids()).
    node_order : dict[int, list[int]]
        {glob_elem_id: [n1, n2, ..., n8]} — Eckknoten pro Element.
    outpath : str
        Zielpfad fuer die Fortran-Stress-Datei.

    Raises
    ------
    ValueError
        Wenn ein Element oder Knoten in presol_data fehlt.
    """
    n_written = 0
    n_missing_warn = 0

    with open(outpath, 'w', encoding='utf-8', newline='') as f:
        for seq_idx, elem_id in enumerate(elem_order):
            # PRESOL-Daten fuer dieses Element holen
            if elem_id not in presol_data:
                raise ValueError(
                    f"[write_fortran_stress_file] Element {elem_id} "
                    f"(seq. Index {seq_idx + 1}) fehlt in PRESOL-Daten."
                )

            # PRESOL-Knoten als dict fuer schnellen Lookup
            presol_nodes = {
                node_id: (s1, s2, s3)
                for node_id, s1, s2, s3 in presol_data[elem_id]
            }

            # Eckknoten in Konnektivitaetsreihenfolge durchgehen
            corner_nodes = node_order[elem_id]
            for node_id in corner_nodes:
                if node_id in presol_nodes:
                    s1, s2, s3 = presol_nodes[node_id]
                else:
                    # Eckknoten nicht in PRESOL gefunden — kritischer Fehler
                    raise ValueError(
                        f"[write_fortran_stress_file] Knoten {node_id} "
                        f"(Element {elem_id}, seq. Index {seq_idx + 1}) "
                        f"fehlt in PRESOL-Daten. "
                        f"Verfuegbare Knoten: {sorted(presol_nodes.keys())}"
                    )

                f.write(f"  {s1:16.8E}  {s2:16.8E}  {s3:16.8E}\n")
                n_written += 1

    print(f"[parse_presol] Fortran-Stress-Datei geschrieben: "
          f"{os.path.basename(outpath)} ({n_written} Zeilen, "
          f"{len(elem_order)} Elemente)")


# ---------------------------------------------------------------------------
# CLI-Entrypoint
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Verwendung: python parse_presol_output.py <pfad/zu/{name}-S.out>")
        print("Ausgabe:    <pfad/zu/{name}-S-clean.csv>  (gleicher Ordner)")
        sys.exit(1)

    in_file = sys.argv[1]
    if not os.path.exists(in_file):
        print(f"FEHLER: Datei nicht gefunden: {in_file}")
        sys.exit(1)

    # Ausgabepfad: gleicher Ordner, Suffix -S-clean.csv
    base = os.path.basename(in_file)
    if base.endswith('-S.out'):
        out_base = base.replace('-S.out', '-S-clean.csv')
    else:
        out_base = os.path.splitext(base)[0] + '-S-clean.csv'

    out_file = os.path.join(os.path.dirname(in_file), out_base)

    print(f"Eingabe:  {in_file}")
    print(f"Ausgabe:  {out_file}")

    parsed = parse_presol_to_csv(in_file, out_file)

    if parsed:
        means = compute_element_means(parsed)
        print(f"Fertig: {len(parsed)} Elemente, {len(means)} Elementmittelwerte.")
    else:
        print("Kein Ergebnis — Datei konnte nicht geparst werden.")
        sys.exit(1)
