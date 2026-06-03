# vtk_exporter.py
# VTK Legacy ASCII Export fuer Risk of Rupture Intensity (RRI) Felder.
#
# Erzeugt VTK Unstructured Grid Dateien (Cell Type 12 = VTK_HEXAHEDRON)
# mit CELL_DATA (RRI_norm, RRI pro m-Wert) und POINT_DATA (FEM-Spannungen,
# optionale analytische Spannungen).
#
# Pure Python — keine externen Dependencies (kein VTK, kein numpy).
#
# Datenquellen:
#   EM-Pfad:    tables/{apdl}_VTK_Nodes.out + tables/{apdl}_VTK_Elemente.out
#               (erzeugt von x_Export_VTK_Geometry.mac)
#               RRI-Stresses aus analyst.elements (korrekt gemittelt per AVG)
#   Gauss-Pfad: tables/VEFF_{apdl}_Nodes.out + tables/VEFF_{apdl}_Elemente.out
#               (identische Dateien wie Fortran nutzt → garantiert konsistent)

import json
import os
import math


# Konfigurierbar: Welche m-Werte exportieren
#M_VTK_DEFAULT = [5, 10, 15, 20, 25, 30, 40, 50]
M_VTK_DEFAULT = [1,2,3,4,5,6,7,8,9,10,
                 11,12,13,14,15,16,17,18,19,20,
                 21,22,23,24,25,26,27,28,29,30,
                 31,32,33,34,35,36,37,38,39,40,
                 41,42,43,44,45,46,47,48,49,50]  # Feinere Abstufung fuer kleine m


# =========================================================
# Geometrie-I/O
# =========================================================

def _read_sparse_nodes(filepath, ncols=6):
    """Liest ANSYS Sparse-Array (Zeilenindex = Knoten-ID).

    ANSYS schreibt Knoten als sparse Array: Zeile N = Knoten-ID N.
    Nicht-existierende Knoten sind Null-Zeilen.

    Args:
        filepath: Pfad zur *.out Datei
        ncols: Erwartete Spaltenanzahl (6 = X,Y,Z,S1,S2,S3 oder 3 = X,Y,Z)

    Returns:
        dict {ansys_node_id: tuple(vals)} — nur Non-Zero-Zeilen
    """
    nodes = {}
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            parts = line.split()
            if len(parts) < ncols:
                continue
            vals = [float(p) for p in parts[:ncols]]
            # Null-Zeile ueberspringen (nicht-selektierter Knoten)
            if all(abs(v) < 1e-30 for v in vals):
                continue
            nodes[line_num] = tuple(vals)
    return nodes


def _read_elements(filepath):
    """Liest Elemente.out (kompaktes Format).

    Erkennt automatisch ob 8 Spalten (NOD: nur Knoten) oder
    9 Spalten (RAW/VTK: ElemID + 8 Knoten).

    Returns:
        list of tuples: [(elem_id_or_none, n1, ..., n8), ...]
        Bei 8 Spalten: elem_id = None (sequentieller Index)
        Bei 9 Spalten: elem_id = int(Spalte 1)
    """
    elements = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            vals = [int(round(float(p))) for p in parts]
            if len(vals) == 9:
                # RAW/VTK-Format: ElemID + 8 Knoten
                elements.append(tuple(vals))
            elif len(vals) == 8:
                # NOD-Format: Nur 8 Knoten (kein ElemID)
                elements.append((None,) + tuple(vals))
            elif len(vals) >= 9:
                # Mehr als 9 Spalten: Nimm erste 9
                elements.append(tuple(vals[:9]))
    return elements


def _build_vtk_mesh(node_dict, elem_list):
    """Baut VTK-kompatibles Mesh aus ANSYS-Daten.

    ANSYS nutzt globale Knoten-IDs (1-basiert, nicht-dicht).
    VTK braucht 0-basierte, dichte Indices.

    Returns:
        points: list of (x, y, z) — sortiert nach ANSYS-ID
        cells: list of [vtk_n1, ..., vtk_n8]
        ansys_to_vtk: dict {ansys_id: vtk_index}
        sorted_node_ids: list of ANSYS-IDs (Reihenfolge = VTK-Index)
    """
    sorted_ids = sorted(node_dict.keys())
    ansys_to_vtk = {aid: idx for idx, aid in enumerate(sorted_ids)}
    points = [node_dict[aid][:3] for aid in sorted_ids]

    cells = []
    for elem in elem_list:
        # Letzte 8 Werte = Eckknoten (bei 9 Spalten: [1:9], bei 8+None: [1:9])
        node_ids = [int(n) for n in elem[1:9]]
        try:
            vtk_cell = [ansys_to_vtk[nid] for nid in node_ids]
            cells.append(vtk_cell)
        except KeyError as e:
            # Knoten fehlt im sparse Array — sollte nicht passieren
            print(f"   [VTK] WARNUNG: Knoten {e} nicht in Node-Dictionary. "
                  f"Element uebersprungen.")
            cells.append(None)

    return points, cells, ansys_to_vtk, sorted_ids


# =========================================================
# Gauss-Pfad: Element-Spannungen aus Knoten berechnen
# =========================================================

def _compute_gauss_element_means(node_dict, elem_list):
    """Element-Mittelwerte aus nodalen Spannungen (8 Eckknoten).

    Fuer Gauss-NOD-Pfad: node_dict hat 6 Spalten (X,Y,Z,S1,S2,S3).

    Returns:
        list of (s1_mean, s2_mean, s3_mean) pro Element
    """
    means = []
    for elem in elem_list:
        node_ids = [int(n) for n in elem[1:9]]
        s1_sum = s2_sum = s3_sum = 0.0
        count = 0
        for nid in node_ids:
            if nid in node_dict and len(node_dict[nid]) >= 6:
                s1_sum += node_dict[nid][3]
                s2_sum += node_dict[nid][4]
                s3_sum += node_dict[nid][5]
                count += 1
        if count > 0:
            means.append((s1_sum / count, s2_sum / count, s3_sum / count))
        else:
            means.append((0.0, 0.0, 0.0))
    return means


# =========================================================
# RRI-Berechnung (EXAKT wie analyze_em() Zeilen 618-637)
# =========================================================

def _compute_rri(s1, s2, s3, smax_norm, sigma0, crit, m):
    """Berechnet RRI_norm und RRI fuer ein Element bei gegebenem m.

    Konsistenz-Garantie: Identische Logik wie ResultAnalyst.analyze_em().

    Args:
        s1, s2, s3: Element-Spannungen (Mittelwerte)
        smax_norm: Normierungswert
        sigma0: Weibull-Skalierungsparameter
        crit: "PIA" oder "S1"
        m: Weibull-Modul

    Returns:
        (rri_norm, rri): Normiertes und physikalisches RRI
    """
    if smax_norm <= 0 or sigma0 <= 0:
        return 0.0, 0.0

    if crit == "S1":
        # Nur erste Hauptspannung
        if s1 > 0:
            ratio = s1 / smax_norm
            rri_norm = ratio ** m
        else:
            rri_norm = 0.0
    else:
        # PIA: Alle positiven Hauptspannungen
        ratios = []
        for s in (s1, s2, s3):
            if s > 0:
                ratios.append(s / smax_norm)
        rri_norm = sum(r ** m for r in ratios)

    # RRI = RRI_norm * (smax_norm / sigma0)^m
    load_factor = (smax_norm / sigma0) ** m
    rri = rri_norm * load_factor

    return rri_norm, rri


# =========================================================
# VTK-Schreiber
# =========================================================

def _write_vtk_file(filepath, points, cells, cell_data, point_data, title):
    """Schreibt VTK Legacy ASCII Unstructured Grid.

    Cell-Type 12 = VTK_HEXAHEDRON (8-Knoten Hex).
    Direkt kompatibel mit ANSYS SOLID186 Corner-Nodes.

    Args:
        filepath: Ausgabepfad (.vtk)
        points: list of (x, y, z)
        cells: list of [n1, ..., n8] (0-basierte VTK-Indices)
        cell_data: dict {name: [values...]} — CELL_DATA Felder
        point_data: dict {name: [values...]} — POINT_DATA Felder
        title: Titel-String
    """
    # Zellen ohne None filtern
    valid_cells = [c for c in cells if c is not None]
    valid_cell_indices = [i for i, c in enumerate(cells) if c is not None]

    npts = len(points)
    ncells = len(valid_cells)

    with open(filepath, 'w') as f:
        # Header
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"{title}\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write("\n")

        # Points
        f.write(f"POINTS {npts} double\n")
        for x, y, z in points:
            f.write(f"{x:.8e} {y:.8e} {z:.8e}\n")
        f.write("\n")

        # Cells (8 Knoten pro Hex)
        total_ints = ncells * 9  # 1 (Knotenanzahl) + 8 (Knoten-IDs) pro Zelle
        f.write(f"CELLS {ncells} {total_ints}\n")
        for cell in valid_cells:
            f.write(f"8 {' '.join(str(n) for n in cell)}\n")
        f.write("\n")

        # Cell Types (alle VTK_HEXAHEDRON = 12)
        f.write(f"CELL_TYPES {ncells}\n")
        for _ in range(ncells):
            f.write("12\n")
        f.write("\n")

        # Cell Data
        if cell_data:
            f.write(f"CELL_DATA {ncells}\n")
            for name, values in cell_data.items():
                f.write(f"SCALARS {name} double 1\n")
                f.write("LOOKUP_TABLE default\n")
                for i in valid_cell_indices:
                    f.write(f"{values[i]:.8e}\n")
                f.write("\n")

        # Point Data
        if point_data:
            f.write(f"POINT_DATA {npts}\n")
            for name, values in point_data.items():
                f.write(f"SCALARS {name} double 1\n")
                f.write("LOOKUP_TABLE default\n")
                for v in values:
                    f.write(f"{v:.8e}\n")
                f.write("\n")


# =========================================================
# Series-Datei fuer ParaView
# =========================================================

def _write_series_file(vtk_dir, safe_name, m_list):
    """Schreibt .vtk.series JSON-Datei fuer ParaView Zeitserie.

    Ausgabe: {vtk_dir}/series/{safe_name}.vtk.series
    Eintraege: relative Pfade ../{safe_name}_RRI_m{m:02d}.vtk mit time=m.
    ParaView erkennt die .series-Datei als Zeitserie (kein manuelles Gruppieren).

    Returns:
        str: Pfad zur erzeugten Series-Datei
    """
    series_dir = os.path.join(vtk_dir, "series")
    os.makedirs(series_dir, exist_ok=True)

    entries = [
        {"name": f"../{safe_name}_RRI_m{m:02d}.vtk", "time": float(m)}
        for m in sorted(m_list)
    ]
    content = {"file-series-version": "1.0", "files": entries}

    series_path = os.path.join(series_dir, f"{safe_name}.vtk.series")
    with open(series_path, 'w', encoding='utf-8') as f:
        json.dump(content, f, indent=2)
    return series_path


# =========================================================
# Oeffentliche API
# =========================================================

def export_vtk_series(
    tables_dir,
    vtk_dir,
    case_name,
    apdl_name,
    int_type,
    avg,
    crit,
    smax_norm,
    sigma0,
    symmetry_factor,
    elements_em=None,
    element_ids_em=None,
    m_values=None,
    case_type="UNKNOWN",
    case_x=0, case_y=0, case_z=0,
    load_n=0,
    export_analytical=True,
):
    """Exportiert VTK-Dateiserie (eine pro m-Wert) mit RRI-Feldern.

    Args:
        tables_dir: Pfad zu tables/ Ordner
        vtk_dir: Ausgabe-Ordner (wird erstellt)
        case_name: full_name (fuer Dateinamen)
        apdl_name: APDL-Name (fuer File-Lookups)
        int_type: "EM" oder "G5" etc.
        avg: "RAW" oder "NOD"
        crit: "PIA" oder "S1"
        smax_norm: Normierungswert
        sigma0: Weibull-Skalierungsparameter
        symmetry_factor: 1.0 oder 4.0
        elements_em: [(vol, s1, s2, s3), ...] fuer EM (aus analyst.elements)
        element_ids_em: [id1, id2, ...] parallel zu elements_em
        m_values: Custom m-Werte (Default: M_VTK_DEFAULT)
        case_type: Lastfall-Typ fuer analytische Spannungen
        case_x, case_y, case_z: Geometrie-Parameter
        load_n: Last [N oder MPa]
        export_analytical: True = analytische Spannungen als POINT_DATA

    Returns:
        bool: True bei Erfolg
    """
    os.makedirs(vtk_dir, exist_ok=True)
    m_list = m_values or M_VTK_DEFAULT

    # --- Geometrie laden ---
    node_dict = None
    elem_list = None
    elem_stresses = None

    if int_type == "EM":
        # EM: Geometrie aus VTK-Export-Files
        nodes_file = os.path.join(tables_dir, f"{apdl_name}_VTK_Nodes.out")
        elems_file = os.path.join(tables_dir, f"{apdl_name}_VTK_Elemente.out")

        if not os.path.exists(nodes_file) or not os.path.exists(elems_file):
            print(f"   [VTK] WARNUNG: Geometrie-Dateien fuer EM nicht gefunden: "
                  f"{nodes_file}")
            return False

        node_dict = _read_sparse_nodes(nodes_file, ncols=6)
        elem_list = _read_elements(elems_file)

        # RRI-Stresses aus analyst.elements via element_ids Mapping
        if elements_em and element_ids_em:
            elem_stress_dict = {}
            for eid, (vol, s1, s2, s3) in zip(element_ids_em, elements_em):
                elem_stress_dict[eid] = (s1, s2, s3)
            # Reihenfolge: elem_list definiert die Reihenfolge fuer CELL_DATA
            elem_stresses = [elem_stress_dict.get(int(e[0]), (0.0, 0.0, 0.0))
                             for e in elem_list]
        else:
            print("   [VTK] WARNUNG: Keine EM-Elementdaten fuer RRI.")
            return False

    elif int_type.startswith("G"):
        # Gauss: Geometrie aus VEFF-Dateien (identisch mit Fortran)
        elems_file = os.path.join(tables_dir, f"VEFF_{apdl_name}_Elemente.out")
        if not os.path.exists(elems_file):
            print(f"   [VTK] WARNUNG: Element-Datei nicht gefunden: {elems_file}")
            return False

        elem_list = _read_elements(elems_file)

        # Koordinaten: Entweder NodeCoords (RAW) oder Nodes (NOD, 6 Spalten)
        coords_file = os.path.join(tables_dir, f"VEFF_{apdl_name}_NodeCoords.out")
        nodes_file = os.path.join(tables_dir, f"VEFF_{apdl_name}_Nodes.out")

        if os.path.exists(coords_file):
            # RAW-Pfad: Nur Koordinaten (3 Spalten)
            node_dict = _read_sparse_nodes(coords_file, ncols=3)
        elif os.path.exists(nodes_file):
            # NOD-Pfad: Koordinaten + Spannungen (6 Spalten)
            node_dict = _read_sparse_nodes(nodes_file, ncols=6)
        else:
            print(f"   [VTK] WARNUNG: Keine Knoten-Datei gefunden fuer Gauss.")
            return False

        # Gauss-Stresses: Element-Mittelwerte aus nodalen Daten
        if len(node_dict) > 0 and len(list(node_dict.values())[0]) >= 6:
            elem_stresses = _compute_gauss_element_means(node_dict, elem_list)
        else:
            # RAW-Pfad: Keine nodalen Spannungen -> Nullen (RRI nur mit analyst)
            elem_stresses = [(0.0, 0.0, 0.0)] * len(elem_list)
            print("   [VTK] HINWEIS: Gauss-RAW — nodale Spannungen nicht verfuegbar, "
                  "RRI aus Element-Mittelwerten nicht moeglich.")

    else:
        print(f"   [VTK] Unbekannter Integrationstyp: {int_type}")
        return False

    if not node_dict or not elem_list:
        print("   [VTK] WARNUNG: Leere Geometrie.")
        return False

    # --- VTK-Mesh aufbauen ---
    points, cells, a2v_map, sorted_nids = _build_vtk_mesh(node_dict, elem_list)

    print(f"   [VTK] Mesh: {len(points)} Knoten, {len(cells)} Elemente")

    # --- Nodale Spannungen fuer POINT_DATA (nur Visualisierung) ---
    nodal_s1, nodal_s2, nodal_s3 = [], [], []
    for nid in sorted_nids:
        nd = node_dict[nid]
        if len(nd) >= 6:
            nodal_s1.append(nd[3])
            nodal_s2.append(nd[4])
            nodal_s3.append(nd[5])
        else:
            nodal_s1.append(0.0)
            nodal_s2.append(0.0)
            nodal_s3.append(0.0)

    # --- Analytische Spannungen (optional) ---
    analytical = None
    if export_analytical:
        try:
            from analytical_stress_fields import compute_analytical_stresses
            analytical = compute_analytical_stresses(
                case_type, points, case_x, case_y, case_z, load_n)
        except ImportError:
            print("   [VTK] HINWEIS: analytical_stress_fields nicht verfuegbar.")
        except Exception as e:
            print(f"   [VTK] HINWEIS: Analytik fehlgeschlagen: {e}")

    # --- Pro m-Wert: RRI berechnen + VTK schreiben ---
    safe_name = case_name.replace(".", "_")
    for m in m_list:
        rri_norm_list = []
        rri_list = []

        for s1, s2, s3 in elem_stresses:
            rn, r = _compute_rri(s1, s2, s3, smax_norm, sigma0, crit, m)
            rri_norm_list.append(rn)
            rri_list.append(r)

        # RRI_norm = RRI / RRI_max → max = 1, σ₀ kürzt sich (σ-Feld-Normierung).
        # RRI bleibt physikalisch: (σ/σ₀)^m, σ₀-abhängig, nicht normiert.
        max_rri = max(rri_norm_list) if rri_norm_list else 0.0
        scale = max_rri if max_rri > 0 else 1.0
        rri_norm_list = [v / scale for v in rri_norm_list]
        # rri_list unverändert: physikalischer RRI-Wert (σ/σ₀)^m

        # CELL_DATA
        cell_data = {
            "RRI_norm": rri_norm_list,
            "RRI": rri_list,
        }

        # POINT_DATA
        point_data = {
            "FEM_S1": nodal_s1,
            "FEM_S2": nodal_s2,
            "FEM_S3": nodal_s3,
        }
        if analytical:
            point_data["Analytical_S1"] = analytical[0]
            point_data["Analytical_S2"] = analytical[1]
            point_data["Analytical_S3"] = analytical[2]

        vtk_path = os.path.join(vtk_dir, f"{safe_name}_RRI_m{m:02d}.vtk")
        _write_vtk_file(
            vtk_path, points, cells, cell_data, point_data,
            title=f"RRI Field: {case_name} m={m}")

    series_path = _write_series_file(vtk_dir, safe_name, m_list)
    print(f"   [VTK] {len(m_list)} VTK-Dateien geschrieben nach {vtk_dir}")
    print(f"   [VTK] Series: {os.path.relpath(series_path, vtk_dir)}")
    return True


# =========================================================
# VTK-Surface-Export (Quad-Cells, Zwei-Datei-Architektur)
# =========================================================

def _write_quad_cells_vtk(filepath, points, quad_cells, cell_data, point_data, title):
    """VTK Legacy ASCII Writer fuer Quad-Cells (CELL_TYPE = 9, VTK_QUAD).

    Analog zu _write_vtk_file, aber mit 4-Knoten-Quads statt 8-Knoten-Hex.

    Args:
        filepath: Ausgabepfad (.vtk)
        points: list of (x, y, z)
        quad_cells: list of [n1, n2, n3, n4] (0-basierte VTK-Indices)
        cell_data: dict {name: [values per cell]}
        point_data: dict {name: [values per point]}
        title: Titel-String
    """
    valid_cells = [c for c in quad_cells if c is not None and len(c) == 4]
    valid_cell_indices = [i for i, c in enumerate(quad_cells) if c is not None and len(c) == 4]

    npts = len(points)
    ncells = len(valid_cells)

    with open(filepath, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"{title}\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write("\n")

        f.write(f"POINTS {npts} double\n")
        for x, y, z in points:
            f.write(f"{x:.8e} {y:.8e} {z:.8e}\n")
        f.write("\n")

        # Cells (4 Knoten pro Quad)
        total_ints = ncells * 5  # 1 (Knotenanzahl=4) + 4 Knoten-IDs
        f.write(f"CELLS {ncells} {total_ints}\n")
        for cell in valid_cells:
            f.write(f"4 {' '.join(str(n) for n in cell)}\n")
        f.write("\n")

        # Cell Types (alle VTK_QUAD = 9)
        f.write(f"CELL_TYPES {ncells}\n")
        for _ in range(ncells):
            f.write("9\n")
        f.write("\n")

        if cell_data:
            f.write(f"CELL_DATA {ncells}\n")
            for name, values in cell_data.items():
                f.write(f"SCALARS {name} double 1\n")
                f.write("LOOKUP_TABLE default\n")
                for i in valid_cell_indices:
                    f.write(f"{values[i]:.8e}\n")
                f.write("\n")

        if point_data:
            f.write(f"POINT_DATA {npts}\n")
            for name, values in point_data.items():
                f.write(f"SCALARS {name} double 1\n")
                f.write("LOOKUP_TABLE default\n")
                for v in values:
                    f.write(f"{v:.8e}\n")
                f.write("\n")


def _compute_surface_rri(s1, s2, s3, smax_norm, sigma0, crit, m):
    """RRI pro Surface-Face. Identische Logik wie _compute_rri (Volumen).

    Returns:
        (rri_norm, rri) Tupel
    """
    return _compute_rri(s1, s2, s3, smax_norm, sigma0, crit, m)


def _write_surface_series_file(vtk_dir, safe_name, m_list):
    """Series-Datei fuer Surface-RRI-Zeitserie."""
    series_dir = os.path.join(vtk_dir, "series")
    os.makedirs(series_dir, exist_ok=True)

    entries = [
        {"name": f"../{safe_name}_surface_RRI_m{m:02d}.vtk", "time": float(m)}
        for m in sorted(m_list)
    ]
    content = {"file-series-version": "1.0", "files": entries}

    series_path = os.path.join(series_dir, f"{safe_name}_surface.vtk.series")
    with open(series_path, 'w', encoding='utf-8') as f:
        json.dump(content, f, indent=2)
    return series_path


def export_vtk_volume_context(tables_dir, vtk_dir, case_name, apdl_name, int_type):
    """v12.1: Exportiert Volumenmesh ohne Cell-Data als ParaView-Kontext (graue Geometrie).

    Genutzt im LOC=A-Modus, damit ParaView die Volumenkontur transluzent darstellen kann
    waehrend die Surface-Faces farbig nach RRI eingefaerbt sind.

    Args:
        tables_dir: Pfad zu tables/ Ordner
        vtk_dir: Ausgabe-Ordner
        case_name: full_name (fuer Dateinamen)
        apdl_name: APDL-Name (fuer File-Lookups)
        int_type: "EM" oder "G5" — bestimmt Geometrie-Quelle

    Returns:
        bool: True bei Erfolg
    """
    os.makedirs(vtk_dir, exist_ok=True)

    if int_type == "EM":
        nodes_file = os.path.join(tables_dir, f"{apdl_name}_VTK_Nodes.out")
        elems_file = os.path.join(tables_dir, f"{apdl_name}_VTK_Elemente.out")
    else:
        elems_file = os.path.join(tables_dir, f"VEFF_{apdl_name}_Elemente.out")
        coords_file = os.path.join(tables_dir, f"VEFF_{apdl_name}_NodeCoords.out")
        nodes_file = coords_file if os.path.exists(coords_file) else (
            os.path.join(tables_dir, f"VEFF_{apdl_name}_Nodes.out"))

    if not os.path.exists(nodes_file) or not os.path.exists(elems_file):
        print(f"   [VTK-Surface] WARNUNG: Volumen-Geometrie nicht gefunden fuer Kontext-Export.")
        return False

    ncols = 3 if "NodeCoords" in nodes_file else 6
    node_dict = _read_sparse_nodes(nodes_file, ncols=ncols)
    elem_list = _read_elements(elems_file)

    points, cells, _, _ = _build_vtk_mesh(node_dict, elem_list)

    safe_name = case_name.replace(".", "_")
    vtk_path = os.path.join(vtk_dir, f"{safe_name}_volume_context.vtk")
    _write_vtk_file(
        vtk_path, points, cells,
        cell_data=None, point_data=None,
        title=f"Volume Context (gray): {case_name}")
    print(f"   [VTK-Surface] Volume-Context: {os.path.basename(vtk_path)}")
    return True


def export_vtk_surface_series(tables_dir, vtk_dir, case_name, apdl_name, int_type,
                              surface_faces, smax_norm, sigma0, crit,
                              m_values=None):
    """v12.1: Exportiert Surface-RRI als Quad-Cell-VTK pro m-Wert + Series-File.

    Args:
        tables_dir: Pfad zu tables/ Ordner
        vtk_dir: Ausgabe-Ordner
        case_name: full_name
        apdl_name: APDL-Name (fuer Knoten-Koordinaten-Lookup)
        int_type: "EM" oder "G5"
        surface_faces: list of (area, s1, s2, s3, [n1, n2, n3, n4])
            (aus analyst.surface_faces nach _load_surface_properties)
        smax_norm: Normierungs-Spannung (Pa oder MPa)
        sigma0: Weibull-Skalierungsparameter
        crit: "PIA" oder "S1"
        m_values: Custom m-Werte (Default: M_VTK_DEFAULT)

    Returns:
        bool: True bei Erfolg
    """
    os.makedirs(vtk_dir, exist_ok=True)
    m_list = m_values or M_VTK_DEFAULT

    # Knoten-Koordinaten laden (fuer Surface-Knoten)
    if int_type == "EM":
        nodes_file = os.path.join(tables_dir, f"{apdl_name}_VTK_Nodes.out")
    else:
        coords_file = os.path.join(tables_dir, f"VEFF_{apdl_name}_NodeCoords.out")
        nodes_file = coords_file if os.path.exists(coords_file) else (
            os.path.join(tables_dir, f"VEFF_{apdl_name}_Nodes.out"))

    if not os.path.exists(nodes_file):
        print(f"   [VTK-Surface] WARNUNG: Knoten-Datei nicht gefunden: {nodes_file}")
        return False

    ncols = 3 if "NodeCoords" in nodes_file else 6
    node_dict = _read_sparse_nodes(nodes_file, ncols=ncols)

    # Filter: nur Faces mit Knoten-IDs (sollte ab v12.1 immer der Fall sein)
    valid_faces = [f for f in surface_faces if len(f) >= 5 and f[4] is not None]
    if not valid_faces:
        print(f"   [VTK-Surface] WARNUNG: Keine Surface-Faces mit Knoten-IDs.")
        return False

    # Eindeutige Knoten-IDs sammeln + auf 0-basierte VTK-Indices mappen
    all_nodes = set()
    for _area, _s1, _s2, _s3, n_list in valid_faces:
        for nid in n_list:
            if nid in node_dict:
                all_nodes.add(nid)

    sorted_nids = sorted(all_nodes)
    a2v = {nid: i for i, nid in enumerate(sorted_nids)}

    # Punkte-Liste (nur die Surface-Knoten)
    points = []
    for nid in sorted_nids:
        nd = node_dict[nid]
        points.append((nd[0], nd[1], nd[2]))

    # Quad-Cells aus Face-Knoten-IDs
    quad_cells = []
    for _area, _s1, _s2, _s3, n_list in valid_faces:
        try:
            quad = [a2v[nid] for nid in n_list]
            quad_cells.append(quad)
        except KeyError:
            quad_cells.append(None)  # wird vom Writer gefiltert

    # Nodale Spannungen fuer POINT_DATA (aus -A.out Face-Mittelwerten interpoliert)
    # Pro Surface-Knoten: Mittel aller anliegenden Face-S1/S2/S3
    node_stress_sum = {nid: [0.0, 0.0, 0.0, 0] for nid in sorted_nids}  # [s1_sum, s2_sum, s3_sum, count]
    for _area, s1, s2, s3, n_list in valid_faces:
        for nid in n_list:
            if nid in node_stress_sum:
                node_stress_sum[nid][0] += s1
                node_stress_sum[nid][1] += s2
                node_stress_sum[nid][2] += s3
                node_stress_sum[nid][3] += 1

    nodal_s1 = []
    nodal_s2 = []
    nodal_s3 = []
    for nid in sorted_nids:
        s1_sum, s2_sum, s3_sum, count = node_stress_sum[nid]
        if count > 0:
            nodal_s1.append(s1_sum / count)
            nodal_s2.append(s2_sum / count)
            nodal_s3.append(s3_sum / count)
        else:
            nodal_s1.append(0.0)
            nodal_s2.append(0.0)
            nodal_s3.append(0.0)

    safe_name = case_name.replace(".", "_")

    # Pro m-Wert RRI berechnen + VTK schreiben
    for m in m_list:
        rri_norm_list = []
        rri_list = []
        for area, s1, s2, s3, _n_list in valid_faces:
            rn, r = _compute_surface_rri(s1, s2, s3, smax_norm, sigma0, crit, m)
            rri_norm_list.append(rn)
            rri_list.append(r)

        # RRI_norm = RRI / RRI_max → max = 1 (analog Volumen-Export)
        # RRI bleibt physikalisch (σ/σ₀)^m, nicht normiert
        max_rri = max(rri_norm_list) if rri_norm_list else 0.0
        scale = max_rri if max_rri > 0 else 1.0
        rri_norm_list = [v / scale for v in rri_norm_list]

        cell_data = {
            "RRI_norm": rri_norm_list,
            "RRI": rri_list,
            "Area": [f[0] for f in valid_faces],
        }
        point_data = {
            "FEM_S1": nodal_s1,
            "FEM_S2": nodal_s2,
            "FEM_S3": nodal_s3,
        }

        vtk_path = os.path.join(vtk_dir, f"{safe_name}_surface_RRI_m{m:02d}.vtk")
        _write_quad_cells_vtk(
            vtk_path, points, quad_cells, cell_data, point_data,
            title=f"Surface RRI: {case_name} m={m}")

    series_path = _write_surface_series_file(vtk_dir, safe_name, m_list)
    print(f"   [VTK-Surface] {len(m_list)} Surface-VTK-Dateien geschrieben (Quad-Cells)")
    print(f"   [VTK-Surface] Series: {os.path.relpath(series_path, vtk_dir)}")
    return True


# =========================================================
# Gauss-Punkt Placeholder (Zukunft)
# =========================================================

def export_vtk_gauss_points(tables_dir, vtk_dir, case_name, **kwargs):
    """ZUKUNFT: RRI an Gauss-Integrationspunkten als vtkPolyData.

    Erfordert Fortran-Export: GP-Koordinaten (X,Y,Z), V_sub = w*det(J), RRI_norm.
    ParaView: Glyph Filter, Wuerfel nach V_sub skaliert, nach RRI eingefaerbt.
    """
    raise NotImplementedError("Gauss-Punkt-VTK-Export noch nicht implementiert.")
