# vtk_refinement_exporter.py
# Mesh Refinement Indicator VTK Export.
#
# Berechnet pro Element einen Refinement-Indikator basierend auf der C3D-sum-
# Methodik aus worknotes/masterarbeit-fehlerindikator-herleitung.md (Variante D
# Proportional Slack-Reuse).
#
# Output (pro m-Wert eine VTK-Datei):
#   - 3 Skalare:  E (Element-Fehler), W (Hazard-Gewicht), I (Indikator)
#   - 2 Vektoren: r_e_local, h_ideal_local
#                 ACHTUNG: Komponenten beziehen sich auf LOKALE Element-Achsen,
#                 NICHT auf globale x/y/z. Bei strukturierten Meshes mit
#                 gekruemmter Topologie (z.B. Lochrand) sind die lokalen Achsen
#                 gegenueber globalen rotiert.
#   - 9 Cell-Data-Felder pro VTK-Datei
#
# Lokale Achsen-Konvention (HEX8):
#   loc1 = von Face-Mitte (Knoten 0,3,4,7) zu Face-Mitte (Knoten 1,2,5,6)
#   loc2 = von Face-Mitte (Knoten 0,1,4,5) zu Face-Mitte (Knoten 2,3,6,7)
#   loc3 = von Bottom-Face (0,1,2,3) zu Top-Face (4,5,6,7)
# Diese Achsen entsprechen der ANSYS-Mesher-Logik. Bei achsparallelen
# Volumen (z.B. POR-Innenwuerfel) sind sie identisch zu globalen x/y/z.
# Bei rotierten Elementen (z.B. PWH-Lochrand) divergieren sie.
#
# LOC=V: Hex-Cell (8 Eckknoten, 3 lokale Richtungen)
# LOC=A: Quad-Cell (4 Eckknoten, 2 lokale Richtungen + 3.-Slot=Dummy)
#
# Pure Python (keine Dependencies). Reuse von vtk_exporter.py-Helfern.

import csv
import math
import os
import json
import sys


# =========================================================
# Face-Index-Sets (0-basiert in Eckknoten-Tupel)
# Verifiziert via Element 1 in Run 2026-04-30_15-36 (POR-V).
# Standard-VTK-Hex8-Konvention (CGNS-konform):
#   Position 1-4: Bottom-Face (z-)
#   Position 5-8: Top-Face (z+)
#   CCW-Reihenfolge in jeder Z-Ebene:
#     Idx 0: (0,0)   Idx 1: (1,0)   Idx 2: (1,1)   Idx 3: (0,1)
#     Idx 4: (0,0)   Idx 5: (1,0)   Idx 6: (1,1)   Idx 7: (0,1)
# =========================================================

HEX8_FACE_INDICES = {
    "x_minus": (0, 3, 4, 7),  # Knoten an x=0
    "x_plus":  (1, 2, 5, 6),  # Knoten an x=h_x
    "y_minus": (0, 1, 4, 5),  # Knoten an y=0
    "y_plus":  (2, 3, 6, 7),  # Knoten an y=h_y
    "z_minus": (0, 1, 2, 3),  # Knoten an z=0
    "z_plus":  (4, 5, 6, 7),  # Knoten an z=h_z
}

# Quad4-Konvention (lokale u-v):
#   Idx 0: (0,0)  Idx 1: (1,0)  Idx 2: (1,1)  Idx 3: (0,1)
QUAD4_EDGE_INDICES = {
    "u_minus": (0, 3),
    "u_plus":  (1, 2),
    "v_minus": (0, 1),
    "v_plus":  (3, 2),
}


# =========================================================
# Math-Helpers
# =========================================================

def _macaulay(s):
    """Macaulay-Bracket: max(s, 0)."""
    return s if s > 0.0 else 0.0


def _face_mean_value(corner_stresses, face_indices):
    """Mittel der Spannungen auf einer Face (mit Macaulay)."""
    vals = [_macaulay(corner_stresses[i]) for i in face_indices]
    n = len(vals)
    return sum(vals) / n if n > 0 else 0.0


def _face_mean_position(corner_coords, face_indices):
    """3D-Mittelpunkt-Position der Face-Knoten."""
    n = len(face_indices)
    if n == 0:
        return (0.0, 0.0, 0.0)
    sx = sum(corner_coords[i][0] for i in face_indices) / n
    sy = sum(corner_coords[i][1] for i in face_indices) / n
    sz = sum(corner_coords[i][2] for i in face_indices) / n
    return (sx, sy, sz)


def _compute_corner_mean(corner_stresses):
    """sigma_e = Mittel aller Eckknoten-Spannungen (mit Macaulay)."""
    pos = [_macaulay(s) for s in corner_stresses]
    n = len(pos)
    return sum(pos) / n if n > 0 else 0.0


def _compute_delta_and_h_hex(corner_stresses, corner_coords):
    """Hex8: Delta_k sigma und h_k fuer jede der 3 Raumrichtungen.

    Returns:
        deltas:    (Delta_x, Delta_y, Delta_z)
        h_lengths: (h_x, h_y, h_z)
    """
    deltas = []
    hs = []
    for face_minus, face_plus in [
        ("x_minus", "x_plus"),
        ("y_minus", "y_plus"),
        ("z_minus", "z_plus"),
    ]:
        idx_minus = HEX8_FACE_INDICES[face_minus]
        idx_plus = HEX8_FACE_INDICES[face_plus]

        sigma_minus = _face_mean_value(corner_stresses, idx_minus)
        sigma_plus = _face_mean_value(corner_stresses, idx_plus)
        deltas.append(sigma_plus - sigma_minus)

        pos_minus = _face_mean_position(corner_coords, idx_minus)
        pos_plus = _face_mean_position(corner_coords, idx_plus)
        h = math.sqrt(sum((pos_plus[k] - pos_minus[k]) ** 2 for k in range(3)))
        hs.append(h)

    return tuple(deltas), tuple(hs)


def _compute_delta_and_h_quad(corner_stresses, corner_coords):
    """Quad4: Delta_k sigma und h_k fuer 2 lokale Richtungen u, v.

    z-Slot wird mit 0 gefuellt fuer konsistente VTK-Vector-Struktur.

    Returns:
        deltas:    (Delta_u, Delta_v, 0.0)
        h_lengths: (h_u, h_v, 0.0)
    """
    deltas = []
    hs = []
    for edge_minus, edge_plus in [
        ("u_minus", "u_plus"),
        ("v_minus", "v_plus"),
    ]:
        idx_minus = QUAD4_EDGE_INDICES[edge_minus]
        idx_plus = QUAD4_EDGE_INDICES[edge_plus]

        sigma_minus = _face_mean_value(corner_stresses, idx_minus)
        sigma_plus = _face_mean_value(corner_stresses, idx_plus)
        deltas.append(sigma_plus - sigma_minus)

        pos_minus = _face_mean_position(corner_coords, idx_minus)
        pos_plus = _face_mean_position(corner_coords, idx_plus)
        h = math.sqrt(sum((pos_plus[k] - pos_minus[k]) ** 2 for k in range(3)))
        hs.append(h)

    return tuple(deltas) + (0.0,), tuple(hs) + (0.0,)


def _compute_E_per_direction(sigma_e, deltas, m, sigma_eps=1e-12):
    """E_e^(k) = m(m-1) / (24 sigma_e^2) * (Delta_k sigma)^2 pro Richtung."""
    if sigma_e <= sigma_eps:
        return tuple(0.0 for _ in deltas)

    factor = m * (m - 1.0) / (24.0 * sigma_e * sigma_e)
    return tuple(factor * (d * d) for d in deltas)


def _compute_hazard_weights(sigma_e_list, volumes, m, sigma_eps=1e-12):
    """Hazard-Gewicht W_e pro Element (sum-normiert: Sigma W_e = 1).

    W_e = (sigma_e^m * V_e) / Sum_j (sigma_j^m * V_j)
    """
    n = len(sigma_e_list)
    raw_hazards = []
    for i in range(n):
        sigma = sigma_e_list[i]
        vol = volumes[i]
        if sigma <= sigma_eps or vol <= 0.0:
            raw_hazards.append(0.0)
        else:
            raw_hazards.append((sigma ** m) * vol)

    total = sum(raw_hazards)
    if total <= 0.0:
        return [0.0] * n

    return [h / total for h in raw_hazards]


# =========================================================
# Variante D — Proportional Slack-Reuse (1D, per Richtung)
# =========================================================

def _apply_variant_D(I_per_element, epsilon, eps_zero=1e-30):
    """Variante D: Proportional Slack-Reuse pro Richtung.

    Algorithmus:
        1. I_target = epsilon / N
        2. excess_e = max(0, I_e - I_target)    # Element ueberschreitet Budget
           slack_e  = max(0, I_target - I_e)    # Element hat Puffer
        3. Wenn Sum slack >= Sum excess: alle r_e = 1 (Slack reicht)
        4. Sonst: alpha = Sum slack / Sum excess (Reuse-Anteil)
           new_target_e = I_target + alpha * excess_e (variabel pro Element)
           r_e = max(1, sqrt(I_e / new_target_e))

    Args:
        I_per_element: list of I_e^(k) fuer eine Richtung
        epsilon: Fehlerbudget fuer DIESE Richtung (z.B. epsilon_global / n_dim)

    Returns:
        list of r_e (>= 1.0) per element
    """
    n = len(I_per_element)
    if n == 0:
        return []

    I_target = epsilon / n
    excess = [max(0.0, I - I_target) for I in I_per_element]
    slack = [max(0.0, I_target - I) for I in I_per_element]

    sum_excess = sum(excess)
    sum_slack = sum(slack)

    if sum_excess <= eps_zero or sum_excess <= sum_slack:
        return [1.0] * n

    alpha = sum_slack / sum_excess

    r_e = []
    for i in range(n):
        if excess[i] <= eps_zero:
            r_e.append(1.0)
        else:
            new_target = max(I_target + alpha * excess[i], eps_zero)
            r = math.sqrt(I_per_element[i] / new_target)
            r_e.append(max(1.0, r))
    return r_e


def _compute_h_ideal(h_lengths, r_factors):
    """h_ideal_k = h_k / r_k (kontinuierlich)."""
    result = []
    for h, r in zip(h_lengths, r_factors):
        if r > 0.0:
            result.append(h / r)
        else:
            result.append(h)
    return result


# =========================================================
# Element-Daten-Pipeline
# =========================================================

def compute_element_indicators(
    elements_corner_stresses,
    elements_corner_coords,
    element_volumes,
    m,
    epsilon,
    n_dim=3,
    sigma_eps=1e-12,
):
    """Berechnet pro Element alle 9 Felder fuer einen festen m-Wert.

    Args:
        elements_corner_stresses: list of [s1_n0, ..., s1_n7] pro Element (8 fuer Hex, 4 fuer Quad)
        elements_corner_coords:   list of [(x,y,z)_n0, ..., (x,y,z)_n7] pro Element
        element_volumes:          list of V_e
        m:        Weibull-Modul
        epsilon:  globales Fehlerbudget (z.B. 0.05)
        n_dim:    3 fuer Hex (LOC=V), 2 fuer Quad (LOC=A)

    Returns: dict mit
        E_e:          list (Skalar)
        W_e:          list (Skalar)
        I_e:          list (Skalar)
        r_e_xyz:      list of (r_x, r_y, r_z)   (z=1.0 bei Quad)
        h_ideal_xyz:  list of (hx, hy, hz)      (z=0.0 bei Quad)
        h_k_xyz:      list of (h_x, h_y, h_z)   Diagnose, original h-Werte
        delta_xyz:    list of (D_x, D_y, D_z)   Diagnose, Spannungs-Delta
        sigma_e:      list (Skalar) Diagnose
    """
    n = len(elements_corner_stresses)
    sigma_e_list = []
    deltas_list = []  # [(Dx, Dy, Dz), ...]
    hs_list = []

    # 1. Pro Element: sigma_e, deltas, h_lengths
    delta_h_fn = _compute_delta_and_h_hex if n_dim == 3 else _compute_delta_and_h_quad
    for i in range(n):
        cs = elements_corner_stresses[i]
        cc = elements_corner_coords[i]
        sigma_e_list.append(_compute_corner_mean(cs))
        d, h = delta_h_fn(cs, cc)
        deltas_list.append(d)
        hs_list.append(h)

    # 2. Hazard-Gewichte
    W_list = _compute_hazard_weights(sigma_e_list, element_volumes, m, sigma_eps)

    # 3. Pro Element: E_e^(k), I_e^(k)
    E_per_dir = []  # [(Ex, Ey, Ez), ...]
    I_per_dir = []  # [(Ix, Iy, Iz), ...]
    for i in range(n):
        E_k = _compute_E_per_direction(sigma_e_list[i], deltas_list[i], m, sigma_eps)
        E_per_dir.append(E_k)
        I_per_dir.append(tuple(W_list[i] * e for e in E_k))

    # Skalare Aggregate
    E_scalar = [sum(e) for e in E_per_dir]
    I_scalar = [sum(i) for i in I_per_dir]

    # 4. Variante D pro aktive Richtung
    epsilon_per_dir = epsilon / float(max(n_dim, 1))
    n_active = n_dim  # 3 oder 2
    r_e_components = []
    for k in range(3):
        if k < n_active:
            I_k_list = [I_per_dir[i][k] for i in range(n)]
            r_k_list = _apply_variant_D(I_k_list, epsilon_per_dir)
        else:
            # z-Achse bei Quad ist Dummy
            r_k_list = [1.0] * n
        r_e_components.append(r_k_list)

    # 5. h_ideal pro Richtung
    r_e_xyz = []
    h_ideal_xyz = []
    for i in range(n):
        r_vec = (r_e_components[0][i], r_e_components[1][i], r_e_components[2][i])
        r_e_xyz.append(r_vec)
        h_ideal_xyz.append(tuple(_compute_h_ideal(hs_list[i], r_vec)))

    # 6. Skalare Aggregate fuer r_e
    # r_e_vol = r_xi * r_eta * r_zeta: volumetrischer Refinement-Faktor
    #   (= Faktor um den Elementvolumen schrumpft; Baseline=1 = kein Refinement)
    # r_e_max = max(r_xi, r_eta, r_zeta): kritischste Einzelrichtung
    r_e_vol = [r[0] * r[1] * r[2] for r in r_e_xyz]
    r_e_max = [max(r[0], r[1], r[2]) for r in r_e_xyz]

    return {
        "E_e":       E_scalar,
        "W_e":       W_list,
        "I_e":       I_scalar,
        "E_per_dir": E_per_dir,    # [(E_xi, E_eta, E_zeta), ...] pro Element
        "I_per_dir": I_per_dir,    # [(I_xi, I_eta, I_zeta), ...] pro Element
        "r_e_xyz":   r_e_xyz,
        "r_e_vol":   r_e_vol,
        "r_e_max":   r_e_max,
        "h_ideal_xyz": h_ideal_xyz,
        "h_k_xyz":   hs_list,
        "delta_xyz": deltas_list,
        "sigma_e":   sigma_e_list,
    }


# =========================================================
# VTK-Writer
# =========================================================

def _write_vtk_refinement_file(
    filepath,
    points,
    cells,
    cell_type,  # 12 = Hex, 9 = Quad
    n_corners,  # 8 fuer Hex, 4 fuer Quad
    cell_scalars,  # dict {name: list of values}
    cell_vectors,  # dict {name: list of (vx, vy, vz)}
    title,
):
    """Schreibt VTK Legacy ASCII Unstructured Grid mit Cell-Vektor-Daten."""
    valid_cells = [c for c in cells if c is not None]
    valid_indices = [i for i, c in enumerate(cells) if c is not None]
    npts = len(points)
    ncells = len(valid_cells)

    with open(filepath, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"{title}\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n\n")

        # Points
        f.write(f"POINTS {npts} double\n")
        for x, y, z in points:
            f.write(f"{x:.8e} {y:.8e} {z:.8e}\n")
        f.write("\n")

        # Cells
        total_ints = ncells * (n_corners + 1)
        f.write(f"CELLS {ncells} {total_ints}\n")
        for cell in valid_cells:
            f.write(f"{n_corners} {' '.join(str(n) for n in cell)}\n")
        f.write("\n")

        # Cell Types
        f.write(f"CELL_TYPES {ncells}\n")
        for _ in range(ncells):
            f.write(f"{cell_type}\n")
        f.write("\n")

        # Cell Data
        f.write(f"CELL_DATA {ncells}\n")

        # Skalare
        for name, values in cell_scalars.items():
            f.write(f"SCALARS {name} double 1\n")
            f.write("LOOKUP_TABLE default\n")
            for i in valid_indices:
                f.write(f"{values[i]:.8e}\n")
            f.write("\n")

        # Vektoren (3 Komponenten)
        for name, values in cell_vectors.items():
            f.write(f"VECTORS {name} double\n")
            for i in valid_indices:
                vx, vy, vz = values[i]
                f.write(f"{vx:.8e} {vy:.8e} {vz:.8e}\n")
            f.write("\n")


# =========================================================
# Series-File (ParaView Time Series)
# =========================================================

def _write_refinement_series_file(vtk_dir, safe_name, m_list, suffix="_REFINE"):
    """Schreibt .vtk.series JSON fuer ParaView-Zeitserie ueber m-Werte."""
    series_dir = os.path.join(vtk_dir, "series")
    os.makedirs(series_dir, exist_ok=True)

    entries = [
        {"name": f"../{safe_name}{suffix}_m{m:02d}.vtk", "time": float(m)}
        for m in sorted(m_list)
    ]
    content = {"file-series-version": "1.0", "files": entries}
    series_path = os.path.join(series_dir, f"{safe_name}_refine.vtk.series")
    with open(series_path, 'w', encoding='utf-8') as f:
        json.dump(content, f, indent=2)
    return series_path


# =========================================================
# Geometrie + Stress Loader (LOC=V Hex-Pfad)
# =========================================================

def _read_sparse_nodes(filepath, ncols=6):
    """Identisch zu vtk_exporter._read_sparse_nodes — Sparse Array Reader."""
    nodes = {}
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            parts = line.split()
            if len(parts) < ncols:
                continue
            try:
                vals = [float(p) for p in parts[:ncols]]
            except ValueError:
                continue
            if all(abs(v) < 1e-30 for v in vals):
                continue
            nodes[line_num] = tuple(vals)
    return nodes


def _read_elements(filepath):
    """Identisch zu vtk_exporter._read_elements — Element-Konnektivitaet."""
    elements = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            try:
                vals = [int(round(float(p))) for p in parts]
            except ValueError:
                continue
            if len(vals) == 9:
                elements.append(tuple(vals))
            elif len(vals) == 8:
                elements.append((None,) + tuple(vals))
            elif len(vals) >= 9:
                elements.append(tuple(vals[:9]))
    return elements


def _read_volumes_from_v_out(filepath):
    """Liest {apdl}-V.out: {elem_id: volume} aus Format mit smax-Header + Elem/Vol-Spalten."""
    volumes = {}
    with open(filepath, 'r') as f:
        in_data = False
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if "Elem" in stripped and "Vol" in stripped:
                in_data = True
                continue
            if "smax" in stripped.lower():
                continue
            if not in_data:
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            try:
                eid = int(round(float(parts[0])))
                vol = float(parts[1])
                volumes[eid] = vol
            except ValueError:
                continue
    return volumes


def _read_corner_stresses_from_csv(csv_path):
    """Liest -S-clean.csv: dict {(elem_id, node_id): s1}.

    Format (Header): element,node,s1,s2,s3
    """
    stress_map = {}
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                eid = int(row['element'])
                nid = int(row['node'])
                s1 = float(row['s1'])
                stress_map[(eid, nid)] = s1
            except (KeyError, ValueError):
                continue
    return stress_map


def load_hex_data(tables_dir, apdl_name, stress_token):
    """Laedt fuer LOC=V (Hex8): corner_stresses, corner_coords, volumes, elem_ids.

    Args:
        tables_dir: Pfad zu .../tables/
        apdl_name:  ANSYS-Job-Name (gekuerzt fuer Dateinamen-Lookup)
        stress_token: "RAW" oder "AVG"

    Returns:
        dict mit:
            corner_stresses_per_elem: list of [s1_0, ..., s1_7]
            corner_coords_per_elem:   list of [(x,y,z)_0, ..., (x,y,z)_7]
            volumes:                  list of V_e
            element_ids:              list of ANSYS Element-IDs
            node_dict:                {nid: (x,y,z,...)}
            elem_list:                [(eid, n1..n8), ...]
    """
    # Geometrie aus VTK-Geometrie-Macro (immer in EM-Pfad vorhanden)
    nodes_file = os.path.join(tables_dir, f"{apdl_name}_VTK_Nodes.out")
    elem_file = os.path.join(tables_dir, f"{apdl_name}_VTK_Elemente.out")

    if not (os.path.exists(nodes_file) and os.path.exists(elem_file)):
        # Gauss-Pfad-Fallback: VEFF_*-Files
        nodes_file = os.path.join(tables_dir, f"VEFF_{apdl_name}_Nodes.out")
        elem_file = os.path.join(tables_dir, f"VEFF_{apdl_name}_Elemente.out")

    if not os.path.exists(nodes_file):
        raise FileNotFoundError(f"Nodes-Datei nicht gefunden: {nodes_file}")
    if not os.path.exists(elem_file):
        raise FileNotFoundError(f"Elemente-Datei nicht gefunden: {elem_file}")

    node_dict = _read_sparse_nodes(nodes_file, ncols=6)
    elem_list = _read_elements(elem_file)

    # Volumen
    vol_file = os.path.join(tables_dir, f"{apdl_name}-V.out")
    volumes_dict = {}
    if os.path.exists(vol_file):
        volumes_dict = _read_volumes_from_v_out(vol_file)

    # Stress-Quelle
    csv_file = os.path.join(tables_dir, f"{apdl_name}-S-clean.csv")
    stress_map = {}
    if os.path.exists(csv_file):
        stress_map = _read_corner_stresses_from_csv(csv_file)

    # Pro Element: 8 Eckknoten extrahieren
    corner_stresses = []
    corner_coords = []
    volumes = []
    element_ids = []

    for elem in elem_list:
        eid = elem[0]
        node_ids = list(elem[1:9])

        # Stresses: aus CSV (RAW) oder aus node_dict (AVG/Gauss)
        elem_stresses = []
        for nid in node_ids:
            if (eid, nid) in stress_map:
                elem_stresses.append(stress_map[(eid, nid)])
            elif nid in node_dict and len(node_dict[nid]) >= 4:
                elem_stresses.append(node_dict[nid][3])  # S1 ist Spalte 4
            else:
                elem_stresses.append(0.0)

        # Coords
        elem_coords = []
        for nid in node_ids:
            if nid in node_dict:
                elem_coords.append(node_dict[nid][:3])
            else:
                elem_coords.append((0.0, 0.0, 0.0))

        corner_stresses.append(elem_stresses)
        corner_coords.append(elem_coords)
        volumes.append(volumes_dict.get(eid, 0.0))
        element_ids.append(eid if eid is not None else len(element_ids) + 1)

    return {
        "corner_stresses": corner_stresses,
        "corner_coords": corner_coords,
        "volumes": volumes,
        "element_ids": element_ids,
        "node_dict": node_dict,
        "elem_list": elem_list,
    }


# =========================================================
# Geometrie + Stress Loader (LOC=A Quad-Pfad)
# =========================================================

def _read_a_out_faces(filepath):
    """Liest {apdl}-A.out: Surface-Face-Daten.

    Format:
        smax = ...
        face_id  parent  area  S1  S2  S3  n1  n2  n3  n4
        1.0      1.0     1.0   911 909 0.04 4.0 1.0 63.0 73.0
        ...

    Returns:
        list of dict: [{"face_id": int, "parent": int, "area": float,
                        "S1_face": float, "n_ids": (n1, n2, n3, n4)}, ...]
        S1_face ist 0.0 bei RAW-Cases (Stresses dort aus CSV).
    """
    faces = []
    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            # Header-Zeilen: "smax = ..." und "face_id  parent ..."
            if "smax" in stripped.lower() or "face_id" in stripped.lower():
                continue
            parts = stripped.split()
            if len(parts) < 10:
                continue
            try:
                face_id = int(round(float(parts[0])))
                parent = int(round(float(parts[1])))
                area = float(parts[2])
                s1_face = float(parts[3])
                n_ids = tuple(int(round(float(parts[i]))) for i in (6, 7, 8, 9))
            except ValueError:
                continue
            faces.append({
                "face_id": face_id,
                "parent": parent,
                "area": area,
                "S1_face": s1_face,
                "n_ids": n_ids,
            })
    return faces


def load_quad_data(tables_dir, apdl_name, stress_token):
    """Laedt fuer LOC=A (Quad4): corner_stresses, corner_coords, areas, face_ids.

    Args:
        tables_dir: Pfad zu .../tables/
        apdl_name:  ANSYS-Job-Name (gekuerzt fuer Dateinamen-Lookup)
        stress_token: "RAW" oder "AVG"

    Returns:
        dict mit:
            corner_stresses: list of [s1_n0, s1_n1, s1_n2, s1_n3]
            corner_coords:   list of [(x,y,z)_n0, ..., (x,y,z)_n3]
            areas:           list of A_e (Face-Flaeche)
            face_ids:        list of Face-IDs
            parent_elems:    list of Parent-Element-IDs
            node_dict:       {nid: (x, y, z, S1, ...)}
    """
    nodes_file = os.path.join(tables_dir, f"{apdl_name}_VTK_Nodes.out")
    if not os.path.exists(nodes_file):
        # Gauss-Pfad-Fallback
        nodes_file = os.path.join(tables_dir, f"VEFF_{apdl_name}_Nodes.out")
    if not os.path.exists(nodes_file):
        raise FileNotFoundError(f"Nodes-Datei nicht gefunden: {nodes_file}")

    a_file = os.path.join(tables_dir, f"{apdl_name}-A.out")
    if not os.path.exists(a_file):
        raise FileNotFoundError(f"-A.out nicht gefunden: {a_file}")

    node_dict = _read_sparse_nodes(nodes_file, ncols=6)
    faces = _read_a_out_faces(a_file)

    # Stress-Quelle (RAW: aus CSV; AVG: aus node_dict[nid][3])
    csv_file = os.path.join(tables_dir, f"{apdl_name}-A-S-clean.csv")
    stress_map = {}
    if os.path.exists(csv_file):
        stress_map = _read_corner_stresses_from_csv(csv_file)

    corner_stresses = []
    corner_coords = []
    areas = []
    face_ids = []
    parent_elems = []

    for face in faces:
        parent = face["parent"]
        n_ids = face["n_ids"]

        # 4 Knoten-Stresses: RAW aus CSV, AVG aus node_dict
        elem_stresses = []
        for nid in n_ids:
            if (parent, nid) in stress_map:
                elem_stresses.append(stress_map[(parent, nid)])
            elif nid in node_dict and len(node_dict[nid]) >= 4:
                elem_stresses.append(node_dict[nid][3])
            else:
                elem_stresses.append(0.0)

        # 4 Knoten-Koordinaten
        elem_coords = []
        for nid in n_ids:
            if nid in node_dict:
                elem_coords.append(node_dict[nid][:3])
            else:
                elem_coords.append((0.0, 0.0, 0.0))

        corner_stresses.append(elem_stresses)
        corner_coords.append(elem_coords)
        areas.append(face["area"])
        face_ids.append(face["face_id"])
        parent_elems.append(parent)

    return {
        "corner_stresses": corner_stresses,
        "corner_coords": corner_coords,
        "areas": areas,
        "face_ids": face_ids,
        "parent_elems": parent_elems,
        "node_dict": node_dict,
        "faces": faces,
    }


# =========================================================
# Top-Level-API
# =========================================================

def export_refinement_series(
    tables_dir,
    vtk_refinement_dir,
    case_name,
    apdl_name,
    stress_token,
    loc,
    epsilon=0.05,
    m_values=None,
):
    """Hauptfunktion: erzeugt VTK-Series mit Refinement-Indikatoren.

    Args:
        tables_dir:          {run}/{case}/tables/
        vtk_refinement_dir:  Output-Verzeichnis (wird angelegt)
        case_name:           voller Case-Name (fuer Dateinamen)
        apdl_name:           gekuerzter ANSYS-Job-Name (fuer File-Lookup)
        stress_token:        "RAW" oder "AVG"
        loc:                 "V" (Hex) oder "A" (Quad)
        epsilon:             globales Fehlerbudget (Default 0.05 = 5%)
        m_values:            list of m-Werte (Default: subset)

    Returns:
        list of geschriebenen VTK-Pfaden.
    """
    if m_values is None:
        m_values = [1, 5, 10, 15, 20, 25, 30, 40, 50]

    os.makedirs(vtk_refinement_dir, exist_ok=True)

    # ----- LOC=V (Hex-Cell, 8 Knoten, 3 Richtungen) -----
    if loc.upper() == "V":
        data = load_hex_data(tables_dir, apdl_name, stress_token)
        n_elem = len(data["element_ids"])
        print(f"   [VTK-Refine] LOC=V: {n_elem} Hex-Elemente geladen, eps={epsilon:.3f}")

        sorted_node_ids = sorted(data["node_dict"].keys())
        ansys_to_vtk = {aid: idx for idx, aid in enumerate(sorted_node_ids)}
        points = [data["node_dict"][aid][:3] for aid in sorted_node_ids]

        cells = []
        for elem in data["elem_list"]:
            node_ids = list(elem[1:9])
            try:
                cells.append([ansys_to_vtk[nid] for nid in node_ids])
            except KeyError:
                cells.append(None)

        cell_type = 12   # VTK_HEXAHEDRON
        n_corners = 8
        n_dim = 3
        sizes = data["volumes"]
        corner_stresses = data["corner_stresses"]
        corner_coords = data["corner_coords"]

    # ----- LOC=A (Quad-Cell, 4 Knoten, 2 Richtungen) -----
    elif loc.upper() == "A":
        data = load_quad_data(tables_dir, apdl_name, stress_token)
        n_face = len(data["face_ids"])
        print(f"   [VTK-Refine] LOC=A: {n_face} Quad-Faces geladen, eps={epsilon:.3f}")

        sorted_node_ids = sorted(data["node_dict"].keys())
        ansys_to_vtk = {aid: idx for idx, aid in enumerate(sorted_node_ids)}
        points = [data["node_dict"][aid][:3] for aid in sorted_node_ids]

        cells = []
        for face in data["faces"]:
            n_ids = face["n_ids"]
            try:
                cells.append([ansys_to_vtk[nid] for nid in n_ids])
            except KeyError:
                cells.append(None)

        cell_type = 9    # VTK_QUAD
        n_corners = 4
        n_dim = 2
        sizes = data["areas"]
        corner_stresses = data["corner_stresses"]
        corner_coords = data["corner_coords"]

    else:
        raise ValueError(f"LOC must be 'V' or 'A', got: {loc!r}")

    # Pro m-Wert: Indikatoren berechnen + VTK schreiben
    written_paths = []
    for m in m_values:
        result = compute_element_indicators(
            corner_stresses,
            corner_coords,
            sizes,
            float(m),
            float(epsilon),
            n_dim=n_dim,
        )

        def _comp(vec_list, k):
            return [v[k] for v in vec_list]

        cell_scalars = {
            # --- Aggregate (Skalare, Baseline: E/I/W >= 0, W: sum=1) ---
            "E":          result["E_e"],
            "W":          result["W_e"],
            "I":          result["I_e"],
            # --- Fehler + Indikator je Elementachse ---
            # LOC=V (Hex): xi/eta/zeta = 3 aktive Richtungen
            # LOC=A (Quad): xi/eta aktiv, zeta = 0.0 (Dummy)
            "E_xi":       _comp(result["E_per_dir"], 0),
            "E_eta":      _comp(result["E_per_dir"], 1),
            "E_zeta":     _comp(result["E_per_dir"], 2),
            "I_xi":       _comp(result["I_per_dir"], 0),
            "I_eta":      _comp(result["I_per_dir"], 1),
            "I_zeta":     _comp(result["I_per_dir"], 2),
            # --- Refinement-Faktoren je Elementachse ---
            # Baseline = 1.0 (kein Refinement noetig)
            # LOC=A: r_e_zeta = 1.0 (Dummy), r_e_vol = r_xi * r_eta
            "r_e_xi":     _comp(result["r_e_xyz"], 0),
            "r_e_eta":    _comp(result["r_e_xyz"], 1),
            "r_e_zeta":   _comp(result["r_e_xyz"], 2),
            # r_e_vol: volumetrischer Faktor r_xi*r_eta*r_zeta
            #   Baseline=1, Wert 8 = Element muss in 8 Sub-Elemente geteilt werden
            "r_e_vol":    result["r_e_vol"],
            # r_e_max: kritischste Einzelrichtung (konservative Abschaetzung)
            "r_e_max":    result["r_e_max"],
            # --- Ideal-Elementgroesse je Elementachse [mm] ---
            # LOC=A: h_ideal_zeta = 0.0 (Dummy)
            "h_ideal_xi":   _comp(result["h_ideal_xyz"], 0),
            "h_ideal_eta":  _comp(result["h_ideal_xyz"], 1),
            "h_ideal_zeta": _comp(result["h_ideal_xyz"], 2),
        }
        cell_vectors = {
            # VECTOR X/Y/Z entspricht Elementachsen xi/eta/zeta (NICHT globalem XYZ).
            # Bei strukturierten Meshes mit gekruemmter Topologie (z.B. PWH-Lochrand)
            # sind die Elementachsen gegenueber dem globalen Frame rotiert.
            # Magnitude = euklidische Norm (Raumdiagonale fuer h_ideal, kein
            # direktes physikalisches Aequivalent fuer r_e — r_e_vol nutzen).
            "r_e_local":     result["r_e_xyz"],
            "h_ideal_local": result["h_ideal_xyz"],
        }

        safe_name = case_name.replace(".", "_")
        fname = f"{safe_name}_REFINE_m{int(m):02d}.vtk"
        fpath = os.path.join(vtk_refinement_dir, fname)
        title = (
            f"Mesh Refinement Indicator (LOC={loc.upper()}) m={float(m):.1f}, "
            f"eps={epsilon:.3f} - {case_name} | "
            f"VECTOR X=xi(loc1) Y=eta(loc2) Z=zeta(loc3) are LOCAL element axes NOT global XYZ | "
            f"r_e_vol=r_xi*r_eta*r_zeta (volume factor, baseline=1)"
        )

        _write_vtk_refinement_file(
            fpath, points, cells,
            cell_type=cell_type, n_corners=n_corners,
            cell_scalars=cell_scalars,
            cell_vectors=cell_vectors,
            title=title,
        )
        written_paths.append(fpath)

    # Series-File
    series_path = _write_refinement_series_file(
        vtk_refinement_dir, case_name.replace(".", "_"), m_values
    )
    print(f"   [VTK-Refine] {len(written_paths)} VTK-Dateien + Series geschrieben")
    print(f"   [VTK-Refine] Series: {series_path}")

    return written_paths


# =========================================================
# CLI / Standalone Usage
# =========================================================

if __name__ == "__main__":
    print("vtk_refinement_exporter.py — Modul-Test")
    print("Nutzung: import und Aufruf von export_refinement_series()")
    print("Oder via 06-export-mesh-refinement.py (Standalone-Entry-Point)")
