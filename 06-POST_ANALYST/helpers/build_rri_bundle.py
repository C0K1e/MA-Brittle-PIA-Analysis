from __future__ import annotations

import argparse
#from importlib.resources import files
import re
import sys
from pathlib import Path

import numpy as np

try:
    import pyvista as pv
except ImportError:
    print("PyVista fehlt. Bitte installieren mit: pip install pyvista vtk")
    sys.exit(1)


RRI_CELL_ARRAYS = ["RRI", "RRI_norm"]
FIXED_POINT_ARRAYS = [
    "FEM_S1",
    "FEM_S2",
    "FEM_S3",
    "Analytical_S1",
    "Analytical_S2",
    "Analytical_S3",
]


def extract_m(path: Path) -> int:
    """
    Liest m aus Dateinamen wie ..._m02.vtk oder ..._m2.vtk
    """
    match = re.search(r"_m(\d+)\.vtk$", path.name, re.IGNORECASE)
    if not match:
        raise ValueError(f"Kein m-Wert im Dateinamen gefunden: {path.name}")
    return int(match.group(1))


def find_input_files(folder: Path, pattern: str) -> list[Path]:
    files = sorted(folder.glob(pattern))
    files = [f for f in files if f.is_file()]

    if not files:
        raise FileNotFoundError(
            f"Keine Dateien gefunden in {folder} mit Pattern {pattern}"
        )

    files = sorted(files, key=extract_m)
    return files


def compare_mesh_geometry_and_topology(mesh_ref: pv.DataSet, mesh_cur: pv.DataSet, file_name: str) -> None:
    """
    Prüft, ob Punkte, Zellen und Zelltypen identisch sind.
    """
    if mesh_ref.n_points != mesh_cur.n_points:
        raise ValueError(
            f"{file_name}: Anzahl Punkte unterschiedlich "
            f"({mesh_cur.n_points} statt {mesh_ref.n_points})"
        )

    if mesh_ref.n_cells != mesh_cur.n_cells:
        raise ValueError(
            f"{file_name}: Anzahl Zellen unterschiedlich "
            f"({mesh_cur.n_cells} statt {mesh_ref.n_cells})"
        )

    if not np.allclose(mesh_ref.points, mesh_cur.points, rtol=0.0, atol=1e-12):
        max_diff = np.max(np.abs(mesh_ref.points - mesh_cur.points))
        raise ValueError(
            f"{file_name}: Punktkoordinaten unterschiedlich, max. Abweichung = {max_diff}"
        )

    if hasattr(mesh_ref, "cells") and hasattr(mesh_cur, "cells"):
        if not np.array_equal(mesh_ref.cells, mesh_cur.cells):
            raise ValueError(f"{file_name}: Zell-Connectivity unterschiedlich")

    if hasattr(mesh_ref, "celltypes") and hasattr(mesh_cur, "celltypes"):
        if not np.array_equal(mesh_ref.celltypes, mesh_cur.celltypes):
            raise ValueError(f"{file_name}: Zelltypen unterschiedlich")


def compare_fixed_point_arrays(mesh_ref: pv.DataSet, mesh_cur: pv.DataSet, file_name: str) -> None:
    """
    Prüft, ob die festen Punktfelder über alle m-Dateien identisch bleiben.
    """
    for name in FIXED_POINT_ARRAYS:
        if name not in mesh_ref.point_data:
            raise KeyError(f"Referenzdatei enthält Punktarray '{name}' nicht")
        if name not in mesh_cur.point_data:
            raise KeyError(f"{file_name}: Punktarray '{name}' fehlt")

        arr_ref = np.asarray(mesh_ref.point_data[name])
        arr_cur = np.asarray(mesh_cur.point_data[name])

        if arr_ref.shape != arr_cur.shape:
            raise ValueError(
                f"{file_name}: Array '{name}' hat andere Form "
                f"{arr_cur.shape} statt {arr_ref.shape}"
            )

        if not np.allclose(arr_ref, arr_cur, rtol=0.0, atol=1e-10):
            max_diff = np.max(np.abs(arr_ref - arr_cur))
            raise ValueError(
                f"{file_name}: fixes Punktarray '{name}' ist nicht identisch, "
                f"max. Abweichung = {max_diff}"
            )


def build_multiarray_vtu(files: list[Path], output_vtu: Path) -> None:
    """
    Baut eine einzige VTU-Datei mit
    - einer Geometrie
    - festen Punktfeldern
    - variablen Zellfeldern RRI_mXX und RRI_norm_mXX
    """
    print("\n[1/2] Erzeuge Multi-Array-VTU ...")

    mesh_ref = pv.read(files[0])

    for name in RRI_CELL_ARRAYS:
        if name not in mesh_ref.cell_data:
            raise KeyError(f"Referenzdatei enthält Zellarray '{name}' nicht")

    for name in FIXED_POINT_ARRAYS:
        if name not in mesh_ref.point_data:
            raise KeyError(f"Referenzdatei enthält Punktarray '{name}' nicht")

    mesh_out = mesh_ref.copy(deep=True)

    # Die variablen Arrays aus der Referenz zunächst entfernen,
    # damit nur die m-spezifisch benannten Versionen drin sind.
    for name in RRI_CELL_ARRAYS:
        if name in mesh_out.cell_data:
            del mesh_out.cell_data[name]

    m_values = []

    for file in files:
        m = extract_m(file)
        print(f"  lese m={m:02d} aus {file.name}")
        mesh_cur = pv.read(file)

        compare_mesh_geometry_and_topology(mesh_ref, mesh_cur, file.name)
        compare_fixed_point_arrays(mesh_ref, mesh_cur, file.name)

        for name in RRI_CELL_ARRAYS:
            if name not in mesh_cur.cell_data:
                raise KeyError(f"{file.name}: Zellarray '{name}' fehlt")

        mesh_out.cell_data[f"RRI_m{m:02d}"] = np.asarray(mesh_cur.cell_data["RRI"]).copy()
        mesh_out.cell_data[f"RRI_norm_m{m:02d}"] = np.asarray(mesh_cur.cell_data["RRI_norm"]).copy()
        m_values.append(m)

    mesh_out.field_data["m_values"] = np.asarray(m_values, dtype=np.int32)

    output_vtu.parent.mkdir(parents=True, exist_ok=True)
    mesh_out.save(output_vtu)
    print(f"  gespeichert: {output_vtu}")


def write_series(files: list[Path], output_series: Path) -> None:
    """
    Schreibt eine JSON-.series-Datei für Legacy .vtk oder XML .vtu/.vtp Dateien.
    ParaView erkennt das als Zeitserie.
    """
    import json
    import os

    print("\n[2/2] Erzeuge SERIES-Zeitserie ...")

    output_series.parent.mkdir(parents=True, exist_ok=True)

    entries = []
    for file in files:
        m = extract_m(file)
        rel = os.path.relpath(file.resolve(), output_series.parent.resolve())
        rel = Path(rel).as_posix()
        entries.append({
            "name": rel,
            "time": float(m)
        })

    content = {
        "file-series-version": "1.0",
        "files": entries
    }

    output_series.write_text(
        json.dumps(content, indent=2),
        encoding="utf-8"
    )
    print(f"  gespeichert: {output_series}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Kombiniert RRI-VTK-Dateien zu Multi-Array-VTU und PVD-Zeitserie"
    )
    parser.add_argument(
        "folder",
        type=Path,
        help="Ordner mit den .vtk Dateien",
    )
    parser.add_argument(
        "--pattern",
        default="*_m*.vtk",
        help="Glob-Pattern für die Eingabedateien, Standard: *_m*.vtk",
    )
    parser.add_argument(
        "--prefix",
        default="RRI_bundle",
        help="Präfix für Ausgabedateien",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Ausgabeordner, Standard: <folder>/combined",
    )

    args = parser.parse_args()

    folder = args.folder.resolve()
    outdir = args.outdir.resolve() if args.outdir else folder / "combined"

    files = find_input_files(folder, args.pattern)

    print("Gefundene Dateien:")
    for f in files:
        print(f"  {f.name}")

    output_vtu = outdir / f"{args.prefix}_multiarray.vtu"
    output_series = outdir / f"{args.prefix}.vtk.series"

    write_series(files, output_series)

    print("\nFertig.")
    print(f"Multi-Array Datei : {output_vtu}")
    print(f"Series Datei      : {output_series}")
    print("\nIn ParaView:")
    print("  - .vtk.series öffnen, wenn m wie eine Zeitachse durchgespielt werden soll")
    print("  - .vtu öffnen, wenn du direkt zwischen Arrays RRI_m01 ... RRI_m50 umschalten willst")


if __name__ == "__main__":
    main()