# Architektur & Ablauf

Detaillierter End-to-End-Durchlauf der Pipeline — von der Fall-Definition bis zur
4-Faktor-Fehlerzerlegung. Ergänzt die [README.md](README.md) (Quickstart + Überblick).

---

## Pipeline auf einen Blick

```text
00-cases.csv
   │  (eine Zeile = ein Fall: Geometrie + Methoden-Baukasten STRESS/INT/REF/CRIT/LOC)
   ▼
01-run_manager.py ──────────────────────────────────────────────┐
   │  legt Batch-Ordner 05-RUNS/{timestamp}/ an,                 │
   │  schreibt pro Fall current_params.inp, startet ANSYS        │
   ▼                                                             │
ANSYS APDL (Batch)                                               │
   02-MAIN/02-main.inp                                           │
      │  liest current_params.inp + Geometrie-Weiche            │
      ▼                                                          │
   03-GEOMETRY/geom_{BEAM1|PWH|POR}.inp   → Modell + Vernetzung  │
      │                                                          │
      ▼  SOLVE (FEM-Lösung)                                      │
   04-POST/post_*.inp                     → Spannungs-Extraktion │
      │                                                          │
      ├─ EM-Pfad:    Element-/Face-Werte  → *.out                │
      └─ Gauss-Pfad: Fortran-Integration  → VEFF_*.out           │
   04-POST/FORTRAN/Effektives_Volumen_Unified.exe                │
      ▼                                                          │
06-POST_ANALYST/06-post_analyst.py (ResultAnalyst) ◄─────────────┘
   │  liest die .out-Dateien, berechnet Veff/Aeff + Pf über m,
   │  vergleicht gegen analytische Referenz
   ▼
05-RUNS/{timestamp}/{case}/tables/ANALYSIS_{case}.csv      ← Hauptausgabe
   │  (optional bei ERR_EXT=1, benötigt Gauss-26-Referenzen)
   ▼
ANALYSIS_{case}_extended.csv                               ← 4-Faktor-Fehlerzerlegung
```

---

## Schritt 1 — Batch- und Fall-Anlage (`01-run_manager.py`)

- Liest `00-cases.csv` (Semikolon-separiert). Jede Zeile ist ein Fall.
- Erzeugt einen Zeitstempel-Batch `05-RUNS/{YYYY-MM-DD_HH-MM}/` und darin pro Fall
  einen Unterordner `{CASE-GEOM-STRESS-INT-REF-CRIT-LOC-MESH}/`.
- Validiert die Methoden-Kombination zentral über
  [`method_naming.py`](06-POST_ANALYST/helpers/method_naming.py) (z. B. `EMX` nur mit
  `INT=EM`, `GPX` nur mit `INT=G*`).
- Schreibt pro Fall ein `current_params.inp` (Geometrie-Parameter, Netz, Methoden-Tokens
  als APDL-Variablen) und startet ANSYS im Batch-Modus.

## Schritt 2 — FEM-Lösung (ANSYS APDL)

[`02-MAIN/02-main.inp`](02-MAIN/02-main.inp) ist die Steuer-Datei:

1. Lädt `current_params.inp`.
2. **Geometrie-Weiche**: lädt anhand von `my_case` das passende
   [`03-GEOMETRY/geom_*.inp`](03-GEOMETRY/) (Modellaufbau, Material, Randbedingungen,
   Vernetzung mit SOLID186, Last).
3. `SOLVE` → FEM-Spannungsfeld.
4. Verzweigt in das passende Post-Processing in [`04-POST/`](04-POST/).

## Schritt 3 — Spannungs-Extraktion & Integration

Hier entscheiden **STRESS** (RAW/AVG) und **INT** (EM/G*) über vier Pfade.

### Elemental Mean (EM)

- **RAW** ([`post_EM_RAW.inp`](04-POST/post_EM_RAW.inp)): schreibt `{name}-V.out`
  (smax-Header + Element/Volumen) und `{name}-S.out` (ungemitteltes `PRESOL`-Listing).
  Python ([`parse_presol_output.py`](06-POST_ANALYST/stress_parsing/parse_presol_output.py))
  parst das PRESOL-Listing, bildet pro Element den Mittelwert der **element-lokalen**
  Knotenspannungen und integriert in Python.
- **AVG** ([`post_EM_NOD.inp`](04-POST/post_EM_NOD.inp)): nutzt die von ANSYS
  **nodal gemittelten** Spannungen.
- Für `LOC=A` (Fläche) gibt es die Surface-Varianten `post_EM_*_surf.inp`.

> EM rechnet die Integration vollständig in **Python** — kein Fortran.

### Gauss-Quadratur (G1…G9, G26)

Die Integration übernimmt das Fortran-Programm
[`Effektives_Volumen_Unified.exe`](04-POST/FORTRAN/) (Quellcode: `overhead_unified.f90` +
`effektivesVol_unified.f90`). Es liest **Knotenkoordinaten, Element-Konnektivität und
Spannungen**, interpoliert die Spannungen über die SOLID186-Formfunktionen auf die
Gauss-Punkte und integriert mit der Jacobi-Determinante. Zwei Aufruf-Wege:

- **AVG (NOD-Pfad)**: [`post_GAUSS_unified.inp`](04-POST/post_GAUSS_unified.inp) ruft das
  Makro `x_Effektives_Volumen_NOD.mac`, das die Eingabe-Dateien schreibt und die EXE direkt
  aus ANSYS via `/SYS` startet. Ergebnis: `VEFF_{name}.out` mit `Veff(m)`.
- **RAW**: [`post_GAUSS_RAW.inp`](04-POST/post_GAUSS_RAW.inp) exportiert nur Geometrie +
  ungemittelte `PRESOL`-Spannungen. **Python** (`_run_gauss_raw_fortran` in
  [`06-post_analyst.py`](06-POST_ANALYST/06-post_analyst.py)) baut daraus
  `effVol_RawStress.out` und ruft die EXE per `subprocess`.

Steuerung der EXE über `AVGModeRAW` (0=NOD/1=RAW) und `GaussPNorm` (0/1, für `REF=GPX`) in
`effVol_Parameter.out`. Die EXE berechnet immer alle vier Spalten (S1/PIA × V/A); Python
wählt anhand `CRIT`×`LOC` die passende aus.

## Schritt 4 — Auswertung (`ResultAnalyst`)

[`06-post_analyst.py`](06-POST_ANALYST/06-post_analyst.py) liest die `.out`-Dateien,
berechnet `Veff`/`Aeff` und `Pf` über das Weibull-Modul `m = 1…50`, bestimmt die drei
Spannungsreferenzen (`smax_nodal`, `smax_norm`, `smax_ref`) und vergleicht gegen die
**analytische Referenz** ([`analytical_helper.py`](06-POST_ANALYST/analytical/analytical_helper.py),
`LOADCASE_REGISTRY`). Hauptausgabe:
`tables/ANALYSIS_{case}.csv` plus JSON-Metadaten und ein Markdown-Report.

---

## Spezial-Workflow: 4-Faktor-Fehlerzerlegung (`ERR_EXT=1`)

Die Zerlegung `ρ_Total = ρ_Disk · ρ_Avg · ρ_Int · ρ_Ref` braucht eine
**exakte Integrations-Referenz** — und zwar in **beiden** Spannungsdarstellungen,
damit der Mittelungs-Faktor (`Avg`, RAW vs. AVG) vom Diskretisierungs-Faktor (`Disk`)
getrennt werden kann. Daher der **Dual-Gauss-26-Referenz-Workflow**:

### a) Referenzläufe rechnen (zweimal G26)

Pro `{CASE, GEOM, CRIT, LOC, MESH}` zwei G26-Läufe mit `REF=NDX`:

| STRESS | INT | REF | Zweck |
| --- | --- | --- | --- |
| `RAW` | `G26` | `NDX` | nahezu-exakte Integration des **ungemittelten** Feldes |
| `AVG` | `G26` | `NDX` | nahezu-exakte Integration des **nodal gemittelten** Feldes |

Vorlage: [`csv-templates/00-cases-g26-refs.csv`](csv-templates/00-cases-g26-refs.csv).

### b) Referenzen aggregieren

```bash
# in 06-generate_gauss26_reference.py die Variable SOURCE_RUN auf den G26-Batch-Ordner
# (05-RUNS/{timestamp}) setzen, dann:
python 06-POST_ANALYST/error_decomposition/06-generate_gauss26_reference.py
```

Das Skript liest alle `ANALYSIS_*.csv` des G26-Laufs (Filter: `INT=G26`, `REF=NDX`,
`STRESS ∈ {RAW, AVG}`) und schreibt die aggregierte
[`06-POST_ANALYST/REFERENCES/Reference_Gauss26.csv`](06-POST_ANALYST/REFERENCES/).

### c) Eigentliche Fälle mit `ERR_EXT=1` rechnen

Sobald für die gewünschten Lastfälle/Netze **beide** G26-Referenzen vorliegen, erzeugt jeder
Fall mit `ERR_EXT=1` zusätzlich `ANALYSIS_{case}_extended.csv` mit den 15 Spalten der
Zerlegung (5 Klassen × `ρ`/`δ`/`Δ`). Rechenkern:
[`calc_differentiated_errors.py`](06-POST_ANALYST/error_decomposition/calc_differentiated_errors.py).

> Fehlt eine passende G26-Referenz, wird der Fall klar als `DEGRADED` markiert
> (Basis-CSV bleibt erhalten, `_extended.csv` fehlt).

**Merge-Schlüssel:** `(Mesh, Loadcase, m)` mit
`Loadcase = {CASE}-{GEOM}-{CRIT}-{LOC}-{STRESS}` — INT/REF gehen bewusst **nicht** in den
Schlüssel ein, weil die Referenz immer G26/NDX ist.

---

## VTK-Export (optional, `VTK=1`)

Bei `VTK=1` erzeugt
[`vtk_exporter.py`](06-POST_ANALYST/output/vtk_exporter.py) pro Weibull-Modul eine
`vtk/{name}_RRI_m{mm}.vtk` (Legacy-ASCII, reines Python) mit dem normierten Risiko-Integral
(RRI) als `CELL_DATA` und den FEM-Hauptspannungen als `POINT_DATA` — für die Visualisierung
in ParaView.

> **VTK ist EM-only.** Im Gauss-Pfad liefert die Pipeline nur das integrierte `Veff(m)`,
> keine pro-Element-Aufschlüsselung — die Datenbasis für das RRI-Feld fehlt dort. Bei
> `INT=G*` wird der Export mit Hinweis übersprungen.

---

## Wichtigste Ausgabedateien

| Datei | Inhalt |
| --- | --- |
| `tables/ANALYSIS_{case}.csv` | Veff/Aeff + Pf über m, analytische Referenz, Basis-Fehler |
| `tables/ANALYSIS_{case}_extended.csv` | 4-Faktor-Fehlerzerlegung (nur mit `ERR_EXT=1` + G26-Refs) |
| `tables/{case}_REPORT.md` / `.json` | Markdown-Report + Metadaten |
| `vtk/{name}_RRI_m{mm}.vtk` | RRI-Feld pro m (nur EM, `VTK=1`) |
| `06-POST_ANALYST/REFERENCES/Reference_Gauss26.csv` | aggregierte Dual-G26-Referenz |
