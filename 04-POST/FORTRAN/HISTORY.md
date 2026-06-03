# Fortran Source-Code History

Historische Versions-Stempel der Fortran-Module, ausgelagert in v17.0 aus den
Source-Headern für bessere Lesbarkeit. Die *aktuellen* Code-Header verweisen nur
noch auf diese Datei.

## `overhead_unified.f90` (MODULE overhead_unified)

| Version | Datum | Beschreibung |
|---|---|---|
| 1.0 | 2016-02-24 | Erste lauffähige Version (dme) |
| 1.1 | 2016-03-03 | PIA-Option gelöscht — berechnet jetzt immer beide Kriterien (S1 + PIA) |
| 2.0 | 2016-03-09 | Arrayoperationen statt elementweiser Berechnung (Performance) |
| 3.0 | 2026-02-06 | GaussNorm-Variante (CK) — separates Modul `overhead_GaussNorm` |
| 4.0 | 2026-03-03 | RAWNodes: Dual-Mode NOD/RAW + integriertes GaussNorm (CK) |
| 5.0 | 2026-04-29 | **Unified** (aktiv): `overhead_unified` mit domain-spezifischer GaussNorm (CK) |

**v5.0 Unified-Erweiterungen** gegenüber RAWNodes (v4.0):

- Modulname `overhead_unified` statt `overhead` (kein `.mod`-Konflikt mit Legacy-Builds)
- Domain-spezifische GaussNorm-Maxima (statt globalem Maximum):
  - `max_sigma_ratio_volume_global`  (3D Volumen, gesetzt in `veff_gausz`)
  - `max_sigma_ratio_surface_global` (2D Fläche,  gesetzt in `seff_gausz`)
  - `max_sigma_ratio_line_global`    (1D Linie,   gesetzt in `leff_gausz`)
- Hauptprogramm korrigiert dann veff/seff/leff jeweils domain-spezifisch
  - `element_type=34` (3D Hex): veff ↔ VOL-Max, seff ↔ SURF-Max
  - `element_type=24` (2D Quad): veff ↔ SURF-Max, seff ↔ LINE-Max

**Bisherige Features** (übernommen aus RAWNodes v4.0):

- `stress_mode`: 0 = NOD (knotengemittelt wie v1.0+), 1 = RAW (elementlokale Spannungen aus PRESOL)
- `do_gaussnorm`: optionales GaussNorm-Tracking (aus `overhead_GaussNorm` v3.0)
- `nodes_coords(nnode, 3)`: nur Koordinaten für RAW-Modus
- `raw_stress(nelem_vol, eckn, 3)`: elementlokale Hauptspannungen (S1, S2, S3)
- `face_parent_elem(nface)`: Parent-Volumen-Element sequenzieller Index pro Face
- Gauss-Ordnungen 1-9 und 26 (tabelliert)

## `effektivesVol_unified.f90` (PROGRAM effektivesVol)

| Version | Datum | Beschreibung |
|---|---|---|
| 1.0 | 2016-02-24 | Erste lauffähige Version (dme) |
| 2.0 | 2016-03-09 | Arrayoperationen statt elementweiser Berechnung |
| 3.0 | 2026-02-06 | GaussNorm-Variante (CK) |
| 4.0 | 2026-03-03 | RAWNodes: Dual-Mode NOD/RAW + integriertes GaussNorm (CK) |
| 5.0 | 2026-04-29 | **Unified** (aktiv): domain-spezifische GaussNorm + erweiterter `GauszInfo`-Output (CK) |

**Methodische Änderungen v5.0** gegenüber GaussNorm-Variante (v3.0):

- veff-Spalten werden mit dem Maximum ihrer Integrationsdomäne korrigiert
  (3D `element_type=34`: VOL-Maximum; 2D `element_type=24`: SURF-Maximum)
- seff-Spalten analog mit dem Maximum ihrer Integrationsdomäne
  (3D: SURF-Maximum; 2D: LINE-Maximum)

**Output-Felder in `{ausgabename}_GauszInfo.out`** (v13.1 — nur domain-spezifisch):

- Metadaten: `SMAX_NODAL, GAUSS_FAKT, MMAX, SYMM, GEOM, BREITE`
- Domain-spezifische Maxima: `SMAX_GAUSS_VOL`, `SMAX_GAUSS_SURF`, `SMAX_GAUSS_LINE`
- Domain-spezifische Ratios: `RATIO_GAUSS_VOL`, `RATIO_GAUSS_SURF`, `RATIO_GAUSS_LINE`

## Legacy-Module (in `ARCHIV/`)

Diese Quelldateien sind nicht mehr aktiv, aber als Disaster-Recovery-Archiv erhalten:

| Datei | Letzter Stand | Ersetzt durch |
|---|---|---|
| `overhead.f90` | v2.0 (2016-03-09) | `overhead_unified` v5.0 |
| `overhead_GaussNorm.f90` | v3.0 (2026-02-06) | `overhead_unified` v5.0 |
| `overhead_highOrder.f90` | mit G26-Quadratur-Tabellen | `overhead_unified` v5.0 (G26 integriert) |
| `overhead_RAWNodes.f90` | v4.0 (2026-03-03) | `overhead_unified` v5.0 (Dual-Mode + Domain-spezifisch) |

**Wichtig**: alle Legacy-Module hatten `MODULE overhead` als Name → Modul-File-Konflikt
beim Build im selben Ordner. Unified hat `MODULE overhead_unified` und kann parallel
existieren.

## v19.0 PIAFix-Variante (parallele EXE)

| Version | Datum | Beschreibung |
|---|---|---|
| 6.0 | 2026-05-01 | **PIAFix** (parallel zu Legacy v5.0): PIA-Auswertung am Gauss-Punkt statt an Knoten (CK) |

### Methodische Differenz Legacy v5.0 vs. PIAFix v6.0

**Legacy** (`overhead_unified.f90` v5.0):

1. An Knoten: `sigma_eq_node = (S1^m + S2^m + S3^m)^(1/m)` (PIA-Aggregation mit Wurzel)
2. Interpoliere `sigma_eq_node` zu Gauss-Punkten via Form-Funktionen
3. Integriere `sigma_eq_gp^m * gram` (Wiederpotenzierung)

**PIAFix** (`overhead_unified_piafix.f90` v6.0):

1. Vor m-Loop: S1, S2, S3 **separat** zu Gauss-Punkten interpolieren
2. Im m-Loop, lokal am Gauss-Punkt:
   - `pia_integrand = MAX(s1, 0)^m + MAX(s2, 0)^m + MAX(s3, 0)^m`
   - `intLenV = pia_integrand * gram`
3. m=0: beide Pfade integrieren das geometrische Mass (User-Spec)

### Anwendungsfaelle

PIAFix ist methodisch konsistenter bei stark **mehrachsigen** Spannungsfeldern:

- **PWH-Lochrand**: biaxial (sigma_theta + sigma_r)
- **POR-Mittelpunkt**: equi-biaxial (sigma_r = sigma_phi)
- **3PB-Auflagerzone**: gemischtes Druck/Zug-Feld

Bei rein **uniaxialen** Cases (BEAM1/BEAM2) ist der Unterschied null (sigma2 = sigma3 = 0 → MAX(0,0)^m = 0).

### Numerische Differenz

Synthetischer Test (`tools/test_pia_fix_synthetic.py`, Two-Node-Line m=10):

- **Node A**: S1=2, S2=2, S3=0 (biaxial)
- **Node B**: S1=3, S2=0, S3=0 (uniaxial)
- **PIA-Legacy** = 1.834e+04, **PIA-Corrected** = 1.601e+04 → **−12.69 %**
- S1-Pfad: identisch in beiden Varianten (S2/S3 ignoriert)

### Limitierung — keine volle Tensor-Interpolation

Diese Variante nutzt **interpolierte Hauptspannungen** S1/S2/S3 — nicht die rohen Tensorkomponenten SX/SY/SZ/SXY/SYZ/SXZ. Das ist eine Vereinfachung, weil die Hauptspannungs-Achsen zwischen Knoten rotieren koennen (an Diskontinuitaeten besonders). Eine voellig saubere Loesung waere:

1. APDL exportiert die 6 Tensorkomponenten an Knoten
2. Fortran interpoliert die 6 Komponenten zum Gauss-Punkt
3. Fortran berechnet die lokalen Hauptspannungen via Eigenwert-Solver
4. PIA-Auswertung am Gauss-Punkt mit lokal berechneten S1/S2/S3

Das ist v20+ Scope — fuer v19.0 reicht der dominante PIA-Ordering-Fix.

### CSV-Schalter PIA_FIX

In `00-cases.csv` (Spalte 18):

- `PIA_FIX = 0` (oder leer): Legacy-EXE `Effektives_Volumen_Unified.exe`
- `PIA_FIX = 1`: PIAFix-EXE `Effektives_Volumen_Unified_PIAFix.exe`

Der Schalter wirkt nur bei `INT=G*` (Gauss-Pfad). Bei `INT=EM` ist er wirkungslos (silent ignore — EM hat keine Gauss-Punkte).

Beide Pfade (NOD via APDL `*IF, my_pia_fix, EQ, 1, ...` und RAW via Python `_run_gauss_raw_fortran()`) propagieren den Schalter.

### Build der PIAFix-EXE

Siehe [BUILD_UNIFIED.md](BUILD_UNIFIED.md) v19.0-Section.

---

## Querverweis

- Aktuelle Build-Anleitung: [BUILD_UNIFIED.md](BUILD_UNIFIED.md)
- Methodische Doku: [../../06-POST_ANALYST/ERROR_THEORY.md](../../06-POST_ANALYST/ERROR_THEORY.md)
- v19.0 Validation: [../../tools/test_pia_fix_synthetic.py](../../tools/test_pia_fix_synthetic.py)
