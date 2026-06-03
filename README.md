# Probabilistische Festigkeitsbewertung spröder Werkstoffe (PIA-Ansatz)

Hybrid-Pipeline aus **Python** (Steuerung & Auswertung), **ANSYS APDL** (FEM) und
**Fortran** (numerische Gauss-Integration) zur Berechnung effektiver Volumina/Flächen
($V_\text{eff}$, $A_\text{eff}$) und Versagenswahrscheinlichkeiten ($P_f$) spröder
Werkstoffe nach dem **Weakest-Link-Modell** und dem **PIA-Kriterium** (Principle of
Independent Action).

Kern der Arbeit ist der **Vergleich verschiedener Berechnungsmethodiken** (Spannungs­darstellung,
Integrationsverfahren, Referenzspannung) und die **Zerlegung des Gesamtfehlers in einzelne,
physikalisch interpretierbare Faktoren**.

> **Kontext:** Dieses Repository ist im Rahmen der Masterarbeit
>
> **„Analyse und Weiterentwicklung von Berechnungsansätzen zur probabilistischen
> Festigkeitsbewertung spröder Materialien unter Verwendung der Finite-Elemente-Methode"**
>
> von **Conrad Kieselberger** an der **HTWK Leipzig** entstanden.
> Zur Herkunft des Fortran-Codes siehe [Herkunft & Attribution](#herkunft--attribution).

---

## Quickstart

### 1. Voraussetzungen

| Komponente | Version | Wofür |
| --- | --- | --- |
| **Python** | 3.11+ | Steuerung + Core-Auswertung (nur Standard-Library) |
| **ANSYS Mechanical** | v25x (Student oder höher) | FEM-Solver |
| **numpy + pandas** | beliebig | Fehlerzerlegung (`06-POST_ANALYST/error_decomposition/`) |
| *(Intel Fortran `ifx`)* | optional | nur falls die `.exe` neu kompiliert werden soll |

Die Fortran-Programme liegen **vorkompiliert** unter `04-POST/FORTRAN/*.exe` vor.

### 2. Lokale Konfiguration anlegen (der einzige Pflicht-Schritt)

Die rechner-spezifischen Pfade (ANSYS-EXE, Kernzahl) liegen **nicht** im Repo, sondern in
einer lokalen `config_local.py`. Einmalig aus der Vorlage erzeugen und anpassen:

```powershell
Copy-Item config_local.example.py config_local.py
# dann config_local.py öffnen und die ANSYS-EXE-Pfade für deinen Rechner eintragen
```

`config_local.py` ist in `.gitignore` und wird nie veröffentlicht.

### 3. Fälle definieren und starten

```bash
# Zu rechnende Fälle in 00-cases.csv definieren (Semikolon-separiert), dann:
python 01-run_manager.py
```

Der Run-Manager liest `00-cases.csv`, legt einen Zeitstempel-Ordner unter `05-RUNS/` an,
ruft für jeden Fall ANSYS auf und führt anschließend automatisch die Python-Auswertung aus.
Die Hauptausgabe liegt danach in `05-RUNS/{timestamp}/{case}/tables/ANALYSIS_*.csv`.

---

## Zielsetzung

Beim Weakest-Link-Modell wird das effektive Volumen

$$ V_\text{eff} = \int_V \left(\frac{\sigma(\vec{x})}{\sigma_\text{max}}\right)^m \, dV $$

numerisch ausgewertet ($m$ = Weibull-Modul). Je nachdem, **wie** die Spannungen aus der
FE-Lösung gewonnen, **wie** integriert und **worauf** normiert wird, ergeben sich
unterschiedliche Werte. Dieses Projekt

1. berechnet $V_\text{eff}$/$A_\text{eff}$ und $P_f$ für eine frei kombinierbare Menge an Methoden,
2. vergleicht sie gegen eine **analytische Referenz** (geschlossene Form bzw. Tabelle), und
3. **zerlegt den Gesamtfehler** in die Beiträge der einzelnen Methoden-Entscheidungen.

---

## Methoden-Baukasten

Jeder Fall in `00-cases.csv` kombiniert vier methodische Achsen (plus Geometrie & Last):

| Achse | Spalte | Werte | Bedeutung |
| --- | --- | --- | --- |
| **Spannungsdarstellung** | `STRESS` | `RAW`, `AVG` | element-lokale (ungemittelte) vs. nodal gemittelte Spannungen |
| **Integration** | `INT` | `EM`, `G1`…`G9`, `G26` | Elemental Mean vs. Gauss-Quadratur der Ordnung *n* |
| **Referenzspannung** | `REF` | `EMX`, `GPX`, `NDX` | Normierung auf Element-Maximum / Gauss-Punkt-Maximum / Knoten-Maximum |
| **Versagenskriterium** | `CRIT` | `PIA`, `S1` | alle Hauptspannungen (PIA) vs. nur $\sigma_1$ |
| **Domäne** | `LOC` | `V`, `A` | effektives Volumen vs. effektive Fläche |

Eine Methode wird eindeutig benannt als
`{CASE}-{GEOM}-{STRESS}-{INT}-{REF}-{CRIT}-{LOC}-{MESH}`, z. B.
`POR-50.8x0x1.5-RAW-G5-GPX-PIA-V-8x12x6`.

Die zulässigen Kombinationen (z. B. `EMX` nur mit `INT=EM`, `GPX` nur mit `INT=G*`) prüft
[`method_naming.py`](06-POST_ANALYST/helpers/method_naming.py) zentral — die *Single Source of
Truth* für Parsing, Validierung und Benennung.

---

## Fehlerzerlegung (4 Faktoren)

Der Gesamtfehler einer Methode gegenüber der analytischen Referenz wird **multiplikativ**
in vier Faktoren zerlegt:

$$ \rho_\text{Total} \;=\; \rho_\text{Disk} \cdot \rho_\text{Avg} \cdot \rho_\text{Int} \cdot \rho_\text{Ref} $$

| Faktor | Quelle des Fehlers | wird 1, wenn … |
| --- | --- | --- |
| **Disk** | Diskretisierung des FE-Spannungsfelds (FE-Linearisierung) | Netz → ∞ |
| **Avg** | nodale Mittelung der Spannungen | `STRESS=RAW` |
| **Int** | Integrationsverfahren (Quadratur-Ordnung) | `INT=G26` (Referenz-Ordnung) |
| **Ref** | Wahl der Referenz-/Normierungsspannung | bei $P_f$ exakt (Kürzung von $\sigma_\text{max}$ in der Weibull-Formel) |

mit $\rho_X = Q_x / Q_{x-1}$ (multiplikativer Faktor je Größe $Q \in \{V_\text{eff}, A_\text{eff}, P_f\}$).
So lässt sich pro Methode quantifizieren, **welcher Modellierungs-Schritt wie viel zum
Gesamtfehler beiträgt**. Implementierung:
[`calc_differentiated_errors.py`](06-POST_ANALYST/error_decomposition/calc_differentiated_errors.py)
(aktiviert per `ERR_EXT=1` in `00-cases.csv`; erzeugt `ANALYSIS_*_extended.csv`).

Die Zerlegung benötigt eine **Gauss-26-Referenz** (`Reference_Gauss26.csv`), die aus
G26-Läufen via [`06-generate_gauss26_reference.py`](06-POST_ANALYST/error_decomposition/06-generate_gauss26_reference.py)
aggregiert wird.

---

## Unterstützte Lastfälle

| Geometrie | Beschreibung | Analytik |
| --- | --- | --- |
| **BEAM1** | Biegebalken, reine Biegung (Moment via `SFGRAD`) | geschlossene Form |
| **PWH** | Zugversuch an quadratischer Lochplatte (Plate With Hole) | tabellenbasiert (diskrete Geometrie-Größen) |
| **POR** | Ringdruckversuch (Pressure on Ring), equi-biaxial | geschlossene Form |

> Frühere Geometrien (Kragarm, 3-Punkt-Biegung, rechteckige Lochplatte) sind in dieser
> öffentlichen Fassung nicht enthalten.

---

## Projektstruktur

```text
00-cases.csv                  Fall-Definitionen (Semikolon-separiert)
01-run_manager.py             Entry Point: Batch-Runs orchestrieren
config_local.example.py       Vorlage für lokale ANSYS-Pfade -> config_local.py
config_local.py               (lokal, gitignored) deine Maschinen-Pfade

02-MAIN/                      ANSYS-Haupt-Pipeline (.inp)
03-GEOMETRY/                  APDL-Geometrien (geom_BEAM1 / geom_PWH / geom_POR)
04-POST/                      Post-Processing
   ├─ post_*.inp / *.mac      APDL: EM- und Gauss-Auswertung
   └─ FORTRAN/                Unified-EXE + Quellcode (effektivesVol/overhead *.f90)
05-RUNS/                      Ergebnis-Ordner (gitignored)

06-POST_ANALYST/             Python-Auswertung
   ├─ 06-post_analyst.py      Hauptlogik (Klasse ResultAnalyst)
   ├─ analytical/             analytische Referenzwerte + LOADCASE_REGISTRY
   ├─ error_decomposition/    4-Faktor-Fehlerzerlegung + G26-Referenz
   ├─ stress_parsing/         PRESOL-Parser (RAW-Spannungen aus MAPDL)
   ├─ helpers/                method_naming (SSOT), pf_helper, …
   └─ output/                 Report- und VTK-Export

csv-templates/               Beispiel-CSVs für typische Studien
```

### Datenfluss (kurz)

```text
00-cases.csv → 01-run_manager.py → ANSYS (02-MAIN → 03-GEOMETRY, 04-POST)
            → Fortran-Integration (04-POST/FORTRAN, nur Gauss-Pfad)
            → 06-POST_ANALYST (Auswertung + Fehlerzerlegung)
            → 05-RUNS/{timestamp}/{case}/tables/ANALYSIS_*.csv
```

> **Detaillierter End-to-End-Ablauf** (FEM → Spannungsextraktion → EM-/Gauss-Integration →
> Fehlerzerlegung, inkl. des Dual-Gauss-26-Referenz-Workflows und VTK-Export):
> siehe **[ARCHITECTURE.md](ARCHITECTURE.md)**.

---

## Herkunft & Attribution

Der **Fortran-Kern** der numerischen Integration geht ursprünglich auf die Diplomarbeit von
**David Mevec** zurück:

> David Mevec. *„Auslegung einer Festigkeitsprüfung von Dentalkeramik mittels des B3B-Tests"*.
> Diplomarbeit, Montanuniversität Leoben, 2016.

Mevecs Original arbeitet mit **knoten-gemittelten Spannungen**. Im Rahmen dieser Masterarbeit
wurde der Fortran-Code methodisch erweitert:

- **Element-reine Spannungen** (`STRESS=RAW`): Auswertung mit ungemittelten, element-lokalen
  Spannungen zusätzlich zur nodal gemittelten Variante.
- **Referenz auf den maximalen Gauss-Punkt** (`REF=GPX`): Normierung wahlweise auf das Maximum
  der interpolierten Spannungen an den Gauss-Punkten statt auf den Knoten-Wert.
- **PIA-Interpolations-Korrektur** (`PIA_FIX`): Mevec bildet die PIA-Leitspannung
  (inkl. Druckspannungs­kürzung) am **Element­knoten** und interpoliert diesen Skalar auf die
  Gauss-Punkte. Hier wird stattdessen **jede Hauptspannung einzeln vom Knoten auf die
  Gauss-Punkte interpoliert** und erst **am Gauss-Punkt** die Druckspannungskürzung und die
  PIA-Summen­bildung durchgeführt — die physikalisch konsistentere Reihenfolge.

Die Python-Steuerung, die analytischen Referenzen, die Methoden-Validierung und die
4-Faktor-Fehlerzerlegung sind im Rahmen dieser Arbeit neu entstanden.

---

## Lizenz / Zitation

Dieses Repository entstand zu wissenschaftlichen Zwecken im Rahmen der oben genannten
Masterarbeit. Bei Nutzung bitte die Arbeit sowie die zugrunde liegende Diplomarbeit von
David Mevec zitieren.

*(Eine formale Open-Source-Lizenz ist noch festzulegen.)*
