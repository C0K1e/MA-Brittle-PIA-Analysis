"""
method_naming.py — Single Source of Truth fuer Methoden-Baukasten-Tokens (v14.0)

Naming-Konvention nach Worknotes "Methoden-Baukasten" (Stress | Int | Ref):
  - Stress  (Spannungsdarstellung):  RAW | AVG       (BND/INR reserviert fuer v15+)
  - Int     (Integrationsverfahren):  EM | G1..G9 | G26
  - Ref     (Referenzspannung):       EMX | GPX | NDX

Kompatibilitaetsregeln:
  R1: REF=EMX nur mit INT=EM       (Element-Mean-Referenz)
  R2: REF=GPX nur mit INT=G*       (Gauss-Punkt-Referenz)
  REF=NDX passt mit allen INT      (Knoten-Referenz, ortsagnostisch)

Zusatzregel (LOC=A spezifisch):
  LOC=A + INT=EM + REF=GPX blockiert (EM-Pfad hat keine Gauss-Punkte)

Case-ID-Format (v14, 8 Tokens):
  {CASE}-{GEOM}-{STRESS}-{INT}-{REF}-{CRIT}-{LOC}-{MESH}

APDL-Bruecke (intern, NICHT Teil der oeffentlichen v14-API):
  STRESS=RAW <-> AVGModeRAW=1, my_avg='RAW'
  STRESS=AVG <-> AVGModeRAW=0, my_avg='NOD'
  REF=EMX    <-> my_norm='EM'
  REF=GPX    <-> my_norm='GP'
  REF=NDX    <-> my_norm='NOD'
"""

# ----------------------------------------------------------------------
# Token-Sets (Single Source of Truth)
# ----------------------------------------------------------------------

STRESS_TOKENS = {"RAW", "AVG"}
"""Aktive Stress-Tokens in v14. Werte fuer CSV-Spalte STRESS."""

STRESS_TOKENS_RESERVED = {"ABN", "AIN"}
"""Reserviert fuer v15+: Avg(BND) / Avg(INR). Aktuell numerisch nicht differenziert."""

INT_TOKENS = {"EM"} | {f"G{n}" for n in range(1, 10)} | {"G26"}
"""Aktive Integrations-Tokens. Werte fuer CSV-Spalte INT."""

REF_TOKENS = {"EMX", "GPX", "NDX"}
"""Aktive Ref-Tokens in v14. Werte fuer CSV-Spalte REF."""

CRIT_TOKENS = {"S1", "PIA"}
LOC_TOKENS = {"V", "A"}


# ----------------------------------------------------------------------
# APDL/Fortran-Bruecke (interne technische Uebersetzung)
# ----------------------------------------------------------------------

_STRESS_TO_LEGACY_AVG = {
    "RAW": "RAW",
    "AVG": "NOD",
}
_REF_TO_LEGACY_NORM = {
    "EMX": "EM",
    "GPX": "GP",
    "NDX": "NOD",
}

_LEGACY_AVG_TO_STRESS = {v: k for k, v in _STRESS_TO_LEGACY_AVG.items()}
_LEGACY_NORM_TO_REF = {v: k for k, v in _REF_TO_LEGACY_NORM.items()}


def legacy_avg_from_stress(stress: str) -> str:
    """STRESS-Token -> APDL my_avg-Wert.

    Brueckt zwischen v14-extern (STRESS=AVG) und v13.1-intern-APDL (my_avg='NOD').
    KEIN Legacy-Fallback fuer Eingangsdaten — strikte Validierung.
    """
    s = stress.strip().upper()
    if s in _STRESS_TO_LEGACY_AVG:
        return _STRESS_TO_LEGACY_AVG[s]
    if s in STRESS_TOKENS_RESERVED:
        raise NotImplementedError(
            f"STRESS={s} (BND/INR) ist in v14 reserviert, aber numerisch nicht implementiert. "
            f"Geplant fuer v15+ (Stress-Field-Postprocessing fuer Avg(BND)/Avg(INR))."
        )
    raise ValueError(f"Unbekanntes STRESS-Token: {s!r}. Erlaubt: {sorted(STRESS_TOKENS)}")


def legacy_norm_from_ref(ref: str) -> str:
    """REF-Token -> APDL my_norm-Wert.

    Brueckt zwischen v14-extern (REF=EMX/GPX/NDX) und v13.1-intern-APDL (my_norm='EM'/'GP'/'NOD').
    """
    r = ref.strip().upper()
    if r in _REF_TO_LEGACY_NORM:
        return _REF_TO_LEGACY_NORM[r]
    raise ValueError(f"Unbekanntes REF-Token: {r!r}. Erlaubt: {sorted(REF_TOKENS)}")


def stress_from_legacy_avg(avg: str) -> str:
    """APDL my_avg-Wert -> STRESS-Token (Reverse-Bruecke fuer Migration alter Daten)."""
    a = (avg or "").strip().upper()
    if a in _LEGACY_AVG_TO_STRESS:
        return _LEGACY_AVG_TO_STRESS[a]
    raise ValueError(f"Unbekanntes legacy AVG-Token: {a!r}. Erwartet: NOD oder RAW.")


def ref_from_legacy_norm(norm: str) -> str:
    """APDL my_norm-Wert -> REF-Token (Reverse-Bruecke fuer Migration alter Daten)."""
    n = (norm or "").strip().upper()
    if n in _LEGACY_NORM_TO_REF:
        return _LEGACY_NORM_TO_REF[n]
    raise ValueError(f"Unbekanntes legacy NORM-Token: {n!r}. Erwartet: EM, GP oder NOD.")


# ----------------------------------------------------------------------
# Validierung
# ----------------------------------------------------------------------

def validate_method_combination(stress: str, int_token: str, ref: str,
                                 loc: str = "V", crit: str = "PIA"):
    """Prueft Methoden-Kombination gegen Kompatibilitaetsregeln.

    Returns:
        (is_valid: bool, reason: str)

    Regeln:
        R1: REF=EMX nur mit INT=EM
        R2: REF=GPX nur mit INT=G*
        Zusatz: LOC=A + INT=EM + REF=GPX blockiert

    v15.3 (Audit-Find #3): CRIT wird jetzt auch validiert. Frueher konnte
    `CRIT=BLAH` silent durchlaufen und in `COL_MAP.get((LOC, CRIT), 2)`
    auf PIA-Spalte zurueckfallen ohne Fehlermeldung.
    """
    s = stress.strip().upper()
    i = int_token.strip().upper()
    r = ref.strip().upper()
    l = loc.strip().upper()
    c = crit.strip().upper()

    if s not in STRESS_TOKENS:
        if s in STRESS_TOKENS_RESERVED:
            return False, f"STRESS={s} reserviert fuer v15+ (BND/INR), nicht implementiert"
        return False, f"Unbekanntes STRESS={s!r}, erlaubt: {sorted(STRESS_TOKENS)}"
    if i not in INT_TOKENS:
        return False, f"Unbekanntes INT={i!r}, erlaubt: EM, G1..G9, G26"
    if r not in REF_TOKENS:
        return False, f"Unbekanntes REF={r!r}, erlaubt: {sorted(REF_TOKENS)}"
    if l not in LOC_TOKENS:
        return False, f"Unbekanntes LOC={l!r}, erlaubt: {sorted(LOC_TOKENS)}"
    # CRIT-Validierung
    if c not in CRIT_TOKENS:
        return False, f"Unbekanntes CRIT={c!r}, erlaubt: {sorted(CRIT_TOKENS)}"

    # R1: EMX nur mit EM
    if r == "EMX" and i != "EM":
        return False, f"R1 verletzt: REF=EMX nur mit INT=EM erlaubt (Element-Mean-Referenz braucht EM-Integration)"

    # R2: GPX nur mit Gauss
    if r == "GPX" and not i.startswith("G"):
        return False, f"R2 verletzt: REF=GPX nur mit INT=G* erlaubt (Gauss-Punkt-Referenz braucht Gauss-Integration)"

    # LOC=A spezifisch: EM hat keine Gauss-Punkte
    if l == "A" and i == "EM" and r == "GPX":
        return False, "LOC=A + INT=EM + REF=GPX nicht moeglich (EM-Pfad hat keine Gauss-Punkte)"

    return True, "OK"


# ----------------------------------------------------------------------
# Display-Labels
# ----------------------------------------------------------------------

STRESS_LABELS = {
    "RAW": "Stress[Raw]",
    "AVG": "Stress[Avg]",
    "ABN": "Stress[Avg(BND)]",
    "AIN": "Stress[Avg(INR)]",
}

REF_LABELS = {
    "EMX": "Ref[EMMax]",
    "GPX": "Ref[GPMax]",
    "NDX": "Ref[Nodal]",
}


def stress_label(stress: str) -> str:
    return STRESS_LABELS.get(stress.strip().upper(), stress)


def int_label(int_token: str) -> str:
    """INT-Token -> 'Int[EMean]' / 'Int[Gauss5]' / 'Int[Gauss26]'."""
    t = int_token.strip().upper()
    if t == "EM":
        return "Int[EMean]"
    if t.startswith("G"):
        return f"Int[Gauss{t[1:]}]"
    return f"Int[{t}]"


def ref_label(ref: str) -> str:
    return REF_LABELS.get(ref.strip().upper(), ref)


def method_code(stress: str, int_token: str, ref: str) -> str:
    """Kurzform fuer Plot-Filter etc.: 'RAW-G5-GPX'."""
    return f"{stress.strip().upper()}-{int_token.strip().upper()}-{ref.strip().upper()}"


def method_label(stress: str, int_token: str, ref: str) -> str:
    """Lange Form fuer Achsentitel: 'Stress[Raw] | Int[Gauss5] | Ref[GPMax]'."""
    return f"{stress_label(stress)} | {int_label(int_token)} | {ref_label(ref)}"


# ----------------------------------------------------------------------
# Case-ID Build / Parse
# ----------------------------------------------------------------------

def build_case_id(case: str, geom: str, stress: str, int_token: str, ref: str,
                  crit: str, loc: str, mesh: str) -> str:
    """Erzeugt v14-Case-ID: {CASE}-{GEOM}-{STRESS}-{INT}-{REF}-{CRIT}-{LOC}-{MESH}."""
    return f"{case}-{geom}-{stress}-{int_token}-{ref}-{crit}-{loc}-{mesh}"


def parse_case_id(case_id: str) -> dict:
    """Zerlegt v14-Case-ID in 8 Tokens.

    Erwartet: {CASE}-{GEOM}-{STRESS}-{INT}-{REF}-{CRIT}-{LOC}-{MESH}

    Cleanup:
        - "ANALYSIS_" Prefix wird entfernt
        - "_extended.csv"/".csv" Suffix wird entfernt

    Returns:
        dict mit Keys: case, geom, stress, int, ref, crit, loc, mesh, full_id, method_str
        Bei Parse-Fehler: dict mit Key 'error'

    method_str = "{STRESS}-{INT}-{REF}" — Lookup-Key fuer METHOD_COLORS etc.
    """
    clean = case_id.replace("ANALYSIS_", "").replace("_extended.csv", "").replace(".csv", "")
    parts = clean.split("-")

    if len(parts) < 8:
        return {"full_id": clean, "error": f"Invalid v14 case-id format: {len(parts)} tokens (expected 8)"}

    keys = ["case", "geom", "stress", "int", "ref", "crit", "loc", "mesh"]
    if len(parts) == 8:
        data = dict(zip(keys, parts))
    else:
        # Mehr als 8 Tokens — letzter ist Mesh, alles vorher zum Geom-Token zusammenfassen waere
        # ein Fehler. v14 erlaubt das nicht — Geometrie darf keine Bindestriche enthalten.
        # Fallback: erste 7 als keys, alles dazwischen verworfen.
        data = dict(zip(keys[:-1], parts[:7]))
        data["mesh"] = parts[-1]

    data["full_id"] = clean
    data["method_str"] = method_code(data["stress"], data["int"], data["ref"])
    return data
