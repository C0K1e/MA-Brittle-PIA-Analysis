# analytical_helper.py
# Analytische Referenzloesungen fuer Veff, Aeff, Smax (Weakest Link / PIA)
#
# LOADCASE_REGISTRY — zentrales Case-Dictionary.
# Geschlossene Formeln fuer BEAM1, BEAM2, 3PB, POR.
# Tabellenbasierte Loesungen (PWH, PWHR) via Getter aus analytical_reference_data.py.

import math

from analytical_reference_data import (
    get_reference_value,
    available_case_x_keys,
)


# =========================================================
# POR (Pressure on Ring) — Analytische Spannungsfelder & Aeff
# =========================================================
# Notation: R_s = Auflagerring-Radius, R_d = Scheiben-Aussenradius,
#           t = Plattendicke, p = Druck, nu = Poisson-Zahl
# Spannungen auf der Zugseite (z=0), r = Radialkoordinate
#
# POR-Parameter-Konvention (geom_POR.inp v22.x):
#   case_x = D_out_base (Basis-Aussendurchmesser, typisch 50.8 mm)
#   case_y = Skalierungsfaktor fuer R_s, R_d, L_sq (1.0 = unskaliert)
#            Bei case_y > 1: alle Laengen x case_y, P_load / case_y^2
#            -> sigma_max bleibt mathematisch invariant (Plattenformel sigma_max ~ p*(R/t)^2)
#            -> Aeff, Veff, Vtot, Atot skalieren mit case_y^2
#   case_z = T_plate (Plattendicke, NICHT skaliert)
# Rueckwaerts-Kompatibilitaet: case_y = 0 wird als case_y = 1.0 interpretiert.

def _por_scale(dim_y):
    """case_y aus CSV -> Skalierungsfaktor.

    Rueckwaerts-Kompat: case_y = 0 (alte CSVs) -> 1.0 (unskaliert).
    case_y > 0 -> echter Skalierungsfaktor.
    """
    return float(dim_y) if dim_y and float(dim_y) > 0.0 else 1.0

def _por_sigma_max(p, R_s, R_d, t, nu):
    """Analytische Maximalspannung bei r=0 auf der Zugseite (equi-biaxial).

    sigma_max = (3p*R_s^2)/(8t^2) * [(1-nu)*R_s^2/R_d^2 + 2(1+nu)]
                + p*(3+nu) / (4*(1-nu))
    """
    bracket = (1 - nu) * R_s**2 / R_d**2 + 2 * (1 + nu)
    return (3 * p * R_s**2) / (8 * t**2) * bracket + p * (3 + nu) / (4 * (1 - nu))


def _por_sigma_rr(r, p, R_s, R_d, t, nu):
    """Radialspannung sigma_rr(r) auf der Zugseite (z=0)."""
    bracket = ((1 - nu) * R_s**2 / R_d**2
               + 2 * (1 + nu)
               - (3 + nu) * r**2 / R_s**2)
    return (3 * p * R_s**2) / (8 * t**2) * bracket + p * (3 + nu) / (4 * (1 - nu))


def _por_sigma_phiphi(r, p, R_s, R_d, t, nu):
    """Tangentialspannung sigma_phiphi(r) auf der Zugseite (z=0)."""
    bracket = ((1 - nu) * R_s**2 / R_d**2
               + 2 * (1 + nu)
               - (1 + 3 * nu) * r**2 / R_s**2)
    return (3 * p * R_s**2) / (8 * t**2) * bracket + p * (3 + nu) / (4 * (1 - nu))


def _compute_aeff_por(m, R_s, R_d, nu):
    """Effektive Flaeche A_e,POR,full (volle Formel, NICHT ASTM-Kurzform).

    A_e = (pi/(m+1)) * (R_s/R_d)^2
          * {4(1+nu)*A - (3+nu)*(1-nu)^(m+1)*(R_d^2+R_s^2)^(m+1) / A^m}
          / ((3+nu)*(1+3*nu))

    wobei A = 2*R_d^2*(1+nu) + R_s^2*(1-nu)

    Args:
        m: Weibull-Modul (int oder float, >= 1)
        R_s: Auflagerring-Radius [mm]
        R_d: Scheiben-Aussenradius [mm]
        nu: Poisson-Zahl

    Returns:
        float: Effektive Flaeche [mm^2]
    """
    m_f = float(m)
    A = 2 * R_d**2 * (1 + nu) + R_s**2 * (1 - nu)

    term1 = 4 * (1 + nu) * A
    term2 = ((3 + nu) * (1 - nu)**(m_f + 1)
             * (R_d**2 + R_s**2)**(m_f + 1) / A**m_f)

    numerator = term1 - term2
    denominator = (3 + nu) * (1 + 3 * nu)

    return (math.pi / (m_f + 1)) * (R_s / R_d)**2 * numerator / denominator


# =========================================================
# Case-spezifische Veff-Funktionen
# =========================================================
# Einheitliche Signatur: (m, vol_base, dim_x, dim_y, dim_z, crit)

def _beam1_veff(m, vol_base, dim_x, dim_y, dim_z, crit):
    """BEAM1 / 3PB: Veff = V / (2*(m+1)). Identisch fuer uniaxiale Faelle mit nu=0."""
    return vol_base / (2.0 * (float(m) + 1.0))


def _beam2_veff(m, vol_base, dim_x, dim_y, dim_z, crit):
    """BEAM2 (Kragarm): Veff = V / (2*(m+1)*(2m+1))."""
    m_f = float(m)
    return vol_base / (2.0 * (m_f + 1.0) * (2.0 * m_f + 1.0))


def _pwh_veff(m, vol_base, dim_x, dim_y, dim_z, crit):
    """PWH Veff via Referenztabelle.

    v15.3 (Audit-Find): Hart fail bei fehlender exakter Tabelle. Frueher gab es
    einen 3-stufigen Fallback (S1→PIA, naechster CASE_X-Key) der wissenschaftlich
    bequem war, aber Methodenname und Datenquelle divergieren liess. Bei einer
    Master-Thesis muss klar sein: entweder exakte Tabelle vorhanden oder Fehler.
    """
    case_x_key = int(dim_x) if dim_x else 0
    crit_upper = crit.upper() if crit else "PIA"

    val = get_reference_value('Veff', 'PWH', crit_upper, case_x_key, m)
    if val is not None:
        return val

    # Hart fail mit klarer Diagnose-Information
    avail = available_case_x_keys('Veff', 'PWH', crit_upper)
    raise RuntimeError(
        f"PWH Veff-Tabelle fehlt fuer CASE_X={case_x_key}, CRIT={crit_upper}, m={m}.\n"
        f"Verfuegbare CASE_X-Keys fuer CRIT={crit_upper}: {sorted(avail) if avail else '[KEINE]'}\n"
        f"Aktion: (1) genaue CASE_X-Tabelle in analytical_reference_data.py erzeugen, "
        f"oder (2) CASE_X aus 00-cases.csv auf einen verfuegbaren Key anpassen.\n"
        f"Kein Fallback auf naechsten Key — wuerde Methodenname und Datenquelle "
        f"divergieren lassen (v15.3-Strikt-Pipeline)."
    )


def _pwhr_veff(m, vol_base, dim_x, dim_y, dim_z, crit):
    """PWHR Veff via Referenztabelle.

    v15.3 (Audit-Find): Hart fail bei fehlender PWHR-Tabelle. Frueher Fallback
    auf PWH (rechteckig != quadratisch — methodisch ungenau), jetzt explizit.
    """
    case_x_key = int(dim_x) if dim_x else 0
    crit_upper = crit.upper() if crit else "PIA"

    val = get_reference_value('Veff', 'PWHR', crit_upper, case_x_key, m)
    if val is not None:
        return val

    # Hart fail
    raise RuntimeError(
        f"PWHR Veff-Tabelle fehlt fuer CASE_X={case_x_key}, CRIT={crit_upper}, m={m}.\n"
        f"PWHR-Daten sind aktuell nicht in analytical_reference_data.py vorhanden.\n"
        f"Aktion: (1) PWHR-Tabelle erzeugen, oder (2) CASE auf PWH umstellen "
        f"(quadratische Lochplatte) wenn die Lastfall-Geometrie das erlaubt.\n"
        f"Kein Fallback auf PWH — rechteckige Lochplatte hat andere Spannungsverteilung "
        f"als quadratische (v15.3-Strikt-Pipeline)."
    )


def _compute_aeff_por_s1(m, R_s, R_d, nu):
    """POR Aeff fuer CRIT=S1: integriert nur sigma_phiphi (groesste positive Hauptspannung).

    v20.0 (Bug #1): _por_aeff hatte den crit-Parameter ignoriert und immer die
    PIA-Summe genutzt. Bei CRIT=S1 muss aber nur die groesste positive
    Hauptspannung integriert werden. Auf der POR-Zugseite (z=0, r in [0, R_s])
    gilt sigma_phiphi >= sigma_rr bei nu < 1 (Koeffizient 1+3nu vor r^2/R_s^2 ist
    kleiner als 3+nu in sigma_rr -> sigma_phiphi faellt langsamer ab).

    Approximation: Membran-Term p(3+nu)/(4(1-nu)) wird vernachlaessigt
    (typisch << Plattenbiegungs-Term P bei duennen Platten R_s^2/t^2 >> 1).
    Bei dieser Approximation kuerzen sich p und t aus dem Verhaeltnis
    sigma_phiphi(r)/sigma_max heraus:

        sigma_phiphi(r) / sigma_max = 1 - alpha * (r/R_s)^2
        alpha = (1+3*nu) / ((1-nu) * (R_s/R_d)^2 + 2*(1+nu))

    Geschlossene Form ueber Substitution y = alpha * (r/R_s)^2:
        Aeff_S1 = pi * R_s^2 / alpha * (1 - (1-alpha)^(m+1)) / (m+1)

    Args:
        m: Weibull-Modul (>=1)
        R_s: Auflagerring-Radius [mm]
        R_d: Scheiben-Aussenradius [mm]
        nu: Poisson-Zahl

    Returns:
        float: Effektive Flaeche Aeff_S1 [mm^2]
    """
    m_f = float(m)
    alpha = (1.0 + 3.0 * nu) / ((1.0 - nu) * (R_s / R_d)**2 + 2.0 * (1.0 + nu))
    return math.pi * R_s**2 / alpha * (1.0 - (1.0 - alpha)**(m_f + 1.0)) / (m_f + 1.0)


def _por_veff(m, vol_base, dim_x, dim_y, dim_z, crit):
    """POR: Veff = Aeff_full * t / (2*(m+1)).

    v20.0 (Bug #1): crit-Parameter wird jetzt respektiert.
    CRIT=PIA: bestehende Volle-Formel via _compute_aeff_por.
    CRIT=S1: aktuell blockiert (NotImplementedError) bis durch-die-Dicke-Integration
             mit Membran-Term sauber abgeleitet ist. Workaround: LOC=A nutzen.
    """
    crit_upper = crit.upper() if crit else "PIA"
    if crit_upper == "S1":
        raise NotImplementedError(
            "POR + CRIT=S1 + LOC=V (Veff): durch-die-Dicke-Integration mit Membran-"
            "Term braucht saubere Ableitung. Aeff(POR+S1) ist verfuegbar (LOC=A); "
            "Veff(POR+S1) deferred bis Formel saniert ist.\n"
            "Workaround: (1) LOC=A nutzen oder (2) CRIT=PIA verwenden."
        )
    m_f = float(m)
    scale = _por_scale(dim_y)          # case_y aus CSV (case_y = 0 -> 1.0)
    R_d = (dim_x / 2.0) * scale         # case_x/2 = D_out_base/2, skaliert
    t = dim_z                            # case_z = T_plate (NICHT skaliert)
    R_s = 22.73 * scale                  # Auflagerring-Radius (Fixture-Basis), skaliert
    nu = 0.22                            # Poisson-Zahl (wie in geom_POR.inp)
    # PIA-Pfad (unveraendert; R_s und R_d gemeinsam skaliert -> R_s/R_d invariant)
    Aeff = _compute_aeff_por(m_f, R_s, R_d, nu)
    return Aeff * t / (2.0 * (m_f + 1.0))


# =========================================================
# Case-spezifische Smax-Referenz-Funktionen
# =========================================================
# Einheitliche Signatur: (smax_nodal, load_n, dim_x, dim_y, dim_z)

def _beam1_smax(smax_nodal, load_n, dim_x, dim_y, dim_z):
    """BEAM1: sigma = 6*M/(b*h^2), M = 1000*load_F."""
    if dim_y > 0 and dim_z > 0:
        return 6000.0 * load_n / (dim_z * dim_y**2)
    return smax_nodal


def _beam2_smax(smax_nodal, load_n, dim_x, dim_y, dim_z):
    """BEAM2: target_sigma=1500 (hardcoded in geom_BEAM2.inp)."""
    return 1500.0


def _three_pb_smax(smax_nodal, load_n, dim_x, dim_y, dim_z):
    """3PB: sigma = 3*F*L_span / (2*b*h^2)."""
    if dim_x > 0 and dim_y > 0 and dim_z > 0:
        return 3.0 * load_n * dim_x / (2.0 * dim_z * dim_y**2)
    return smax_nodal


def _pwh_smax(smax_nodal, load_n, dim_x, dim_y, dim_z):
    """PWH/PWHR: K_t = 3 fuer kleine Bohrung in unendlicher Platte."""
    return 3.0 * load_n


def _por_smax(smax_nodal, load_n, dim_x, dim_y, dim_z):
    """POR: sigma_max bei r=0 auf der Zugseite (equi-biaxial).

    case_y-Skalierung: R_s, R_d skalieren mit dim_y; load_n / dim_y^2 da
    geom_POR.inp den effektiven Druck als P_load = load_F / case_y^2 setzt.
    Damit bleibt sigma_max mathematisch invariant unter Skalierung
    (sigma_max ~ p * (R/t)^2; (p/scale^2) * (scale*R/t)^2 = p * (R/t)^2).
    """
    scale = _por_scale(dim_y)
    R_d = (dim_x / 2.0) * scale if dim_x > 0 else 25.4 * scale
    R_s = 22.73 * scale
    t = dim_z if dim_z > 0 else 1.5
    p_eff = load_n / (scale * scale)         # entspricht P_load im APDL
    return _por_sigma_max(p_eff, R_s, R_d, t, 0.22)


# =========================================================
# Case-spezifische Aeff-Funktionen (Effektive Flaeche)
# =========================================================
# Einheitliche Signatur: (m, area_base, dim_x, dim_y, dim_z, crit)

def _beam1_aeff(m, area_base, dim_x, dim_y, dim_z, crit):
    """BEAM1: Konstantes Moment -> sigma uniform auf Zugflaeche -> Aeff = A_total."""
    return area_base


def _beam2_aeff(m, area_base, dim_x, dim_y, dim_z, crit):
    """BEAM2: Quadratisches Moment sigma/smax=(x/L)^2 -> Aeff = A_total/(2m+1)."""
    return area_base / (2.0 * float(m) + 1.0)


def _por_aeff(m, area_base, dim_x, dim_y, dim_z, crit):
    """POR Aeff: crit-aware Dispatch (v20.0 Bug #1) mit case_y-Skalierung (v22.4).

    CRIT=PIA: _compute_aeff_por (volle Formel ueber sigma_rr^m + sigma_phiphi^m)
    CRIT=S1:  _compute_aeff_por_s1 (nur sigma_phiphi^m, S1-Approximation ohne Membran-Term)

    Skalierung: R_s und R_d gemeinsam x case_y -> Aeff skaliert mit case_y^2
    (weil Vorfaktor R_s^2 in der Plattenformel).
    """
    scale = _por_scale(dim_y)
    R_d = (dim_x / 2.0) * scale
    R_s = 22.73 * scale
    nu = 0.22
    crit_upper = crit.upper() if crit else "PIA"
    if crit_upper == "S1":
        return _compute_aeff_por_s1(float(m), R_s, R_d, nu)
    elif crit_upper == "PIA":
        return _compute_aeff_por(float(m), R_s, R_d, nu)
    else:
        raise ValueError(f"POR aeff: unknown CRIT={crit_upper!r} (erwartet 'PIA' oder 'S1')")


# =========================================================
# LOADCASE_REGISTRY — Zentrales Case-Dictionary
# =========================================================
# Jeder Eintrag definiert alle case-spezifischen Eigenschaften.
# Neue Lastfaelle: Eintrag hier + CASE_DESCRIPTIONS (report_generator.py)
#                  + CASE_LABELS (case_id_legend.py).

LOADCASE_REGISTRY = {
    "BEAM1": {
        "symmetry_factor": 1.0,
        "vtot_exact": lambda cx, cy, cz: cx * cy * cz,
        "atot_exact": lambda cx, cy, cz: cx * cz,           # Zugflaeche Y=H: L*T
        "veff": _beam1_veff,
        "aeff": _beam1_aeff,                                  # = A_total (sigma uniform)
        "smax_ref": _beam1_smax,
        "constants": {},
    },
    "BEAM2": {
        "symmetry_factor": 1.0,
        "vtot_exact": lambda cx, cy, cz: cx * cy * cz,
        "atot_exact": lambda cx, cy, cz: cx * cz,           # Zugflaeche Y=H: L*T
        "veff": _beam2_veff,
        "aeff": _beam2_aeff,                                  # = A_total / (2m+1)
        "smax_ref": _beam2_smax,
        "constants": {"target_stress": 1500.0},
    },
    "PWH": {
        "symmetry_factor": 4.0,
        "vtot_exact": lambda cx, cy, cz: (cx**2 - math.pi * cz**2) * 1.0,
        "atot_exact": lambda cx, cy, cz: 2.0 * (cx**2 - math.pi * cz**2),
        "veff": _pwh_veff,
        "aeff": None,                                          # TODO: Tabelle (biaxiales Feld)
        "smax_ref": _pwh_smax,
        "constants": {"K_t": 3.0, "T_plate": 1.0},
    },
    "PWHR": {
        "symmetry_factor": 4.0,
        "vtot_exact": lambda cx, cy, cz: (cx * cy - math.pi * cz**2) * 1.0,
        "atot_exact": lambda cx, cy, cz: 2.0 * (cx * cy - math.pi * cz**2),
        "veff": _pwhr_veff,
        "aeff": None,                                          # TODO: Tabelle (wie PWH)
        "smax_ref": _pwh_smax,
        "constants": {"K_t": 3.0, "T_plate": 1.0},
    },
    "3PB": {
        "symmetry_factor": 4.0,
        "vtot_exact": lambda cx, cy, cz: cx * cy * cz,
        "atot_exact": None,                                    # TODO: Komplex (Auflager-Ausschluss)
        "veff": _beam1_veff,                                   # Identisch mit BEAM1 (nu=0)
        "aeff": None,                                          # TODO: Lineares Moment + Ausschlusszone
        "smax_ref": _three_pb_smax,
        "constants": {},
    },
    "POR": {
        "symmetry_factor": 4.0,
        # case_y-Skalierung: Vtot, Atot skalieren mit case_y^2
        # weil Radius x case_y -> Flaeche x case_y^2.
        "vtot_exact": lambda cx, cy, cz: math.pi * ((cx / 2.0) * _por_scale(cy))**2 * cz,
        "atot_exact": lambda cx, cy, cz: math.pi * ((cx / 2.0) * _por_scale(cy))**2,
        "veff": _por_veff,
        "aeff": _por_aeff,
        "smax_ref": _por_smax,
        # R_s_base = Fixture-Basisradius. Effektiver R_s = R_s_base * case_y.
        "constants": {"R_s_base": 22.73, "nu": 0.22, "E": 607000},
    },
}


# =========================================================
# Oeffentliche Dispatcher-Funktionen
# =========================================================

def get_analytical_veff(case_type, m, vol_base,
                        dim_x=0.0, dim_y=0.0, dim_z=0.0,
                        crit="PIA"):
    """
    Gibt den analytischen Veff-Wert fuer gegebenen Lastfall und Weibull-Modul zurueck.

    Args:
        case_type: Geometrie-Typ ("BEAM1", "BEAM2", "PWH", "PWHR", "3PB", "POR")
        m: Weibull-Modul (int, 1..50)
        vol_base: Basisvolumen (analytisch oder numerisch)
        dim_x, dim_y, dim_z: Geometrie-Parameter
        crit: Kriterium ("PIA" oder "S1")

    Returns:
        float: Analytischer Veff-Wert, oder 0.0 bei Fehler
    """
    if vol_base == 0:
        return 0.0
    entry = LOADCASE_REGISTRY.get(case_type)
    if entry is None:
        print(f"   [analytical_helper] WARNUNG: Unbekannter Case '{case_type}'. "
              f"Nutze BEAM1-Formel als Fallback.")
        return vol_base / (2.0 * (float(m) + 1.0))
    return entry["veff"](m, vol_base, dim_x, dim_y, dim_z, crit)


def get_analytical_smax_ref(case_type, smax_nodal, load_n,
                            dim_x=0.0, dim_y=0.0, dim_z=0.0):
    """
    Bestimmt die analytische Referenz-Maximalspannung je Lastfall.

    Args:
        case_type: Geometrie-Typ ("BEAM1", "BEAM2", "PWH", "PWHR", "3PB", "POR")
        smax_nodal: Nodale Smax aus FEM (Fallback)
        load_n: Aufgebrachte Last (N oder MPa, lastfallabhaengig)
        dim_x, dim_y, dim_z: Geometrie-Parameter (case_x/y/z aus CSV)

    Returns:
        float: Analytische Referenz-Smax
    """
    entry = LOADCASE_REGISTRY.get(case_type)
    if entry is None:
        print(f"   [analytical_helper] WARNUNG: Unbekannter Case '{case_type}'. "
              f"Analytik nutzt Smax_nodal.")
        return smax_nodal
    return entry["smax_ref"](smax_nodal, load_n, dim_x, dim_y, dim_z)


def get_analytical_aeff(case_type, m, area_base,
                        dim_x=0.0, dim_y=0.0, dim_z=0.0,
                        crit="PIA"):
    """Analytische effektive Flaeche. Gibt None zurueck wenn keine Formel vorhanden.

    Args:
        case_type: Geometrie-Typ
        m: Weibull-Modul
        area_base: Gesamtflaeche (analytisch oder numerisch)
        dim_x, dim_y, dim_z: Geometrie-Parameter
        crit: Kriterium ("PIA" oder "S1")

    Returns:
        float oder None: Analytischer Aeff-Wert, oder None wenn keine Formel
    """
    entry = LOADCASE_REGISTRY.get(case_type)
    if entry is None:
        return None
    aeff_func = entry.get("aeff")
    if aeff_func is None:
        return None
    return aeff_func(m, area_base, dim_x, dim_y, dim_z, crit)
