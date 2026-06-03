# analytical_stress_fields.py
# Analytische Spannungsfelder an beliebigen Knotenkoordinaten fuer VTK-Export.
#
# Stellt compute_analytical_stresses() bereit, das fuer jeden Lastfall
# die analytischen Hauptspannungen (S1, S2, S3) an gegebenen Punkten berechnet.
#
# Formeln aus:
# - 08-PLOTS/spezial/08-PLOT-Pfad-Spannungen.py (Kirsch-Loesung, PWH)
# - 06-POST_ANALYST/analytical/analytical_helper.py (POR-Spannungsfelder)
# - Analytische Balkentheorie (BEAM1, BEAM2, 3PB)

import math


EPS = 1e-10  # Singularitaets-Schutz (statt harter Maske)


# =========================================================
# PWH: Kirsch-Loesung (2D, unendliche Platte mit Loch)
# =========================================================

def _kirsch_stresses(x, y, a, sigma_inf):
    """Erste Hauptspannung nach Kirsch fuer Platte mit Loch.

    Aus 08-PLOT-Pfad-Spannungen.py:99-116 (Mathematica-Formel).
    Singularitaets-Schutz via EPS statt harter Maske.

    Args:
        x, y: Koordinaten [mm]
        a: Lochradius [mm]
        sigma_inf: Fernfeld-Zugspannung [MPa]

    Returns:
        (s1, s2): Erste und zweite Hauptspannung [MPa]
    """
    r2 = max(x**2 + y**2, EPS**2)

    term1 = (2 * a**2 * (-x**2 + y**2)) / r2**2

    under_sqrt = (9 * a**8
                  + r2**4
                  - 6 * a**6 * (3 * x**2 + y**2)
                  - 2 * a**2 * r2 * (3 * x**4 - 12 * x**2 * y**2 + y**4)
                  + a**4 * (15 * x**4 - 26 * x**2 * y**2 + 7 * y**4))
    under_sqrt = max(under_sqrt, 0.0)

    term2 = math.sqrt(under_sqrt) / r2**2

    s1 = 0.5 * sigma_inf * (1 + term1 + term2)
    s2 = 0.5 * sigma_inf * (1 + term1 - term2)
    return s1, s2


# =========================================================
# POR: Ringdruckversuch (Spannungsfeld mit Dickenabhaengigkeit)
# =========================================================
# Importiert Formeln aus analytical_helper.py

def _por_sigma_rr(r, z, p, R_s, R_d, t, nu):
    """Radialspannung sigma_rr(r, z) mit linearer Dickenabhaengigkeit.

    z=0: Zugseite (max. Biegung positiv), z=t: Druckseite (Biegung negativ).
    """
    bracket = ((1 - nu) * R_s**2 / R_d**2
               + 2 * (1 + nu)
               - (3 + nu) * r**2 / R_s**2)
    sig_bend = (3 * p * R_s**2) / (8 * t**2) * bracket
    sig_memb = p * (3 + nu) / (4 * (1 - nu))
    return sig_bend * (1.0 - 2.0 * z / t) + sig_memb


def _por_sigma_phiphi(r, z, p, R_s, R_d, t, nu):
    """Tangentialspannung sigma_phiphi(r, z) mit linearer Dickenabhaengigkeit.

    z=0: Zugseite (max. Biegung positiv), z=t: Druckseite (Biegung negativ).
    """
    bracket = ((1 - nu) * R_s**2 / R_d**2
               + 2 * (1 + nu)
               - (1 + 3 * nu) * r**2 / R_s**2)
    sig_bend = (3 * p * R_s**2) / (8 * t**2) * bracket
    sig_memb = p * (3 + nu) / (4 * (1 - nu))
    return sig_bend * (1.0 - 2.0 * z / t) + sig_memb


# =========================================================
# Dispatcher
# =========================================================

def compute_analytical_stresses(case_type, points, case_x, case_y, case_z, load_n):
    """Analytische Hauptspannungen an Knotenpositionen.

    Args:
        case_type: Lastfall-Typ ("PWH", "POR", "BEAM1", "BEAM2", "3PB", ...)
        points: list of (x, y, z) tuples [mm]
        case_x, case_y, case_z: Geometrie-Parameter aus CSV
        load_n: Last (N oder MPa, lastfallabhaengig)

    Returns:
        (s1_list, s2_list, s3_list) — Listen mit je len(points) Eintraegen.
        None wenn keine analytische Formel vorhanden.
    """
    n = len(points)

    if case_type == "PWH":
        a = case_z             # Lochradius = CASE_Z
        sigma_inf = load_n     # Fernfeld-Zug = LOAD_N [MPa]
        s1_list, s2_list, s3_list = [], [], []
        for x, y, z in points:
            s1, s2 = _kirsch_stresses(x, y, a, sigma_inf)
            s1_list.append(s1)
            s2_list.append(s2)
            s3_list.append(0.0)  # Ebener Spannungszustand
        return s1_list, s2_list, s3_list

    elif case_type == "POR":
        R_d = case_x / 2.0    # Aussendurchmesser / 2
        t = case_z             # Plattendicke
        R_s = 22.73            # Auflagerring-Radius (hardcoded)
        nu = 0.22              # Poisson-Zahl (wie geom_POR.inp)
        p = load_n             # Druck [MPa]
        s1_list, s2_list, s3_list = [], [], []
        for x, y, z in points:
            r = max(math.sqrt(x**2 + y**2), EPS)
            sig_rr = _por_sigma_rr(r, z, p, R_s, R_d, t, nu)
            sig_pp = _por_sigma_phiphi(r, z, p, R_s, R_d, t, nu)
            # Hauptspannungen sortieren (s1 >= s2 >= s3).
            # sigma_zz = 0 (freie Ober-/Unterflaeche). Muss als dritter Kandidat
            # mitsortiert werden: Auf der Druckseite (z=t) sind sig_rr und sig_pp
            # negativ -> S1=0, S3=min(sig_rr, sig_pp) ~ -sigma_max.
            all_s = sorted([sig_rr, sig_pp, 0.0], reverse=True)
            s1_list.append(all_s[0])
            s2_list.append(all_s[1])
            s3_list.append(all_s[2])
        return s1_list, s2_list, s3_list

    elif case_type == "BEAM1":
        # Reines Biegemoment: sigma_x = sigma_max * (2y/H - 1)
        # sigma_max = 6000 * load_F / (T * H^2), mit M = 1000*F
        if case_y <= 0 or case_z <= 0:
            return None
        sigma_max = 6000.0 * load_n / (case_z * case_y**2)
        s1_list, s2_list, s3_list = [], [], []
        for x, y, z in points:
            sigma_x = sigma_max * (2.0 * y / case_y - 1.0)
            s1_list.append(max(sigma_x, 0.0))
            s2_list.append(0.0)
            s3_list.append(min(sigma_x, 0.0))
        return s1_list, s2_list, s3_list

    elif case_type == "BEAM2":
        # Kragarm: sigma_x = 1500 * (x/L)^2 * (2y/H - 1)
        if case_x <= 0 or case_y <= 0:
            return None
        s1_list, s2_list, s3_list = [], [], []
        for x, y, z in points:
            sigma_x = 1500.0 * (x / case_x)**2 * (2.0 * y / case_y - 1.0)
            s1_list.append(max(sigma_x, 0.0))
            s2_list.append(0.0)
            s3_list.append(min(sigma_x, 0.0))
        return s1_list, s2_list, s3_list

    elif case_type == "3PB":
        # 3-Punkt-Biegung: sigma_x = sigma_max * (1 - 2|x|/L) * (1 - 2y/H)
        # fuer |x| < L/2 (Stuetzweite), sonst 0
        if case_x <= 0 or case_y <= 0 or case_z <= 0:
            return None
        sigma_max = 3.0 * load_n * case_x / (2.0 * case_z * case_y**2)
        half_span = case_x / 2.0
        s1_list, s2_list, s3_list = [], [], []
        for x, y, z in points:
            if abs(x) <= half_span:
                sigma_x = sigma_max * (1.0 - 2.0 * abs(x) / case_x) * (1.0 - 2.0 * y / case_y)
            else:
                sigma_x = 0.0
            s1_list.append(max(sigma_x, 0.0))
            s2_list.append(0.0)
            s3_list.append(min(sigma_x, 0.0))
        return s1_list, s2_list, s3_list

    else:
        # PWHR, unbekannte Cases: Keine geschlossene Formel
        return None
