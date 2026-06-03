# pf_helper.py
"""v20.0 (Bug #4): Numerisch robuste Pf-Berechnung im Log-Raum.

Ersetzt die fragile Inline-Formel `1 - exp(-Veff * (smax/sigma0)^m)`:
  - Auslöschung bei sehr kleinen Hazards (1 - 0.999... verliert Stellen)
  - Overflow bei sehr grossen Hazards (Power-Term ueberlaeuft double precision)

Mathematik:
    Pf = 1 - exp(-H)    mit Hazard H = Veff * (smax/sigma0)^m

Numerische Strategie:
  1. Im Log-Raum rechnen: log_h = log(Veff) + m * log(smax/sigma0)
  2. Bei log_h > log(36) -> Pf = 1.0 (Pf in double precision exakt 1)
  3. Bei log_h < -745 -> Pf ~ exp(log_h) (exp(log_h) selbst nahe Underflow)
  4. Sonst: -expm1(-h) statt 1 - exp(-h) (numerisch stabil bei kleinem h)
"""
import math


def pf_from_hazard(v_eff, smax, sigma0, m):
    """Robuste Pf-Berechnung im Log-Raum mit Clamping fuer Extremwerte.

    Args:
        v_eff: Effektives Volumen (oder Aeff bei LOC=A)
        smax:  Maximalspannung (Referenz-Norm)
        sigma0: Sigma0 (Weibull-Charakteristik)
        m: Weibull-Modul (int oder float)

    Returns:
        Pf in [0, 1]. Numerisch korrekt auch bei sehr kleinen oder grossen Hazards.
    """
    if v_eff <= 0 or smax <= 0 or sigma0 <= 0:
        return 0.0

    log_h = math.log(v_eff) + float(m) * math.log(smax / sigma0)

    # Hazard so gross, dass Pf in double precision exakt 1.0 ist
    # (1 - exp(-36) ~ 1 - 2.3e-16 ~ 1.0 in double precision)
    if log_h > math.log(36.0):
        return 1.0

    # Hazard so klein, dass exp(log_h) selbst Underflow erleidet (< 5e-324)
    # Pf ~ Hazard fuer sehr kleine Hazards (Taylor: 1-exp(-h) ~ h - h^2/2 + ...)
    if log_h < -745.0:
        return math.exp(log_h)

    h = math.exp(log_h)
    return -math.expm1(-h)   # = 1 - exp(-h), numerisch genau bei kleinen h
