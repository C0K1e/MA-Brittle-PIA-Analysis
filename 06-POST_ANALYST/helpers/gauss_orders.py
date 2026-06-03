import numpy as np

# Generiert die exakten Gauss-Legendre Punkte und Gewichte für n=...
x, w = np.polynomial.legendre.leggauss(26)

print("  CASE (26)")
for i in range(26):
    print(f"    gausz_r({i+1}) = {x[i]: .25f}")
print("")
for i in range(26):
    print(f"    gausz_w({i+1}) = {w[i]: .25f}")