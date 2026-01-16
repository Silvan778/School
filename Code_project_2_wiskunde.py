# modulen
import numpy as np
from scipy.signal import place_poles

# vaste waarden
ρ = 700 # kg/m3
d = 0.02 # m
l = 0.6 # m
g = 9.81 # m/s2

# berekende waarden 1
mp = ρ * d**2 * l # kg
Icm = (1/12) * mp * (d**2 + l**2)
I = Icm + mp * (l/2)**2

# berekende waarden 2
a = I + mp * l**2
b = (mp*l**2) / 2
c = Icm + (1/4) * mp * l**2

# berekende waarden van A en B matrix
A = np.array([[0, 0, 1, 0],
             [0, 0, 0, 1],
             [20.99, -20.97, 0, 0],
             [-28.03, 76.88, 0, 0]])

B = np.array([[0],
             [0],
             [28.32],
             [-70.89]])

# gain matrix bepalen met polen
Polen = np.array([-0.5,-10,-11,-100]) # zelf gekozen polen / eigenwaarde
Geplaatste_Polen = place_poles(A, B, Polen)

K = Geplaatste_Polen.gain_matrix
print(f"Gain matrix K (1x4):\n {K}")
