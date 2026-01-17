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
a = (3/2) * mp * g * l
b = I + mp * l**2
c = (mp*l**2) / 2
d = Icm + (1/4) * mp * l**2

# Matrix M1, M2, M3 berekenen
M1 = np.array([[1, 0, 0, 0],
              [0,b+c, 0, c],
              [0,0,1,0],
              [0,c+d,0,d]
])

M2 = np.array([[0, 1, 0, 0],
              [-a, 0, 0, 0],
              [0,0,0,1],
              [-c,0,-c,0]
])

M3 = np.array([[0],
              [1],
              [0],
              [0]
])

# Matrix A en B berekenen
A = np.linalg.inv(M1) @ M2
B = np.linalg.inv(M1) @ M3

# gain matrix bepalen met polen
Polen = np.array([-0.5,-10,-11,-100]) # zelf gekozen polen / eigenwaarde
Geplaatste_Polen = place_poles(A, B, Polen)

K = Geplaatste_Polen.gain_matrix
print(f"Gain matrix K (1x4):\n {K}")
