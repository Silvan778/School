import wis_2_2_utilities as util
import wis_2_2_systems as systems

from scipy.linalg import solve_continuous_are
from scipy.signal import place_poles
from scipy import linalg

import pandas as pd
import matplotlib.pyplot as plt

#import random
import numpy as np

timestep = 2e-3

A = np.array([
    [ 0,  1, 0, 0],
    [ 0, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0,  24.63, 0],
])

B = np.array([
    [ 0],
    [ 0.09],
    [0],
    [-0.22],
])

macht = 6

Q = np.diag([10**macht, 10**-macht, 10**macht, 10**-macht])
R = np.array([[1]])

P = solve_continuous_are(A, B, Q, R)

K = np.linalg.inv(R) @ B.T @ P

#desired_poles = [-12,-13,-14,-51]
#result = place_poles(A, B, desired_poles)
#K = result.gain_matrix

class PID_controller():
  def __init__(self, target=0):
    self.integral1=0
    self.integral2=0
    self.K_P1=0
    self.K_I1=0
    self.K_D1=0
    self.K_P2=0
    self.K_I2=0
    self.K_D2=0
   
  def feedBack(self, observe):
    self.integral1+=observe[0]
    self.integral2+=observe[2]
    u=self.K_P1*observe[0]+\
      self.K_I1*self.integral1+\
      self.K_D1*observe[1]+\
      self.K_P2*observe[2]+\
      self.K_I2*self.integral2+\
      self.K_D2*observe[3]
    return u
 
class pp_controller():
  def __init__(self, target=0):
    self.matrix_gain=np.array(K)
   
  def feedBack(self, observe):
    u= -self.matrix_gain @ observe
    return u  
 
class controller():
  def __init__(self, target=0):
    pass
   
  def feedBack(self, observe):
    u=0
    return u


def main():
  model=systems.cart_inverted_pendulum()
  control = pp_controller()
  simulation = util.simulation(model=model,timestep=timestep)
  simulation.setCost()
  simulation.max_duration = 2 #seconde
  simulation.GIF_toggle = False #set to false to avoid frame and GIF creation

  while simulation.vis.Run():
      if simulation.time<simulation.max_duration:
        simulation.step()
        u = control.feedBack(simulation.observe())
        simulation.control(u)
        simulation.log()
        simulation.refreshTime()
      else:
        print('Ending visualisation...')
        simulation.vis.GetDevice().closeDevice()
       
  simulation.writeData()

if __name__ == "__main__":
  main()
  
plt.rcParams["font.family"] = "Calibri"

# =======================
# CSV-bestand inlezen
# =======================
bestand = "Cart_Inverted_Pendulum.csv"
df = pd.read_csv(bestand, header=None)

tijd = df.iloc[:, 0]

# Kolommen (0-based indexing)
kwad_toestanden = df.iloc[:, 1]
kwad_input = df.iloc[:, 2]
positie_kar = df.iloc[:, 3]
snelheid_kar = df.iloc[:, 4]
hoek_slinger_1 = df.iloc[:, 5]
hoeksnelheid_1 = df.iloc[:, 6]
input_signaal = df.iloc[:, 7]

# Lijst van grootheden om te plotten
plots = [
    ("Kwadratische toestandskosten", kwad_toestanden, "Kosten"),
    ("Kwadratische inputkosten", kwad_input, "Kosten"),
    ("Positie kar", positie_kar, "Positie (m)"),
    ("Snelheid kar", snelheid_kar, "Snelheid (m/s)"),
    ("Hoek slinger 1", hoek_slinger_1, "Hoek (rad)"),
    ("Hoeksnelheid slinger 1", hoeksnelheid_1, "Hoeksnelheid (rad/s)"),
    ("Input", input_signaal, "Input")
]

# =======================
# Plotten met maximale markers
# =======================
for titel, data, ylabel in plots:
    plt.figure(figsize=(8, 5))
    plt.plot(tijd, data, '--', label=titel)

    plt.xlabel("Tijd (s)", fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.title(titel, fontsize=25, pad=30, fontweight='bold')
    plt.grid(True)
    plt.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.25),
        ncol=1,
        fontsize=14
    )
    plt.tight_layout()
    plt.show()
