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

A = np.array([[ 4.43635595e-04,  1.00244364e+00, -8.14827503e-04,  4.39449202e-04,
   6.01014193e-04,  4.44158813e-04],
 [ 8.89053930e-02,  8.89053930e-02, -1.62840055e-01,  8.76468528e-02,
   1.20445624e-01,  8.90628019e-02],
 [-9.50805927e-04, -9.50805927e-04,  1.07079566e-01,  1.00140852e+00,
  -1.06586769e-01, -1.30217349e-03],
 [-1.90620945e-01, -1.90620945e-01,  2.14244099e+01, -8.25885612e-02,
  -2.13352019e+01, -2.96277190e-01],
 [ 1.26790183e-03,  1.26790183e-03, -1.42791072e-01,  7.88814258e-04,
   3.87803643e-01,  1.00455361e+00],
 [ 2.54296989e-01,  2.54296989e-01, -2.85812275e+01,  1.10235328e-01,
   7.76162422e+01,  6.40893341e-01]
])

B = np.array([[ 0.00044541],
 [ 0.08908321],
 [-0.00095462],
 [-0.19100282],
 [ 0.00127299],
 [ 0.25480725]
])

macht = 6

Q = np.diag([10**macht, 10**-macht, 10**macht, 10**-macht, 10**macht, 10**-macht])
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
  model=systems.cart_inverted_pendulum(second_pendulum=True)
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

# =======================
# Maximale uitwijking & hoeken berekenen
# =======================
positie_kar = df.iloc[:, 3]           # kolom 4 = positie kar
hoek_slinger_1 = df.iloc[:, 5]        # kolom 6 = hoek slinger 1
hoek_slinger_2 = df.iloc[:, 7] if df.shape[1] > 7 else None

# Maxima
max_uitwijking_kar = positie_kar.abs().max()

if hoek_slinger_2 is not None:
    max_hoek_slingers = max(hoek_slinger_1.abs().max(), hoek_slinger_2.abs().max())
else:
    max_hoek_slingers = hoek_slinger_1.abs().max()

# Printen
print(f"• Maximale uitwijking car: {max_uitwijking_kar:.3f} m")
print(f"• Maximale hoek slingers: {max_hoek_slingers:.3f} rad")

# =======================
# Enkelvoudige plots
# =======================
single_plots = [
    ("Kwadratische toestandskosten", 1, "Kosten"),
    ("Kwadratische inputkosten", 2, "Kosten"),
    ("Positie kar", 3, "Positie (m)"),
    ("Snelheid kar", 4, "Snelheid (m/s)"),
    ("Input", 9, "Input")
]

for titel, kolom, ylabel in single_plots:
    if kolom < df.shape[1]:
        data = df.iloc[:, kolom]

        plt.figure(figsize=(8, 5))
        plt.plot(tijd, data, '--', label=titel)

        plt.xlabel("Tijd (s)", fontsize=16)
        plt.ylabel(ylabel, fontsize=16)
        plt.title(titel, fontsize=25, pad=30, fontweight='bold')
        plt.tick_params(axis='both', labelsize=12)
        plt.grid(True)

        plt.legend(
            loc='upper center',
            bbox_to_anchor=(0.5, -0.25),
            ncol=1,
            fontsize=14
        )

        plt.tight_layout()
        plt.show()

# =======================
# Hoeken samen
# =======================
if df.shape[1] > 7:
    plt.figure(figsize=(8, 5))
    plt.plot(tijd, df.iloc[:, 5], '--', label='Hoek slinger 1')
    plt.plot(tijd, df.iloc[:, 7], '--', label='Hoek slinger 2')

    plt.xlabel("Tijd (s)", fontsize=16)
    plt.ylabel("Hoek (rad)", fontsize=16)
    plt.title("Hoek slinger 1 en 2", fontsize=25, pad=30, fontweight='bold')
    plt.tick_params(axis='both', labelsize=12)
    plt.grid(True)

    plt.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.25),
        ncol=2,
        fontsize=14
    )

    plt.tight_layout()
    plt.show()

# =======================
# Hoeksnelheden samen
# =======================
if df.shape[1] > 8:
    plt.figure(figsize=(8, 5))
    plt.plot(tijd, df.iloc[:, 6], '--', label='Hoeksnelheid slinger 1')
    plt.plot(tijd, df.iloc[:, 8], '--', label='Hoeksnelheid slinger 2')

    plt.xlabel("Tijd (s)", fontsize=16)
    plt.ylabel("Hoeksnelheid (rad/s)", fontsize=16)
    plt.title("Hoeksnelheid slinger 1 en 2", fontsize=25, pad=30, fontweight='bold')
    plt.tick_params(axis='both', labelsize=12)
    plt.grid(True)

    plt.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.25),
        ncol=2,
        fontsize=14
    )

    plt.tight_layout()
    plt.show()
