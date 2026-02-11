import wis_2_2_utilities_nochrono as util
import wis_2_2_systems_nochrono as systems
#import random
import numpy as np

#set timestep
timestep = 2e-5
second_pendulum = False

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
    self.matrix_gain=np.array([[0, 0, 0, 0]])
   
  def feedBack(self, observe):
    u= -self.matrix_gain @ observe
    return u  
 
class controller():
  def __init__(self, target=0):
    pass
   
  def feedBack(self, observe):
    u=0
    return u

def run(state0, u0, duration, eps, second_bool):
    model = systems.cart_inverted_pendulum(state0=state0,u0=u0, second_pendulum = second_bool)

    control = controller()
    simulation = util.simulation(model=model,timestep=timestep)
    simulation.setCost()
    simulation.max_duration = duration #seconde
    simulation.GIF_toggle = False #set to false to avoid frame and GIF creation

    while simulation.vis.Run():
        if simulation.time<simulation.max_duration:
          simulation.step()
          u = control.feedBack(simulation.observe())
          simulation.control(eps)
          simulation.log()
          simulation.refreshTime()
        else:
          print('Ending visualisation...')
          simulation.vis.GetDevice().closeDevice()
         
    simulation.writeData()
    data = simulation.observe()
    return data

def Main():
    size = 6 if second_pendulum else 4

    uitwijking = 0.01
    duratie = 0.01

    A = np.zeros((size, size))
    B = np.zeros((size, 1))

    # --- A matrix ---
    for i in range(size):
        state0 = np.zeros(size)
        state0[i] = uitwijking

        result = run(state0, 0, duratie, uitwijking, second_pendulum)

        # Finite difference linearisation
        derivative = (result - state0) / (uitwijking * duratie)

        A[:, i] = derivative

    # --- B matrix ---
    state0 = np.zeros(size)
    result = run(state0, uitwijking, duratie, uitwijking, second_pendulum)

    B[:, 0] = result / (uitwijking * duratie)

    # --- Print results ---
    print("A =")
    print(np.array2string(A, separator=", "))

    print("\nB =")
    print(np.array2string(B, separator=", "))

if __name__ == "__main__":
  Main()
