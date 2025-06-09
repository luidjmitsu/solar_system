import numpy as np 
import matplotlib.pyplot as plt
from planets import *

def energy(t, planet):
    plt.style.use('default')
    plt.figure()
    plt.plot(t, planet.Ec, label='Ec')
    plt.plot(t, planet.Ep, label='Ep')
    plt.plot(t, np.array(planet.Ec) + np.array(planet.Ep), label='Etot')
    plt.grid(True)
    plt.xlabel('t (s)')
    plt.ylabel('E (J)')
    plt.legend()
    plt.show()
    return None


