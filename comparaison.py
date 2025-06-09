import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import functions as f 
from planets import *
import energy as e


planets.append(Planets("Earth", Me, 1.017*AU, 0, 0, 0, e_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Sun", Ms, 0, 0, 0, 0, 0, 0, [], [], [], [], []))

t = 0.0
dt = 1 * daysec
T = []

result_rk4, result_euler = [[],[]], [[],[]]
for i in range(1, 150):
    result_rk4[0].append(i)
    result_rk4[1].append(f.efficiency(f.RK4, f.F_planet, t, planets[0].vx, planets[0].vy, planets[0].vz, planets[0].x, planets[0].y, planets[0].z, planets[0].mass, planets[0].name, planets, b, dt, i))

    result_euler[0].append(i)
    result_euler[1].append(f.efficiency(f.euler, f.F_planet, t, planets[0].vx, planets[0].vy, planets[0].vz, planets[0].x, planets[0].y, planets[0].z, planets[0].mass, planets[0].name, planets, b, dt, i))

plt.figure()
plt.xlabel('Step')
plt.ylabel('Speed (step/s)')
plt.title('Efficiency')
plt.plot(result_rk4[0], result_rk4[1], label='RK4')
plt.plot(result_euler[0], result_euler[1], label='Euler')
plt.legend()
plt.show()