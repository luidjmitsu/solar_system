import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import functions as f 
from planets import *
import energy as e

#------------------------------------------------------------
# Changement de la constante beta
#------------------------------------------------------------
b = 3.0

#------------------------------------------------------------
# Etude de Mercure
#------------------------------------------------------------
planets.append(Planets("Mercury", Mm, 0.387*AU, 0, 0, 0, m_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Sun", Ms, 0, 0, 0, 0, 0, 0, [], [], [], [], []))

t = 0.0
dt = 1 * daysec
T = []

while t < yearsec:
    for p in planets: 
        p.vx, p.vy, p.vz, p.x, p.y, p.z = f.RK4(f.F_planet, t, p.vx, p.vy, p.vz, p.x, p.y, p.z, p.mass, p.name, planets, b, dt)
        p.x_list.append(p.x)
        p.y_list.append(p.y)
        p.z_list.append(p.z)

        modr = (p.x**2 + p.y**2 + p.z**2)**0.5
        modv = (p.vx**2 + p.vy**2 + p.vz**2)**0.5
        p.Ec.append(0.5*p.mass*modv**2)  
        p.Ep.append(-G*Ms*p.mass/((b-1)*modr**(b-1)) )

    T.append(t)    
    t += 2*dt

print('Data ready')

#------------------------------------------------------------
# Plot
#------------------------------------------------------------
fig, ax = plt.subplots()
ax.set_xlim(-1*AU, 1*AU)
ax.set_ylim(-1*AU, 1*AU)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title("Mercury's orbit (" + '\u03B2' + ' = ' + str(b) + ' ) ' )
ax.grid()

# Plot Mercury and its orbit
mercury, = ax.plot([], [], 'o', color='red', markersize=5, label=planets[0].name)
mercury_orbit, = ax.plot([], [], color='red')

# Plot the Sun and its orbit
sun, = ax.plot([], [], 'o', color='yellow', markersize=20, label=planets[1].name)
sun_orbit, = ax.plot([], [], color='yellow')

def init():
    mercury.set_data([], [])
    mercury_orbit.set_data([], [])
    sun.set_data([], [])
    sun_orbit.set_data([], [])
    return mercury, mercury_orbit, sun, sun_orbit

def animate(i):
    mercury.set_data([planets[0].x_list[i]], [planets[0].y_list[i]])
    mercury_orbit.set_data(planets[0].x_list[:i], planets[0].y_list[:i])
    sun.set_data([planets[1].x_list[i]], [planets[1].y_list[i]])
    sun_orbit.set_data(planets[1].x_list[:i], planets[1].y_list[:i])
    return mercury, mercury_orbit, sun, sun_orbit

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(planets[0].x_list), interval=1, blit=True, repeat=False)
plt.legend()
plt.show()

#------------------------------------------------------------
# 3D plot
#------------------------------------------------------------
plt.style.use('dark_background')

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.grid(False)
ax.axis('off')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title("mercury's orbit (" + '\u03B2' + ' = ' + str(b) + ' ) ' )
ax.grid()
ax.set_xlim(-1*AU, 1*AU)
ax.set_ylim(-1*AU, 1*AU)
ax.set_zlim(-1*AU, 1*AU)

# Plot Mercury and its orbit
mercury, = ax.plot([], [], [], 'o', color='red', markersize=5, label=planets[0].name)
mercury_orbit, = ax.plot([], [], [], color='red')

# Plot the Sun and its orbit
sun, = ax.plot([], [], [], 'o', color='yellow', markersize=20, label=planets[1].name)
sun_orbit, = ax.plot([], [], [], color='yellow')

def init():

    mercury.set_data([], [])
    mercury.set_3d_properties([])
    mercury_orbit.set_data([], [])
    mercury_orbit.set_3d_properties([])

    sun.set_data([], [])
    sun.set_3d_properties([])
    sun_orbit.set_data([], [])
    sun_orbit.set_3d_properties([])

    return mercury, mercury_orbit, sun, sun_orbit

def animate(i):

    mercury.set_data([planets[0].x_list[i]], [planets[0].y_list[i]])
    mercury.set_3d_properties([planets[0].z_list[i]])
    mercury_orbit.set_data(planets[0].x_list[:i], planets[0].y_list[:i])
    mercury_orbit.set_3d_properties(planets[0].z_list[:i])

    sun.set_data([planets[1].x_list[i]], [planets[1].y_list[i]])
    sun.set_3d_properties([planets[1].z_list[i]])
    sun_orbit.set_data(planets[1].x_list[:i], planets[1].y_list[:i])
    sun_orbit.set_3d_properties(planets[1].z_list[:i])

    return mercury, mercury_orbit, sun, sun_orbit

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(planets[0].x_list), interval=1, blit=True, repeat=False)
plt.legend()
plt.show()
