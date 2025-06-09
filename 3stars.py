import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import functions as f 
from planets import *
import energy as e

#------------------------------------------------------------
# Problème à 3 astres
#------------------------------------------------------------
Mj = 1.9e27

planets.append(Planets("Earth", Me, 1.017*AU, 0, 0, 0, e_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Jupiter", Mj, (5.2 + 0.048)*AU, 0, 0, 0, j_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Sun", Ms, 0, 0, 0, 0, 0, 0, [], [], [], [], []))

t = 0.0
dt = 1 * daysec
T = []

while t < 10 * yearsec:
    for p in planets: 
        p.vx, p.vy, p.vz, p.x, p.y, p.z = f.RK4(f.F_planet, t, p.vx, p.vy, p.vz, p.x, p.y, p.z, p.mass, p.name, planets, b, dt)
        p.x_list.append(p.x)
        p.y_list.append(p.y)
        p.z_list.append(p.z)

        modr = (p.x**2 + p.y**2 + p.z**2)**0.5
        modv = (p.vx**2 + p.vy**2 + p.vz**2)**0.5
        p.Ec.append(0.5*p.mass*modv**2)
        p.Ep.append(-G*Ms*p.mass/modr)

    T.append(t) 
    t += 4*dt 

print('Data ready')

#------------------------------------------------------------
# Plot
#------------------------------------------------------------
fig, ax = plt.subplots()
ax.set_xlim(-6*AU, 6*AU)
ax.set_ylim(-6*AU, 6*AU)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Solar system ; 3 stars (' + '\u03B2' + ' = ' + str(b) + ' ) ')
ax.grid()

# Plot the Earth and its orbit
earth, = ax.plot([], [], 'o', color='blue', markersize=5, label=planets[0].name)
earth_orbit, = ax.plot([], [], color='blue')

# Plot the Sun and its orbit
sun, = ax.plot([], [], 'o', color='yellow', markersize=20, label=planets[2].name)
sun_orbit, = ax.plot([], [], color='yellow')

# Plot Jupiter and its orbit
jupiter, = ax.plot([], [], 'o', color='orange', markersize=10, label='Jupiter')
jupiter_orbit, = ax.plot([], [], color='orange')

def init():
    earth.set_data([], [])
    earth_orbit.set_data([], [])
    sun.set_data([], [])
    sun_orbit.set_data([], [])
    jupiter.set_data([], [])
    jupiter_orbit.set_data([], [])
    return earth, earth_orbit, sun, sun_orbit, jupiter, jupiter_orbit

def animate(i):
    earth.set_data([planets[0].x_list[i]], [planets[0].y_list[i]])
    earth_orbit.set_data(planets[0].x_list[:i], planets[0].y_list[:i])
    sun.set_data([planets[2].x_list[i]], [planets[2].y_list[i]])
    sun_orbit.set_data(planets[2].x_list[:i], planets[2].y_list[:i])
    jupiter.set_data([planets[1].x_list[i]], [planets[1].y_list[i]])
    jupiter_orbit.set_data(planets[1].x_list[:i], planets[1].y_list[:i])
    return earth, earth_orbit, sun, sun_orbit, jupiter, jupiter_orbit

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
ax.set_title('Solar system ; 3 stars (' + '\u03B2' + ' = ' + str(b) + ' ) ')
ax.set_xlim(-6*AU, 6*AU)
ax.set_ylim(-6*AU, 6*AU)
ax.set_zlim(-6*AU, 6*AU)

# Plot the Earth and its orbit
earth, = ax.plot([], [], [], 'o', color='blue', markersize=5, label=planets[0].name)
earth_orbit, = ax.plot([], [], [], color='blue')

# Plot the Sun and its orbit
sun, = ax.plot([], [], [], 'o', color='yellow', markersize=20, label=planets[2].name)
sun_orbit, = ax.plot([], [], [], color='yellow')

# Plot Jupiter and its orbit
jupiter, = ax.plot([], [], [], 'o', color='orange', markersize=10, label=planets[1].name)
jupiter_orbit, = ax.plot([], [], [], color='orange')

def init():

    earth.set_data([], [])
    earth.set_3d_properties([])
    earth_orbit.set_data([], [])
    earth_orbit.set_3d_properties([])

    sun.set_data([], [])
    sun.set_3d_properties([])
    sun_orbit.set_data([], [])
    sun_orbit.set_3d_properties([])

    jupiter.set_data([], [])
    jupiter.set_3d_properties([])
    jupiter_orbit.set_data([], [])
    jupiter_orbit.set_3d_properties([])

    return earth, earth_orbit, sun, sun_orbit, jupiter, jupiter_orbit

def animate(i):

    earth.set_data([planets[0].x_list[i]], [planets[0].y_list[i]])
    earth.set_3d_properties([planets[0].z_list[i]])
    earth_orbit.set_data(planets[0].x_list[:i], planets[0].y_list[:i])
    earth_orbit.set_3d_properties(planets[0].z_list[:i])

    sun.set_data([planets[2].x_list[i]], [planets[2].y_list[i]])
    sun.set_3d_properties([planets[2].z_list[i]])
    sun_orbit.set_data(planets[2].x_list[:i], planets[2].y_list[:i])
    sun_orbit.set_3d_properties(planets[2].z_list[:i])

    jupiter.set_data([planets[1].x_list[i]], [planets[1].y_list[i]])
    jupiter.set_3d_properties([planets[1].z_list[i]])
    jupiter_orbit.set_data(planets[1].x_list[:i], planets[1].y_list[:i])
    jupiter_orbit.set_3d_properties(planets[1].z_list[:i])

    return earth, earth_orbit, sun, sun_orbit, jupiter, jupiter_orbit

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(planets[0].x_list), interval=1, blit=True, repeat=False)
plt.legend()
plt.show()
