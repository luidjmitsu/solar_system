import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import functions as f 
from planets import *

#------------------------------------------------------------
# Toutes les planètes du système solaire
#------------------------------------------------------------
planets.append(Planets("Mercury", Mm, 0.387*AU, 0, 0, 0, m_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Venus", Mv, 0.723*AU, 0, 0, 0, v_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Earth", Me, 1.017*AU, 0, 0, 0, e_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Mars", MM, 1.666*AU, 0, 0, 0, M_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Jupiter", Mj, (5.2 + 0.048)*AU, 0, 0, 0, j_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Saturn", MS, (9.54 + 0.056)*AU, 0, 0, 0, S_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Uranus", Mu, (19.19 + 0.046)*AU, 0, 0, 0, u_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Neptune", Mn, (30.06 + 0.010)*AU, 0, 0, 0, n_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Pluto", Mp, (39.53 + 0.248)*AU, 0, 0, 0, p_ap_v, 0, [], [], [], [], []))
planets.append(Planets("Sun", Ms, 0, 0, 0, 0, 0, 0, [], [], [], [], []))

t=0.0
dt = 1 * daysec
T = []

while t < 100 * yearsec:
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
    t += 4*dt

print('Data ready')

#------------------------------------------------------------
# Plot
#------------------------------------------------------------
fig, ax = plt.subplots()
ax.set_xlim(-50*AU, 50*AU)
ax.set_ylim(-50*AU, 50*AU)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Solar system (' + '\u03B2' + ' = ' + str(b) + ' ) ')
ax.grid()

# Plot the Earth and its orbit
earth, = ax.plot([], [], 'o', color='blue', markersize=5, label=planets[2].name)
earth_orbit, = ax.plot([], [], color='blue')

# Plot the Sun and its orbit
sun, = ax.plot([], [], 'o', color='yellow', markersize=20, label=planets[9].name)
sun_orbit, = ax.plot([], [], color='yellow')

# Plot the other planets and their orbits
mercury, = ax.plot([], [], 'o', color='red', markersize=5, label=planets[0].name)
mercury_orbit, = ax.plot([], [], color='red')
venus, = ax.plot([], [], 'o', color='orange', markersize=5, label=planets[1].name)
venus_orbit, = ax.plot([], [], color='orange')
mars, = ax.plot([], [], 'o', color='brown', markersize=5, label=planets[3].name)
mars_orbit, = ax.plot([], [], color='brown')
jupiter, = ax.plot([], [], 'o', color='purple', markersize=10, label=planets[4].name)
jupiter_orbit, = ax.plot([], [], color='purple')
saturn, = ax.plot([], [], 'o', color='green', markersize=10, label=planets[5].name)
saturn_orbit, = ax.plot([], [], color='green')
uranus, = ax.plot([], [], 'o', color='cyan', markersize=10, label=planets[6].name)
uranus_orbit, = ax.plot([], [], color='cyan')
neptune, = ax.plot([], [], 'o', color='blue', markersize=10, label=planets[7].name)
neptune_orbit, = ax.plot([], [], color='blue')
pluto, = ax.plot([], [], 'o', color='pink', markersize=5, label=planets[8].name)
pluto_orbit, = ax.plot([], [], color='pink')

def init():
    earth.set_data([], [])
    earth_orbit.set_data([], [])
    sun.set_data([], [])
    sun_orbit.set_data([], [])
    mercury.set_data([], [])
    mercury_orbit.set_data([], [])
    venus.set_data([], [])
    venus_orbit.set_data([], [])
    mars.set_data([], [])
    mars_orbit.set_data([], [])
    jupiter.set_data([], [])
    jupiter_orbit.set_data([], [])
    saturn.set_data([], [])
    saturn_orbit.set_data([], [])
    uranus.set_data([], [])
    uranus_orbit.set_data([], [])
    neptune.set_data([], [])
    neptune_orbit.set_data([], [])
    pluto.set_data([], [])
    pluto_orbit.set_data([], [])
    return earth, earth_orbit, sun, sun_orbit, mercury, mercury_orbit, venus, venus_orbit, mars, mars_orbit, jupiter, jupiter_orbit, saturn, saturn_orbit, uranus, uranus_orbit, neptune, neptune_orbit, pluto, pluto_orbit

def animate(i):
    earth.set_data([planets[2].x_list[i]], [planets[2].y_list[i]])
    earth_orbit.set_data(planets[2].x_list[:i], planets[2].y_list[:i])
    sun.set_data([planets[9].x_list[i]], [planets[9].y_list[i]])
    sun_orbit.set_data(planets[9].x_list[:i], planets[9].y_list[:i])
    mercury.set_data([planets[0].x_list[i]], [planets[0].y_list[i]])
    mercury_orbit.set_data(planets[0].x_list[:i], planets[0].y_list[:i])
    venus.set_data([planets[1].x_list[i]], [planets[1].y_list[i]])
    venus_orbit.set_data(planets[1].x_list[:i], planets[1].y_list[:i])
    mars.set_data([planets[3].x_list[i]], [planets[3].y_list[i]])
    mars_orbit.set_data(planets[3].x_list[:i], planets[3].y_list[:i])
    jupiter.set_data([planets[4].x_list[i]], [planets[4].y_list[i]])
    jupiter_orbit.set_data(planets[4].x_list[:i], planets[4].y_list[:i])
    saturn.set_data([planets[5].x_list[i]], [planets[5].y_list[i]])
    saturn_orbit.set_data(planets[5].x_list[:i], planets[5].y_list[:i])
    uranus.set_data([planets[6].x_list[i]], [planets[6].y_list[i]])
    uranus_orbit.set_data(planets[6].x_list[:i], planets[6].y_list[:i])
    neptune.set_data([planets[7].x_list[i]], [planets[7].y_list[i]])
    neptune_orbit.set_data(planets[7].x_list[:i], planets[7].y_list[:i])
    pluto.set_data([planets[8].x_list[i]], [planets[8].y_list[i]])
    pluto_orbit.set_data(planets[8].x_list[:i], planets[8].y_list[:i])
    return earth, earth_orbit, sun, sun_orbit, mercury, mercury_orbit, venus, venus_orbit, mars, mars_orbit, jupiter, jupiter_orbit, saturn, saturn_orbit, uranus, uranus_orbit, neptune, neptune_orbit, pluto, pluto_orbit

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(planets[2].x_list), interval=1, blit=True, repeat=False)
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
ax.set_title('Solar system (' + '\u03B2' + ' = ' + str(b) + ' ) ')
ax.set_xlim(-50*AU, 50*AU)
ax.set_ylim(-50*AU, 50*AU)
ax.set_zlim(-50*AU, 50*AU)

# Plot the Earth and its orbit
earth, = ax.plot([], [], [], 'o', color='blue', markersize=5, label=planets[2].name)
earth_orbit, = ax.plot([], [], [], color='blue')

# Plot the Sun and its orbit
sun, = ax.plot([], [], [], 'o', color='yellow', markersize=20, label=planets[9].name)
sun_orbit, = ax.plot([], [], [], color='yellow')

# Plot the other planets and their orbits
mercury, = ax.plot([], [], [], 'o', color='red', markersize=5, label=planets[0].name)
mercury_orbit, = ax.plot([], [], [], color='red')
venus, = ax.plot([], [], [], 'o', color='orange', markersize=5, label=planets[1].name)
venus_orbit, = ax.plot([], [], [], color='orange')
mars, = ax.plot([], [], [], 'o', color='brown', markersize=5, label=planets[3].name)
mars_orbit, = ax.plot([], [], [], color='brown')
jupiter, = ax.plot([], [], [], 'o', color='purple', markersize=10, label=planets[4].name)
jupiter_orbit, = ax.plot([], [], [], color='purple')
saturn, = ax.plot([], [], [], 'o', color='green', markersize=10, label=planets[5].name)
saturn_orbit, = ax.plot([], [], [], color='green')
uranus, = ax.plot([], [], [], 'o', color='cyan', markersize=10, label=planets[6].name)
uranus_orbit, = ax.plot([], [], [], color='cyan')
neptune, = ax.plot([], [], [], 'o', color='blue', markersize=10, label=planets[7].name)
neptune_orbit, = ax.plot([], [], [], color='blue')
pluto, = ax.plot([], [], [], 'o', color='pink', markersize=5, label=planets[8].name)
pluto_orbit, = ax.plot([], [], [], color='pink')

def init():
    earth.set_data([], [])
    earth_orbit.set_data([], [])
    sun.set_data([], [])
    sun_orbit.set_data([], [])
    mercury.set_data([], [])
    mercury_orbit.set_data([], [])
    venus.set_data([], [])
    venus_orbit.set_data([], [])
    mars.set_data([], [])
    mars_orbit.set_data([], [])
    jupiter.set_data([], [])
    jupiter_orbit.set_data([], [])
    saturn.set_data([], [])
    saturn_orbit.set_data([], [])
    uranus.set_data([], [])
    uranus_orbit.set_data([], [])
    neptune.set_data([], [])
    neptune_orbit.set_data([], [])
    pluto.set_data([], [])
    pluto_orbit.set_data([], [])
    return earth, earth_orbit, sun, sun_orbit, mercury, mercury_orbit, venus, venus_orbit, mars, mars_orbit, jupiter, jupiter_orbit, saturn, saturn_orbit, uranus, uranus_orbit, neptune, neptune_orbit, pluto, pluto_orbit

def animate(i):
    earth.set_data([planets[2].x_list[i]], [planets[2].y_list[i]])
    earth.set_3d_properties([planets[2].z_list[i]])
    earth_orbit.set_data(planets[2].x_list[:i], planets[2].y_list[:i])
    earth_orbit.set_3d_properties(planets[2].z_list[:i])
    sun.set_data([planets[9].x_list[i]], [planets[9].y_list[i]])
    sun.set_3d_properties([planets[9].z_list[i]])
    sun_orbit.set_data(planets[9].x_list[:i], planets[9].y_list[:i])
    sun_orbit.set_3d_properties(planets[9].z_list[:i])
    mercury.set_data([planets[0].x_list[i]], [planets[0].y_list[i]])
    mercury.set_3d_properties([planets[0].z_list[i]])
    mercury_orbit.set_data(planets[0].x_list[:i], planets[0].y_list[:i])
    mercury_orbit.set_3d_properties(planets[0].z_list[:i])
    venus.set_data([planets[1].x_list[i]], [planets[1].y_list[i]])
    venus.set_3d_properties([planets[1].z_list[i]])
    venus_orbit.set_data(planets[1].x_list[:i], planets[1].y_list[:i])
    venus_orbit.set_3d_properties(planets[1].z_list[:i])
    mars.set_data([planets[3].x_list[i]], [planets[3].y_list[i]])
    mars.set_3d_properties([planets[3].z_list[i]])
    mars_orbit.set_data(planets[3].x_list[:i], planets[3].y_list[:i])
    mars_orbit.set_3d_properties(planets[3].z_list[:i])
    jupiter.set_data([planets[4].x_list[i]], [planets[4].y_list[i]])
    jupiter.set_3d_properties([planets[4].z_list[i]])
    jupiter_orbit.set_data(planets[4].x_list[:i], planets[4].y_list[:i])
    jupiter_orbit.set_3d_properties(planets[4].z_list[:i])
    saturn.set_data([planets[5].x_list[i]], [planets[5].y_list[i]])
    saturn.set_3d_properties([planets[5].z_list[i]])
    saturn_orbit.set_data(planets[5].x_list[:i], planets[5].y_list[:i])
    saturn_orbit.set_3d_properties(planets[5].z_list[:i])
    uranus.set_data([planets[6].x_list[i]], [planets[6].y_list[i]])
    uranus.set_3d_properties([planets[6].z_list[i]])
    uranus_orbit.set_data(planets[6].x_list[:i], planets[6].y_list[:i])
    uranus_orbit.set_3d_properties(planets[6].z_list[:i])
    neptune.set_data([planets[7].x_list[i]], [planets[7].y_list[i]])
    neptune.set_3d_properties([planets[7].z_list[i]])
    neptune_orbit.set_data(planets[7].x_list[:i], planets[7].y_list[:i])
    neptune_orbit.set_3d_properties(planets[7].z_list[:i])
    pluto.set_data([planets[8].x_list[i]], [planets[8].y_list[i]])
    pluto.set_3d_properties([planets[8].z_list[i]])
    pluto_orbit.set_data(planets[8].x_list[:i], planets[8].y_list[:i])
    pluto_orbit.set_3d_properties(planets[8].z_list[:i])
    return earth, earth_orbit, sun, sun_orbit, mercury, mercury_orbit, venus, venus_orbit, mars, mars_orbit, jupiter, jupiter_orbit, saturn, saturn_orbit, uranus, uranus_orbit, neptune, neptune_orbit, pluto, pluto_orbit

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(planets[2].x_list), interval=1, blit=True, repeat=False)
plt.legend()
plt.show()

