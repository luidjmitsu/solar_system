import time
from planets import *


def F_planet(t, vx, vy, vz, x, y, z, M, name, planets, beta):
    fx, fy, fz = 0, 0, 0

    for p in planets:
        if name==p.name:
            continue

        rx, ry, rz = x - p.x, y - p.y, z - p.z
        modr_beta = (rx**2 + ry**2 + rz**2)**((beta+1)/2)
        K = G * M * p.mass

        fx += -K * rx / modr_beta
        fy += -K * ry / modr_beta
        fz += -K * rz / modr_beta
        
    dvx, dvy, dvz = fx/M, fy/M, fz/M
    dx, dy, dz = vx, vy, vz
    return dvx, dvy, dvz, dx, dy, dz


def RK4(system, t, vx, vy, vz, x, y, z, M, name, planets, beta, dt):
    k1, l1, m1, n1, o1, p1 = system(t, vx, vy, vz, x, y, z, M, name, planets, beta)
    k2, l2, m2, n2, o2, p2 = system(t + dt/2, vx + dt*k1/2, vy + dt*l1/2, vz + dt*m1/2, x + dt*n1/2, y + dt*o1/2, z + dt*p1/2, M, name, planets, beta)
    k3, l3, m3, n3, o3, p3 = system(t + dt/2, vx + dt*k2/2, vy + dt*l2/2, vz + dt*m2/2, x + dt*n2/2, y + dt*o2/2, z + dt*p2/2, M, name, planets, beta)
    k4, l4, m4, n4, o4, p4 = system(t + dt, vx + dt*k3, vy + dt*l3, vz + dt*m3, x + dt*n3, y + dt*o3, z + dt*p3, M, name, planets, beta)
    return vx + dt*(k1 + 2*k2 + 2*k3 + k4)/6, vy + dt*(l1 + 2*l2 + 2*l3 + l4)/6, vz + dt*(m1 + 2*m2 + 2*m3 + m4)/6, x + dt*(n1 + 2*n2 + 2*n3 + n4)/6, y + dt*(o1 + 2*o2 + 2*o3 + o4)/6, z + dt*(p1 + 2*p2 + 2*p3 + p4)/6


def euler(system, t, vx, vy, vz, x, y, z, M, name, planets, beta, dt):
    dvx, dvy, dvz, dx, dy, dz = system(t, vx, vy, vz, x, y, z, M, name, planets, beta)
    return vx + dt*dvx, vy + dt*dvy, vz + dt*dvz, x + dt*dx, y + dt*dy, z + dt*dz


def efficiency(method, system, t, vx, vy, vz, x, y, z, M, name, planets, beta, dt, N):
    result = [[],[]]
    t0 = time.time()
    for i in range(N):
        method(system, t, vx, vy, vz, x, y, z, M, name, planets, beta, i*dt)
    t1 = time.time()
    return (t1-t0)/N