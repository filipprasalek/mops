import numpy as np
import matplotlib.pyplot as plt
from math import pow

length = 100
width = 5
depth = 1
injection_point = 10
measurment_point = 90
tracer_amount = 1  # kg

U = 0.1  # m/s
D = 0.01  # m^2/delta

time = 1000  # s
dt = 1
dx = 1
nt = time // dt
nx = length // dx


def simulation(U, D):
    Ca = U * dt / dx
    Cd = D * dt / pow(dx, 2)

    c = np.zeros((nt, nx))

    for x in range(nx):
        if x == injection_point:
            c[0, x] = tracer_amount / (dx * width * depth)
        else:
            c[0, x] = 0

    for t in range(nt - 1):
        for x in range(2, nx - 1):
            component1 = (Cd * (1 - Ca) - Ca / 6 * (pow(Ca, 2) - 3 * Ca + 2)) * c[t, x + 1]
            component2 = (Cd * (2 - 3 * Ca) - Ca / 2 * (pow(Ca, 2) - 2 * Ca - 1)) * c[t, x]
            component3 = (Cd * (1 - 3 * Ca) - Ca / 2 * (pow(Ca, 2) - Ca - 2)) * c[t, x - 1]
            component4 = (Cd * Ca + Ca / 6 * (pow(Ca, 2) - 1)) * c[t, x - 2]
            c[t + 1, x] = c[t, x] + component1 - component2 + component3 + component4
    return c

def simulation_mass_conservation(U, D):
    Ca = U * dt / dx
    Cd = D * dt / pow(dx, 2)
    c = np.zeros((nt, nx))

    for x in range(nx):
        if x == injection_point:
            c[0, x] = tracer_amount / (dx * width * depth)
        else:
            c[0, x] = 0

    total_mass = np.zeros(nt)
    for t in range(nt - 1):
        for x in range(2, nx - 1):
            component1 = (Cd * (1 - Ca) - Ca / 6 * (pow(Ca, 2) - 3 * Ca + 2)) * c[t, x + 1]
            component2 = (Cd * (2 - 3 * Ca) - Ca / 2 * (pow(Ca, 2) - 2 * Ca - 1)) * c[t, x]
            component3 = (Cd * (1 - 3 * Ca) - Ca / 2 * (pow(Ca, 2) - Ca - 2)) * c[t, x - 1]
            component4 = (Cd * Ca + Ca / 6 * (pow(Ca, 2) - 1)) * c[t, x - 2]
            c[t + 1, x] = c[t, x] + component1 - component2 + component3 + component4
        total_mass[t] = sum(c[t, :])
    return total_mass

c = simulation(U, D)

plt.plot(range(nx), c[0, :])
plt.title('[Figure 1] Temporal evolution of tracer concentration - beginning of simulation')
plt.xlabel('step')
plt.ylabel('c(t)')
plt.show()

plt.plot(range(nx), c[500, :])
plt.title('[Figure 2] Temporal evolution of tracer concentration - middle of simulation')
plt.xlabel('step')
plt.ylabel('c(t)')
plt.show()

plt.plot(range(nx), c[nt - 1, :])
plt.title('[Figure 3] Temporal evolution of tracer concentration - end of simulation')
plt.xlabel('step')
plt.ylabel('c(t)')
plt.show()
#
plt.plot(range(nt), c[:, 20])
plt.title('[Figure 4] Spacial evolution of tracer concentration - beginning of simulation')
plt.xlabel('step')
plt.ylabel('c(x)')
plt.show()

plt.plot(range(nt), c[:, 50])
plt.title('[Figure 5] Spacial evolution of tracer concentration - middle of simulation')
plt.xlabel('step')
plt.ylabel('c(x)')
plt.show()

plt.plot(range(nt), c[:, nx - 2])
plt.title('[Figure 6] Spacial evolution of tracer concentration - end of simulation')
plt.xlabel('step')
plt.ylabel('c(x)')
plt.show()

# Numerical stability
c = simulation(0.001, 1)
plt.plot(range(nx), c[500, :])
plt.title('[Figure 7] Temporal evolution - testing numerical stability for U=%f and D=%d' % (0.001, 1))
plt.xlabel('step')
plt.ylabel('c(t)')
plt.show()

c = simulation(3, 10)
plt.plot(range(nx), c[50, :])
plt.title('[Figure 8] Temporal evolution - testing numerical stability for U=%d and D=%d' % (3, 10))
plt.xlabel('step')
plt.ylabel('c(t)')
plt.show()

# Constant amount

total_mass = simulation_mass_conservation(U, D)
plt.plot(range(nt), total_mass)
plt.title('[Figure 9] Tracer concentration over time - mass conservation law test')
plt.xlabel('timestep')
plt.ylabel('tracer concentration')
plt.show()


