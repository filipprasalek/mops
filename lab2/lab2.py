import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import cm


def is_border(i, j, area_size):
    return i == 0 or i == area_size - 1 or j == 0 or j == area_size - 1


def is_heater(i, j, area_size, heater_size):
    low_boundary = area_size / 2 - heater_size / 2
    high_boundary = area_size / 2 + heater_size / 2
    return (low_boundary <= i < high_boundary) and (low_boundary <= j < high_boundary)


def initialize_area(area_size, heater_size, border_temp, heater_temp):
    area = np.zeros((area_size, area_size))
    for i in range(area_size):
        for j in range(area_size):
            if is_border(i, j, area_size):
                area[i, j] = border_temp
            elif is_heater(i, j, area_size, heater_size):
                area[i, j] = heater_temp
            else:
                area[i, j] = 20
    return area


def new_x_temp(x, y, previous_area, metal, denominator, dt):
    return (metal["K"] * dt * (
            previous_area[x + 1, y] - 2 * previous_area[x, y] + previous_area[x - 1, y])) / denominator


def new_y_temp(x, y, previous_area, metal, denominator, dt):
    return (metal["K"] * dt * (
            previous_area[x, y + 1] - 2 * previous_area[x, y] + previous_area[x, y - 1])) / denominator


def draw_3d_plot(area, dx, dy, time, size):
    x = np.multiply(np.arange(0, size), dx)
    y = np.multiply(np.arange(0, size), dy)
    xx, yy = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_title("Simulation time: {0}s, spatial step dx = dy = {1}m".format(time, dx, dy))
    surf = ax.plot_surface(xx, yy, area, cmap=cm.coolwarm, linewidth=0, antialiased=True)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    colorbar = fig.colorbar(surf, shrink=0.5, aspect=5)
    colorbar.set_label("temp [C]")
    plt.show()


def first_boundary_conndition_simulation(steps, initial_area, metal, dx, dy, dt, size):
    denom_x = metal["cw"] * metal["ro"] * dx ** 2
    denom_y = metal["cw"] * metal["ro"] * dy ** 2
    area_in_step = initial_area

    for i in range(1, steps):
        prev_area = area_in_step[i - 1]
        for x in range(1, size - 1):
            for y in range(1, size - 1):
                if is_heater(x, y, size, h_size):
                    continue
                x_temp = new_x_temp(x, y, prev_area, metal, denom_x, dt)
                y_temp = new_y_temp(x, y, prev_area, metal, denom_y, dt)
                area_in_step[i][x, y] = prev_area[x, y] + x_temp + y_temp
    draw_3d_plot(area_in_step[steps - 1], dx, dy, dt * steps, size)


metal_params = {
    "alumina": {"ro": 2700, "cw": 900, "K": 237},
    "cooper": {"ro": 8920, "cw": 380, "K": 401},
    "steel": {"ro": 7860, "cw": 450, "K": 58}
}

a_size = 6
h_size = 2
border_t = 10  # [C]
heater_t = 80  # [C]

steps = 100
dt = 0.01  # [s]
dx = dy = 0.01  # [m]

area_states = [initialize_area(a_size, h_size, border_t, heater_t)] * steps
first_boundary_conndition_simulation(steps, area_states, metal_params["alumina"], dx, dy, dt, a_size)
first_boundary_conndition_simulation(steps, area_states, metal_params["cooper"], dx, dy, dt, a_size)
first_boundary_conndition_simulation(steps, area_states, metal_params["steel"], dx, dy, dt, a_size)
