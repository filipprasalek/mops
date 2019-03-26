import matplotlib.pyplot as plt
import numpy as np
import random


def get_random_number():
    return random.normalvariate(mu=0, sigma=1)


def draw_trajectory_plot(particles, steps):
    x = np.zeros((steps, particles))
    y = np.zeros((steps, particles))
    for step in range(1, steps):
        for particle in range(0, particles):
            x[step, particle] = x[step - 1, particle] + get_random_number()
            y[step, particle] = y[step - 1, particle] + get_random_number()

    for particle in range(0, particles):
        plt.plot(x[:, particle], y[:, particle], linewidth=0.5)
    plt.title('Trajectory for {0} particles in {1} steps'.format(particles, steps))
    plt.xlabel('x coordinate')
    plt.ylabel('y coordinate')
    plt.show()


def draw_mean_time_displacement(particles, steps):
    mean_squares = np.zeros(steps)
    x_in_step = np.zeros((steps, particles))
    y_in_step = np.zeros((steps, particles))
    for step in range(1, steps):
        for particle in range(0, particles):
            x_in_step[step, particle] = x_in_step[step - 1, particle] + get_random_number()
            y_in_step[step, particle] = y_in_step[step - 1, particle] + get_random_number()
        mean_squares[step] = np.mean(np.power(x_in_step[step], 2) + np.power(y_in_step[step], 2))
    plt.plot(range(0, steps), mean_squares)
    plt.title('Relationship between mean square of displacement and time for {0} particles'.format(particles))
    plt.xlabel('timestep')
    plt.ylabel('mean square displacement')
    plt.show()


def draw_evolution_of_particles_density(particles, steps):
    x = np.zeros(particles)
    y = np.zeros(particles)
    for step in range(1, steps):
        for particle in range(0, particles):
            x[particle] = x[particle] + get_random_number()
            y[particle] = y[particle] + get_random_number()

    fig, (hist, hist2d) = plt.subplots(1, 2)
    hist.hist(x, bins=20)
    hist.set_title('Density over time for {0} particles in {1} steps - 1D'.format(particles, steps))
    hist.set_xlabel('x coordinate')
    hist.set_ylabel('particles')
    h = hist2d.hist2d(x, y, bins=20)
    hist2d.set_title('Density over time for {0} particles in {1} steps (2D)'.format(particles, steps))
    hist2d.set_xlabel("x coordinate")
    hist2d.set_ylabel("y coordinate")
    colorbar = plt.colorbar(h[3], ax=hist2d)
    colorbar.set_label("particles")
    plt.show()


draw_trajectory_plot(particles=1, steps=10000)
draw_trajectory_plot(particles=1000, steps=100)
draw_mean_time_displacement(particles=10000, steps=100)
draw_evolution_of_particles_density(particles=10000, steps=100)
