"""
Solving the two-dimensional diffusion equation

Example acquired from https://scipython.com/book/chapter-7-matplotlib/examples/the-two-dimensional-diffusion-equation/
"""

import numpy as np
import matplotlib.pyplot as plt
from output import create_plot, output_plots


def do_timestep(u_nm1, u, D, dt, dx2, dy2):
    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1] = u_nm1[1:-1, 1:-1] + D * dt * (
            (u_nm1[2:, 1:-1] - 2 * u_nm1[1:-1, 1:-1] + u_nm1[:-2, 1:-1]) / dx2
            + (u_nm1[1:-1, 2:] - 2 * u_nm1[1:-1, 1:-1] + u_nm1[1:-1, :-2]) / dy2)

    u_nm1 = u.copy()
    return u_nm1, u

def solve(dx=0.1, dy=0.1, D=4.0):
    # Plate size, mm
    w = h = 10.

    # Thermal diffusivity of steel, mm^2/s
    T_cold = 300
    T_hot = 700

    # Number of mesh points in X and Y directions
    nx, ny = int(w / dx), int(h / dy)

    # Time step calculations
    dx2, dy2 = dx * dx, dy * dy
    dt = dx2 * dy2 / (2 * D * (dx2 + dy2))

    print("dt = {}".format(dt))

    # Initial temperature distribution
    u0 = T_cold * np.ones((nx, ny))
    u = u0.copy()

    # Initial conditions - circular disc at center
    r, cx, cy = 2, 5, 5
    r2 = r ** 2
    for i in range(nx):
        for j in range(ny):
            p2 = (i * dx - cx) ** 2 + (j * dy - cy) ** 2
            if p2 < r2:
                u0[i, j] = T_hot

    # Number of timesteps
    nsteps = 101
    n_output = [0, 10, 50, 100]

    # Initialize figure
    fig = plt.figure()
    fig_counter = 0
    im = None

    # Time loop
    for n in range(nsteps):
        u0, u = do_timestep(u0, u, D, dt, dx2, dy2)

        # Plot at specific timesteps
        if n in n_output:
            fig_counter += 1
            im = create_plot(u, T_cold, T_hot, n, dt, fig_counter, fig)

    # Output all plots with color bar
    output_plots(fig, im)
