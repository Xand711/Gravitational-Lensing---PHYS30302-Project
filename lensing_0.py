# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:50:02 2026

@author: emily
"""

import numpy as np
import matplotlib.pyplot as plt


# Adjustable parameter
theta_E = 0.75        # Einstein radius, can be calculated if need
N = 500              # grid resolution
extent = 2.0         # image plane spans [-extent, extent]

# Lens position (in image plane)
lens_x = -0.4
lens_y = 0.2

# Source position (in source plane)
src_x = 0.5
src_y = 0.1

# Source size
sigma = 0.1


# Lensing functions
def make_image_plane(N, extent):
    x = np.linspace(-extent, extent, N)
    y = np.linspace(-extent, extent, N)
    return np.meshgrid(x, y)


def lens_equation(theta_x, theta_y, theta_E, lens_x, lens_y):
    # Shift coordinates so lens is not forced to be at (0,0)
    dx = theta_x - lens_x
    dy = theta_y - lens_y

    r2 = dx**2 + dy**2

    beta_x = theta_x - (theta_E**2) * dx / r2
    beta_y = theta_y - (theta_E**2) * dy / r2

    return beta_x, beta_y


def source_brightness(beta_x, beta_y, src_x, src_y, sigma):
    return np.exp(-((beta_x - src_x)**2 + (beta_y - src_y)**2) / (2 * sigma**2))


def ray_shoot(theta_x, theta_y, theta_E, lens_x, lens_y, src_x, src_y, sigma):
    beta_x, beta_y = lens_equation(theta_x, theta_y, theta_E, lens_x, lens_y)
    return source_brightness(beta_x, beta_y, src_x, src_y, sigma)



# Plotting
def plot_lensed_image(image, lens_x, lens_y, src_x, src_y, extent):
    plt.figure(figsize=(6, 6))
    plt.imshow(
        image,
        cmap='inferno',
        origin='lower',
        extent=[-extent, extent, -extent, extent]
    )

    # Plot lens (cross)
    plt.scatter(lens_x, lens_y, color='cyan', s=80, marker='x', label='Lens')

    # Plot source (circle)
    plt.scatter(src_x, src_y, color='white', s=60, marker='o', label='Source')

    plt.legend(loc='upper right')
    plt.title("Gravitational Lensing (Manual Lens & Source Positions)")
    plt.xlabel("θ_x")
    plt.ylabel("θ_y")
    plt.show()



# Run simulation
theta_x, theta_y = make_image_plane(N, extent)

image = ray_shoot(
    theta_x, theta_y,
    theta_E,
    lens_x, lens_y,
    src_x, src_y,
    sigma
)

plot_lensed_image(image, lens_x, lens_y, src_x, src_y, extent)