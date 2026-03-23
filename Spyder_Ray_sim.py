# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 12:20:40 2026

@author: xandb
"""

# Enter " %matplotlib qt " into terminal before running to make interactive plots.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.random import uniform as rand
from numpy.linalg import norm as mag


global G, c, R_solar
G  = 4 * np.pi**2 # In AU, solar masses and yr
c = 63241.1 # AU yr^-1
R_solar = 0.046504 # In AU

# Just normalises vectors (nothing special)
def normalise(vec):
    return vec/mag(vec)

# Using object classes for all elements of the simulation
class Transmitter:
    def __init__(self, position, wavlen, direction):
       self.pos = np.array([position[0],position[1],position[2]], dtype=np.float64)
       self.wavlen = wavlen
       self.direction = normalise(direction - self.pos)
       
    # Initialises the rays from the transmitter
    def emit_signal(self, n_ray, spread_angle, lens):
        rays = [None] * n_ray
        for i in range(n_ray):
            direc = self.random_direc(spread_angle, self.direction)
            rays[i] = Ray(self.wavlen, self.pos.copy(), direc, lens.pos.copy())
        return rays

    # Generate a random direction in a cone
    def random_direc(self, spread, d):
        # Sample random direction within cone
        theta = rand(0, 2*np.pi)
        phi = np.arccos(1 - rand(0,1)*(1 - np.cos(spread)))

        x = np.sin(phi)*np.cos(theta)
        y = np.sin(phi)*np.sin(theta)
        z = np.cos(phi)
        
        # Next we need to point this cone in the direction of the transmitter
        # Generate arbitrary perpendicular vectors to d
        if abs(d[2]) < 0.999:
            temp = np.array([0, 0, 1])
        else:
            temp = np.array([1, 0, 0])

        u = np.cross(temp, d)
        u /= mag(u)

        v = np.cross(u, d) 

        # Return random unit vector within defined cone
        return normalise(u*x + v*y + d*z)


class Ray:
    def __init__(self, wavlen, position, direction, lens_pos):
        self.wavlen = wavlen
        self.pos = np.array([position[0],position[1],position[2]], dtype=np.float64)
        self.direction = normalise(direction)
        self.history = [self.pos.copy()]
        self.lens_pos = lens_pos
        self.calc_lens_dist()
        self.path_len = 0
    
    # Just updates the position of the ray to simplify code elsewhere
    def update_pos(self, change):
        self.pos += change
        self.history.append(self.pos.copy())
        self.path_len += change
        self.calc_lens_dist()

    # Calculates distance to the lens from the rays current position
    def calc_lens_dist(self):
        self.r = self.pos - self.lens_pos
        self.r_mag = mag(self.r)

    # Simple straight line propagation
    def propagate(self, dt):
        self.update_pos(self.direction * dt * c)
    
    # More complex propagation that updates the position and direction of the ray based on the gravitational field of the star
    def step(self, dt, lens):
        acc = -2 * 10 * (G* lens.mass / self.r_mag**3) * self.r

        self.direction += (acc / c) * dt
        self.direction = normalise(self.direction)
        self.update_pos(self.direction * dt * c)
    # Calculates phase of ray at a given point (not currently functional)
    def get_phase(self):
        return 2 * np.pi * self.path_length / self.wavlen
        
class Lensing_object:
    def __init__(self, mass, position):
        self.mass = mass
        self.pos = position

class Receiver:
    def __init__(self, position, radius):
        self.pos = np.array([position[0],position[1],position[2]], dtype=np.float64)
        self.radius = radius
        self.detected_rays = []
    
    # Adds ray to a list of detected rays if in range
    def detect(self, ray):
        if np.linalg.norm(ray.pos - self.pos) < self.radius:
            self.detected_rays.append(ray)
            
# This function actually initialises and runs the simulation
def run_sim(trans, rec, lens, n_rays, spread):
    # Initialise rays
    rays = trans.emit_signal(n_rays, spread, lens)
    # Time step within the influence of the sun
    dt_lens = 1e-6

    for ray in rays:
        # Propagates rays towards the lens whilst outside of the its infulence
        while ray.r_mag > 100 * R_solar:
            dt = min(0.01, (ray.r_mag - 100*R_solar)/(2*c))
            ray.propagate(dt)

        ray.entry_dir = ray.direction
        
        # Propagates the ray within the influence of the lens along curved path
        while ray.r_mag <= 100 * R_solar:
            ray.step(dt_lens, lens)

        ray.exit_dir = ray.direction
        ray.interstellar_idx = len(ray.history)
    
    # Propagates rays from the lens directly to the receiver (Direction doesn't change here so we don't need to simulate everything inbetween)
    for ray in rays:
        dt = (abs(rec.pos[1] - ray.pos[1])/c) * ray.direction[1]
        ray.propagate(dt)
        rec.detect(ray)

    return rays

# Calculates the difference in angle between two vectors (Used mainly to calculate the deflection of rays by the sun)
def angle_calc(v_1, v_2):
    cos_theta = np.dot(v_1, v_2) / (np.linalg.norm(v_1) * np.linalg.norm(v_2))
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    return angle_rad

sat_pos = np.array([0,-30,0])
sun_pos = np.array([0,5,0])
n_ray = 500
detection_range = 0.3
# Spread is calculated to limit the rays to the area of influence around the sun. (Don't want to simulate rays that aren't bent by grav lensing)
spread = np.arctan(10*R_solar/mag(sat_pos - sun_pos))

home = Receiver([0,550,0],detection_range)
sun = Lensing_object(1,sun_pos)
sat = Transmitter(sat_pos,5,sun_pos)


# This runs the simulation
r = run_sim(sat,home,sun,n_ray,spread)

n = np.random.randint(0,n_ray)

theta = angle_calc(r[n].entry_dir, r[n].exit_dir)
print(r[n].entry_dir, r[n].exit_dir)
print(theta)

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection='3d')

# Plots rays
for ray in r[::40]:
    prev = ray.history[0]
    for vec in ray.history[1::15]:
        ax.plot([vec[0],prev[0]],[vec[1],prev[1]],[vec[2],prev[2]], 'r')
        prev = vec

ax.scatter(sat.pos[0],sat.pos[1],sat.pos[2])
#ax.scatter(home.pos[0],home.pos[1],home.pos[2])
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")


u = np.linspace(0, 2*np.pi, 50)
v = np.linspace(0, np.pi, 50)

# Uncomment section to view the part of the simulation that you want to see

'''
# View Satellite

vol = 0.05

ax.set_xlim(sat.pos[0]-vol, sat.pos[0]+vol)
ax.set_ylim(sat.pos[1]-vol, sat.pos[1]+vol)
ax.set_zlim(sat.pos[2]-vol, sat.pos[2]+vol)
'''


# View Sun

X = sun_pos[0] + R_solar * np.outer(np.cos(u), np.sin(v))
Y = sun_pos[1] + R_solar * np.outer(np.sin(u), np.sin(v))
Z = sun_pos[2] + R_solar * np.outer(np.ones_like(u), np.cos(v))

ax.plot_surface(X, Y, Z, color='orange', alpha=0.8)

vol = 0.3

ax.set_xlim(sun.pos[0]-vol, sun.pos[0]+vol)
ax.set_ylim(sun.pos[1]-vol, sun.pos[1]+vol)
ax.set_zlim(sun.pos[2]-vol, sun.pos[2]+vol)


'''
# View receiver

u = np.linspace(0, 2*np.pi, 50)
v = np.linspace(0, np.pi, 50)

X = home.pos[0] + detection_range * np.outer(np.cos(u), np.sin(v))
Y = home.pos[1] + detection_range * np.outer(np.sin(u), np.sin(v))
Z = home.pos[2] + detection_range * np.outer(np.ones_like(u), np.cos(v))

ax.plot_surface(X, Y, Z, color='orange', alpha=0.3)
vol = 1

ax.set_xlim(home.pos[0]-vol,home.pos[0]+vol)
ax.set_ylim(home.pos[1]-vol,home.pos[1]+vol)
ax.set_zlim(home.pos[2]-vol,home.pos[2]+vol)
'''

plt.show()