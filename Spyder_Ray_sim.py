# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 12:20:40 2026

@author: xandb
"""

# Enter %matplotlib qt into terminal before running to make interactive plots.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.random import uniform as rand
from numpy.linalg import norm as mag


global G, c, R_solar
G  = 4 * np.pi**2 # In AU, solar masses and yr
c = 63241.1 # AU yr^-1
R_solar = 0.046504 # In AU

def normalise(vec):
    return vec/mag(vec)

class Transmitter:
    def __init__(self, position, wavlen, direction):
       self.pos = np.array([position[0],position[1],position[2]], dtype=np.float64)
       self.wavlen = wavlen
       self.direction = normalise(direction - self.pos)
        
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

    def update_pos(self, change):
        self.pos += change
        self.history.append(self.pos.copy())
        self.path_len += change
        self.calc_lens_dist()

    def calc_lens_dist(self):
        self.r = self.pos - self.lens_pos
        self.r_mag = mag(self.r)

    # Simple straight line propagation
    def propagate(self, dt):
        self.update_pos(self.direction * dt * c)
        
    def step(self, dt, lens):
        acc = -2 * 10 * (G* lens.mass / self.r_mag**3) * self.r

        self.direction += (acc / c) * dt
        self.direction = normalise(self.direction)
        self.update_pos(self.direction * dt * c)

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

    def detect(self, ray):
        if np.linalg.norm(ray.pos - self.pos) < self.radius:
            self.detected_rays.append(ray)

    def signal_strength(self):
        return len(self.detected_rays)
    
def run_sim(trans, rec, lens, n_rays, spread):
    rays = trans.emit_signal(n_rays, spread, lens)
    dt_lens = 1e-6
    dt_post_lens = 0.1

    for ray in rays:
        while ray.r_mag > 100 * R_solar:
            dt = min(0.01, (ray.r_mag - 100*R_solar)/(2*c))
            ray.propagate(dt)

        ray.entry_dir = ray.direction
            
        while ray.r_mag <= 100 * R_solar:
            ray.step(dt_lens, lens)

        ray.exit_dir = ray.direction
        ray.interstellar_idx = len(ray.history)
                    
    for ray in rays:
        while ray.pos[1] < rec.pos[1]:
            ray.propagate(dt_post_lens)
            #rec.detect(ray)

    return rays

def angle_calc(v_1, v_2):
    cos_theta = np.dot(v_1, v_2) / (np.linalg.norm(v_1) * np.linalg.norm(v_2))
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    return angle_rad

sat_pos = np.array([0,-30,0])
sun_pos = np.array([0,5,0])

n_ray = 20

spread = np.arctan(10*R_solar/mag(sat_pos - sun_pos))

rec = Receiver([0,20000,0],2)
sun = Lensing_object(1,sun_pos)

sat = Transmitter(sat_pos,5,sun_pos)

r = run_sim(sat,rec,sun,n_ray,spread)

n = np.random.randint(0,n_ray)

theta = angle_calc(r[n].entry_dir, r[n].exit_dir)
print(r[n].entry_dir, r[n].exit_dir)
print(theta)

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection='3d')

for ray in r:
    prev = ray.history[0]
    for vec in ray.history[1:ray.interstellar_idx:20]:
        ax.plot([vec[0],prev[0]],[vec[1],prev[1]],[vec[2],prev[2]], 'r')
        prev = vec

ax.scatter(sat.pos[0],sat.pos[1],sat.pos[2])
#ax.scatter(rec.pos[0],rec.pos[1],rec.pos[2])
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")


u = np.linspace(0, 2*np.pi, 50)
v = np.linspace(0, np.pi, 50)

X = sun_pos[0] + R_solar * np.outer(np.cos(u), np.sin(v))
Y = sun_pos[1] + R_solar * np.outer(np.sin(u), np.sin(v))
Z = sun_pos[2] + R_solar * np.outer(np.ones_like(u), np.cos(v))

ax.set_xlim(-0.3, 0.3)
ax.set_ylim(4.7, 5.3)
ax.set_zlim(-0.3, 0.3)

ax.plot_surface(X, Y, Z, color='orange', alpha=0.8)

X = sun_pos[0] + 100*R_solar * np.outer(np.cos(u), np.sin(v))
Y = sun_pos[1] + 100*R_solar * np.outer(np.sin(u), np.sin(v))
Z = sun_pos[2] + 100*R_solar * np.outer(np.ones_like(u), np.cos(v))

#ax.plot_surface(X, Y, Z, color='orange', alpha=0.3)

plt.show()