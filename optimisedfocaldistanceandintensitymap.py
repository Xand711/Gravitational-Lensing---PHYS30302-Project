# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 16:48:27 2026

@author: annab
"""

"""
Newton Raphson method to work out best distance to reciever intensity wise (on axis so aligned)

Expecting around z = 550AU?  bc of focal length calculations. or well 550 is the min it should be.. 
Attempting to comment this so u guys understand what im trying to do lol 

Best z value can vary a lot but think thats bc of the monte carlo - w more N would be more accurate but 
my laptop is ass and i dont want it to explode sorry - ok is better with higher max impact paramm (umax) but can get unrealistic.. see codeinfo.pdf

+++ assumes the source is coming from far enough away that the rays can be assumed parallel (small angle approx etc)

also to be clear about the netwon raph, bc finding max of the intensity, the equation is the 
derivitive of I(z) = 0, hence the new raph is znew = z - I'/I''
but bc no formula for I specifically, using central difference formlae which are
I' = I(z + dz) = I(z-dz)/2dz
I'' = (I(z+dz) - 2I(z) + I(z-dz))/(dz)^2
and have chosen a value of dz. 
"""
# important stuff first!
import numpy as np
import matplotlib.pyplot as plt
G = 6.674e-11 #m^3kg^-1s^-2
Msun = 1.989e30 #sun mass in kg 
c = 3.0e8 #m/s
AU = 1.496e11 #m
Rsun = 6.96e8 #sun radius in m

def rayangle(u):
    """
    angle by which light ray bends w impact param u (m)
    angle is 4 G M / (c^2 u), radians i think #wikipedia
    """
    return 4 * G * Msun / (c**2 * u)

def rays(N, umax, zobs):
    """
    allegedly monte carlo,,, google helped me expeditiously here but think it works..?
    only works if rays parallel in z (in long distances works i think)
    N, umax, zobs see below (sorry i started w newton raph and worked upwards)
    """
    # r^2 uniform (i think), so r = umax*sqrt(U), where U is uniform(0,1) distribution..?
    r = umax*np.sqrt(np.random.rand(N)) # making random points 
    phi = 2*np.pi*np.random.rand(N) # making random angular distribution - hopefully uniform?

    #cartesian inital positions
    x0 = r*np.cos(phi) 
    y0 = r*np.sin(phi)
    u = r  # impact parammag

    # bending angle around sun for rays 
    angle = rayangle(u)

    # deflection towards origin so get alpha along -x and -y directions ACCORDING TO GOOGLE LOL
    # small angle approx allowed bc of large distances yay
    thetax = -angle*(x0/u)
    thetay = -angle*(y0/u)

    # follow rays z=0->z_obs:
    x = x0 + thetax* zobs
    y = y0 + thetay * zobs

    return x, y # array of where rays hit observer plane

def intensity(zobs, N, umax, aperture):
    """
    on-axis intensity - within aperture , no angular deviations (an approx!)
    zobs = observer distance in m
    N = number of rays simulating
    umax see below
    aperture see below. 
    """
    x, y = rays(N, umax, zobs)
    r = np.sqrt(x**2 + y**2) # pythagoras i guess ?
    I = np.sum(r<aperture) # rays only collected inside size of collector!
    return I

def intensity_avg(z, N, umax, aperture, repeats=20):
    return np.mean([intensity(z, N, umax, aperture) for _ in range(repeats)])


def newtonraphson(z0,dz,tol,maxattempt,N,umax,aperture):
    """
    z0 = inital guess for focal distance
    dz = bottom of differentials... small change in z basically
    tol = how close 2 guesses can be for it to be accepted (tolerance)
    maxattempt = numer of iterations before stopping so my computer doesnt crash lol
    N = number of rays being simulated
    umax = max impact paramter (= perpendicular dis between ray and centre of sun) - cant be less that sun radiys too
    aperture = radius of aperture - so size of reciever i guess - random... 
    
    approximating derivs (found this online! so hope it works lol)
    dI/dz  ->  (I(z+dz) - I(z-dz)) / (2 dz)
    d2I/dz2 -> (I(z+dz) - 2I(z) + I(z-dz)) / dz^2
    z_new = z - (dI/dz) / (d2I/dz2) as per 
    """
    z = z0 # inital guess to start!
    history = [] # so we can plot iteration 

    for attempt in range(maxattempt):
        #intensities at z-dz, z, z+dz
        #Ilow = intensity(z-dz, N, umax, aperture)
        #I0 = intensity(z, N, umax, aperture)
        #Ihigh = intensity(z+dz, N, umax, aperture)
        
        Ilow = intensity_avg(z-dz, N, umax, aperture,20)
        I0 = intensity_avg(z, N, umax, aperture,20)
        Ihigh = intensity_avg(z+dz, N, umax, aperture,20)
        
        
        
        history.append((z, I0)) # for later teehee
        dI = (Ihigh-Ilow)/(2*dz)
        d2I = (Ihigh - 2 * I0 + Ilow)/(dz**2)

        if d2I == 0:
            print("#flop bc divide by zero error god save me")
            break

        znew = z-(dI/d2I)

        print(
            f"Iter {attempt}: z = {z/AU:.1f} AU, I = {I0}, "
            f"dI = {dI:.3e}, d2I = {d2I:.3e}, znew = {znew/AU:.1f} AU"
        )

        # if in tolerance we're done ig
        if abs(znew-z) < tol:
            z = znew
            history.append((z, intensity(z, N, umax, aperture)))
            print("YIPPEE")
            break

        z = znew

    return z, history


#ACTUAL CODE STARTS HERE - have vaguely had a fiddle w these params but idk. 
zopt = 0
z0guess = 550 * AU
dz = 10*AU
tol = 0.1*AU
maxattempt = 10
N = 150000
umax = 3*Rsun # originally had this as 2, but if u increaqse it it gets more consistent results
aperture = 1e9


zopt, history = newtonraphson(z0guess,dz,tol,maxattempt,N,umax,aperture)
print(f"\nOptimised focal distance: {zopt/AU:.1f} AU") # in AU

# plotting each iteration
z = np.array([a[0] for a in history]) / AU
I = np.array([a[1] for a in history]) # 2D array hate club this took me forever

plt.figure(figsize=(6, 4))
plt.plot(z, I, marker='x',color='purple')
plt.xlabel("z (AU)")
plt.ylabel("Intensity (counts)")
plt.title("Newton–Raphson iterations")
plt.show()

#going to try and make a map like in emily's code bc i feel like with all that ray tracing nonsense 
# i may as well try and use it visually bc i hated coding that #releaseme


def intensitymap(x, y, bins, extent, zobs):
    """
     this is dubious ngl. 
    x,y = where ray ends up
    bins = number of spots in x and y
    extent = size map covers in x and y from - to + 
    """
    H, xedges, yedges = np.histogram2d(x, y,bins,range=[[-extent, extent], [-extent, extent]])
    # H is 2D histogram of rays, edges are that of bins
    
    extentplot = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    plt.figure(figsize=(6, 5))
    plt.imshow(H.T, origin='lower', extent=np.array(extentplot) / 1e9, # Gm
        cmap='inferno',aspect='equal')
    plt.colorbar(label="Relative intensity (counts)")
    plt.xlabel("x (Gm)")
    plt.ylabel("y (Gm)")
    plt.title(f"Intensity map at z = {zobs / AU:.1f} AU")
    plt.tight_layout()
    plt.show()


finalz = zopt # from part 1. can test w 550 AU to check works also
#finalz = 550
xmap, ymap = rays(N, umax, finalz)
bins = 250
extent = 2e10
intensitymap(xmap, ymap, bins, extent, finalz)
