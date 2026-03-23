# Gravitational-Lensing---PHYS30302-Project

## About the 3D simulation
1. The transmitter creates a group of rays distributed randomly in a cone pointed towards the sun. (Initiated by transmitter.emit_signal())
2. The rays are propagated towards the sun in a straight line using ray.propagte .
3. When the ray is within 100 solar radii of the sun the path of the ray is calculated using euler method (will update to runge kutta). Rays are treated as particles in a potential.
4. Once the ray is outside of 100 solar radii it is propagated direclty to the reciever (Nothing interesiting happens in this span the ray so the ray just moves straight)

Distance units are in AU, mass in solar masses, and time in years
