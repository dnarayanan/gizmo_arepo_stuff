from astropy.cosmology import Planck15
from astropy import units as u
import numpy as np

Omega_b = 0.048
rho_dens = Omega_b*Planck15.critical_density(0)

L = 12.5*u.Mpc
N_particles = (2.**11)**3.


mass_in_box = rho_dens*(L**3)
mass_per_particle = mass_in_box/N_particles

print(mass_per_particle.to(u.Msun))

