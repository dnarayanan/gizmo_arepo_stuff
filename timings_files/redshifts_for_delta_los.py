# -*- coding: utf-8 -*-
"""
Created on Fri May 20 14:44:49 2016

@author: jaguirre
"""

import numpy as np
import pylab as plt
from astropy import units as u
from astropy import constants as c
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value

#MODIFIABLE HEADER
BOXSIZE = 50./cosmo.h*u.Mpc

# Start at redshift 0 and find the redshifts for each BOXSIZE comoving Mpc increase
# thereafter
zs = [0]
zlast = zs[0]
dc = cosmo.comoving_distance(zlast)
dcs = [dc.value]

while (zlast < 30):
    dc += BOXSIZE 
    dcs.append(dc.value)
    zlast = z_at_value(cosmo.comoving_distance,dc)
    zs.append(zlast)
    print zlast

zs = np.array(zs)
dcs = np.array(dcs)*u.Mpc
a = 1./(1+zs)

a.sort()

f = open('scale_factors_for_sims.txt','w')
for a1 in a:
    f.write('%10.5f \n' % (a1))
    
f.close()
