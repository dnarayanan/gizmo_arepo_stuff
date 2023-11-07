import caesar
import numpy as np
import pdb

def contamination_check(obj,ad,halonum):

    pos = obj.halos[halonum].pos
    radius = obj.halos[halonum].radii['total_r20']

    
    lowres_particle_pos = ad['PartType2','Coordinates'].in_units('kpccm')
    lowres_particle_mass = ad['PartType2','Masses'].in_units('Msun')
    
    lowres_rad = np.sqrt( (pos[0]-lowres_particle_pos[:,0])**2 + (pos[1]-lowres_particle_pos[:,1])**2 + (pos[2]-lowres_particle_pos[:,2])**2)
    w = np.where(lowres_rad <= radius)[0]
    total_contam_mass = np.sum(lowres_particle_mass[w])

    return(total_contam_mass)

