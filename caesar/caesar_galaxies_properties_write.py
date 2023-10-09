# import SPHGR, numpy, pylab, and glob2
import yt
import numpy as np
from glob2 import glob
import ipdb,pdb
from scipy.optimize import newton
from astropy.cosmology import Planck13
from astropy.io import ascii
from astropy.table import Table
import caesar
import sys
sys.path.insert(0,'inputfiles/')





#from caesar_galaxies_properties_write_input import *
from scipy import spatial



def get_cgf(obj,index,dd,radius):
        #radius is in kpc
    slist = obj.galhalo[INDEX].slist
    glist = obj.galhalo[INDEX].glist

    if len(glist) > 0:

        gpos = dd[ ('PartType0', 'particle_position')][glist].in_units('kpc').d
        spos = dd[ ('PartType4', 'particle_position')][slist].in_units('kpc').d
        
        sages = dd[('PartType4', 'StellarFormationTime')]
        scalefactor = sages[obj.galhalo[INDEX].slist]
        formation_z = (1./scalefactor)-1.
        yt_cosmo = yt.utilities.cosmology.Cosmology()
        stellar_formation_age = yt_cosmo.t_from_z(formation_z).in_units('Gyr')
        simtime = yt_cosmo.t_from_z(obj.simulation.redshift).in_units('Gyr') #what is the age of the Universe right now?
        stellar_ages = (simtime -stellar_formation_age).in_units('Myr')


        gtree = spatial.cKDTree(gpos)
        stree = spatial.cKDTree(spos)

        nearby_young = []
        naked_young = 0.
        nearby_old = []
        naked_old = 0.
        

        for i in range(len(spos[stellar_ages <= 50])):
            neighbors = gtree.query_ball_point(spos[i],radius)
            nearby_young.append(len(neighbors))

            if (len(neighbors) == 0) or (np.sum(dd[('PartType0','Metallicity_00')][neighbors]) == 0):
                naked_young += 1

        #old stars
        for i in range(len(spos[stellar_ages >= 1000])):
            neighbors = gtree.query_ball_point(spos[i],radius)
            nearby_old.append(len(neighbors))

            if (len(neighbors) == 0) or (np.sum(dd[('PartType0','Metallicity_00')][neighbors]) == 0):
                naked_old += 1



        #cgf = np.sum(nearby)/float(np.sum(slist))
        cgf1 = float(naked_young)/len(slist)
        cgf2 = float(naked_old)/len(slist)


    else:
        cgf1,cgf2 = 0,0

    return cgf1,cgf2 #np.sum(nearby)/float(np.sum(slist))

def get_sfr(obj,index):

    #ages
    yt_cosmo = yt.utilities.cosmology.Cosmology()
    simtime = yt_cosmo.t_from_z(obj.simulation.redshift).in_units('Gyr') #what is the age of the Universe right now?

    
    galaxy_gas = obj.galhalo[INDEX].glist
    sindexes = obj.galhalo[INDEX].slist
        

    #needs to be fixed
    
    dd = obj.yt_dataset.all_data()
    sages = dd[('PartType4', 'StellarFormationTime')]
    scalefactor = sages[obj.galhalo[INDEX].slist]


    formation_z = (1./scalefactor)-1.

    smass = dd[('PartType4', 'Masses')]
    stellar_mass = smass[obj.galhalo[INDEX].slist].in_units('Msun')
    
    #xz = newton_raphson(f,simtime,obj.redshift+0.5,tol=0.01)
    #np.savez('junk.npz',simtime=simtime)
    

    yt_cosmo = yt.utilities.cosmology.Cosmology()
    stellar_formation_age = yt_cosmo.t_from_z(formation_z).in_units('Gyr')
    stellar_ages = simtime -stellar_formation_age
    w50 = np.where(stellar_ages.in_units('Myr').value < 50)[0]


    '''
    #calculate the SFR for 50 Myr intervals
    xz_50 = newton(f2_50,obj.redshift+1)
    

    print '=========='
    print 'redshift is ',obj.redshift
    print 'newton_raphson derived redshift is for xz_50: ',xz_50
    print "newton_raphson derived delta t (for xz_50) = ",Planck13.age(obj.redshift).value-Planck13.age(xz_50).value
    w = np.where(formation_z <= xz_50)[0]
    '''

    sfr_50 = np.sum(stellar_mass[w50])/50.e6



   
    return float(sfr_50.value),np.median(stellar_ages).value
    

def f2_50(formation_z):
    simtime = np.load('junk.npz')
    simtime = float(simtime['simtime'])

    print formation_z
    print (simtime-Planck13.age(formation_z).value)
    return (simtime-Planck13.age(formation_z).value)-0.05


def get_metal_mass_radius(obj,dd,INDEX):
    galaxy_gas_ids = obj.galhalo[INDEX].glist

    metal_mass = dd[('PartType0','Masses')][galaxy_gas_ids].in_units('Msun')* dd[('PartType0','Metallicity_00')][galaxy_gas_ids]

    gas_mass = dd[('PartType0','Masses')][galaxy_gas_ids].in_units('Msun')
    
    coords = dd[("PartType0","Coordinates")].in_units('kpccm')
    galaxy_center = obj.galhalo[INDEX].pos.in_units('kpccm')
    radii = np.sqrt( (coords[:,0]-galaxy_center[0])**2. +
                     (coords[:,1]-galaxy_center[1])**2. +
                     (coords[:,2]-galaxy_center[2])**2. )


    metal_mass = gas_mass
    r = 0.
    cumulative_mass = 0
    half_mass = 0.5*np.sum(metal_mass)
    for i in range(0,len(metal_mass)):
        cumulative_mass += metal_mass[i]
        if cumulative_mass >= half_mass:
            r = radii[i]
            break



    return r
    





#========================================
#MAIN CODE
#========================================


# query all available member files
MEMBERS = np.sort(glob('%s/%s/output/Groups/caesar*.hdf5' % (BASEDIR,SNAPPRE)))

#this try/except sequence to determine the lastsnapshot is in place
#because the first option is for caesar files that were run in
#driver.py mode, while the second is for files that were just CLI
#caesar'd (which oddly have different file name structures)
try:  LASTSNAP = int(MEMBERS[-1][-16:-12])
except ValueError,e: LASTSNAP = int(MEMBERS[-1][-8:-5])
# set the galaxy of interest starting index

num_local_halos = len(LOCAL_HALOS)

for lh in LOCAL_HALOS:
    INDEX = lh

# create empty lists to hold Mstar and z
    stellar_masses, gas_masses, h2_masses,hI_masses,halo_masses = [],[],[],[],[]
    redshifts, hmr, fmr, instsfr, sfr_50,metallicity = [],[],[],[],[],[]
    r_baryon,r_h_baryon,r_dm,r_h_dm,r_gas,r_h_gas,r200,r200c,rstellar,r_h_stellar,r_virial = [],[],[],[],[],[],[],[],[],[],[]
    snapshotnames,snaps = [],[]
    metal_radius_list = []
    dm_masses = []
    xpos,ypos,zpos = [],[],[]
    stellar_ages = []
    central_stellar_masses=[] #this is the stellar mass of the central; if
                          #we haven't engaged HALOS=True, then this
                          #will obviously just be the same thing as
                          #'stellar_masses'
    central_sfr = []
    cgf1,cgf2 = [],[]

# cycle through all MEMBERS in reverse
    for i in reversed(range(0,len(MEMBERS))):
        
#for i in reversed(range(180,192))
        # load the current member file
        print 'loading %s' % MEMBERS[i]
        obj = caesar.load(MEMBERS[i])
        ds = yt.load(BASEDIR+SNAPPRE+'output/'+obj.simulation.basename)
        obj.yt_dataset = ds
        

        dd = ds.all_data()
        sages = dd[('PartType4', 'StellarFormationTime')]
        
        snapshotnames.append(obj.simulation.basename)

        if HALOS == False:
            obj.galhalo = obj.galaxies
            halo_masses.append(-1)
        else:
            obj.galhalo = obj.halos
            halo_masses.append(obj.galhalo[INDEX].masses['dm'])

        #temporary holder until we get it sorted
        metal_radius = ds.quan(-1,'kpc')
        #metal_radius = get_metal_mass_radius(obj,dd,INDEX)
        metal_radius_list.append(metal_radius.value)
            
    # get the different masses
        stellar_masses.append(obj.galhalo[INDEX].masses['stellar'])
        gas_masses.append(obj.galhalo[INDEX].masses['gas'])
        h2_masses.append(obj.galhalo[INDEX].masses['H2'])
        hI_masses.append(obj.galhalo[INDEX].masses['HI'])
        
        if HALOS == True:

            try: central_stellar_masses.append(obj.galhalo[INDEX].galaxies[0].masses['stellar'])
            except IndexError:
                 central_stellar_masses.append(-1)
        else:
            try: central_stellar_masses.append(obj.galhalo[INDEX].masses['stellar'])
            except IndexError:
                central_stellar_masses.append(-1)
        

        redshifts.append(obj.simulation.redshift)
        metallicity.append(obj.galhalo[INDEX].metallicity)
        
        hmr.append(obj.galhalo[INDEX].radii['stellar_half_mass'])
        fmr.append(obj.galhalo[INDEX].radii['stellar'])

        r_baryon.append(obj.galhalo[INDEX].radii['baryon'])
        r_h_baryon.append(obj.galhalo[INDEX].radii['baryon_half_mass'])
        r_gas.append(obj.galhalo[INDEX].radii['gas'])
        r_h_gas.append(obj.galhalo[INDEX].radii['gas_half_mass'])
        r200.append(obj.galhalo[INDEX].radii['r200'])
        r200c.append(obj.galhalo[INDEX].radii['r200c'])
        rstellar.append(obj.galhalo[INDEX].radii['stellar'])
        r_h_stellar.append(obj.galhalo[INDEX].radii['stellar_half_mass'])
        r_virial.append(obj.galhalo[INDEX].radii['virial'])


        instsfr.append(obj.galhalo[INDEX].sfr)

        dum_sfr_50,dum_stellar_ages = get_sfr(obj,INDEX)
        sfr_50.append(dum_sfr_50)
        stellar_ages.append(dum_stellar_ages)
        
        dum_cgf1,dum_cgf2 = get_cgf(obj,INDEX,dd,0.25)
        cgf1.append(dum_cgf1)
        cgf2.append(dum_cgf2)
        

        xpos.append(obj.galhalo[INDEX].pos[0].in_units('code_length').value)
        ypos.append(obj.galhalo[INDEX].pos[1].in_units('code_length').value)
        zpos.append(obj.galhalo[INDEX].pos[2].in_units('code_length').value)

        #find the snapnum
        end = obj.simulation.basename.find('.hdf5')
        snap = int(obj.simulation.basename[end-3:end])
        snaps.append(snap)

 #snaps.append(LASTSNAP-(len(redshifts)-1))
        
        # set the new INDEX value to be used in the previous snapshot.  If
        # a -1 is encountered that means there was no progenitor and we
        # can break from the loop.
        
        if FORCE_PROGEN_INDEX == True: 
            INDEX = FORCED_PROGEN_INDEX_VALUE 
        else:
            try: INDEX = obj.galhalo[INDEX].progen_index
            except AttributeError:
                INDEX=FORCED_PROGEN_INDEX_VALUE


        phys_table = Table([snapshotnames[::-1],
                            snaps[::-1],
                            redshifts[::-1],
                            instsfr[::-1],
                            sfr_50[::-1],
                            stellar_masses[::-1],
                            gas_masses[::-1],
                            h2_masses[::-1],
                            hI_masses[::-1],
                            hmr[::-1],
                            fmr[::-1],
                            metallicity[::-1],
                            xpos[::-1],
                            ypos[::-1],
                            zpos[::-1],
                            halo_masses[::-1],
                            central_stellar_masses[::-1],
                            metal_radius_list[::-1],
                            cgf1[::-1],
                            cgf2[::-1]],
                           names=['snapname',
                                  'snap',
                                  'redshift',
                                  'instsfr',
                                  'sfr_50',
                                  'M*',
                                  'Mgas',
                                  'MH2',
                                  'MHI',
                                  'HMR',
                                  'FMR',
                                  'metallicity',
                                  'xpos',
                                  'ypos',
                                  'zpos',
                                  'mhalo',
                                  'central_stellar_mass',
                                  'metal_radius_list',
                                  'cgf1',
                                  'cgf2'])
        
    
        if HALOS == False:
            outfile = BASEDIR+SNAPPRE+'/output/Groups/caesar_physical_properties.local.'+str(lh)+'.dat'
            outsavefile = BASEDIR+SNAPPRE+'/output/Groups/caesar_physical_properties.local.'+str(lh)+'.npz'
        else:
            outfile = BASEDIR+SNAPPRE+'/output/Groups/caesar_physical_properties.halos.local.'+str(lh)+'.dat'
            outsavefile = BASEDIR+SNAPPRE+'/output/Groups/caesar_physical_properties.halos.local.'+str(lh)+'.npz'
            

        ascii.write(phys_table,outfile)
        np.savez(outsavefile,snapshotnames=snapshotnames,
                 snaps=snaps,
                 redshifts=redshifts,
                 instsfr=instsfr,
                 sfr_50=sfr_50,
                 stellar_masses=stellar_masses,
                 gas_masses=gas_masses,
                 h2_masses=h2_masses,
                 hI_masses=hI_masses,
                 halo_masses=halo_masses,
                 hmr=hmr,
                 fmr=fmr,
                 metallicity=metallicity,
                 xpos=xpos,
                 ypos=ypos,
                 zpos=zpos,
                 central_stellar_masses = central_stellar_masses,
                 metal_radius_list = metal_radius_list,
                 r_baryon = r_baryon,
                 r_h_baryon = r_h_baryon,
                 r_gas = r_gas,
                 r_h_gas = r_h_gas,
                 r200 = r200,
                 r200c = r200c,
                 rstellar = rstellar,
                 r_h_stellar = r_h_stellar,
                 r_virial = r_virial,
                 stellar_ages = stellar_ages,
                 cgf1 = cgf1,
                 cgf2 = cgf2)
        
