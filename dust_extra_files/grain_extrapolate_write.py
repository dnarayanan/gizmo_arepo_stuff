import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import interp
from scipy.optimize import curve_fit


import numpy as np
import h5py
from astropy.modeling import models
from scipy.signal import savgol_filter
from astropy.convolution import convolve,Box1DKernel

def fx(x,a,b):
    return a*x+b

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()

    return idx

magic_number=33 #this is the index in mckinnon's grain size distribution where we approximate the powerlaw break at - roughly 1.25 cm

#load the mckinnon GSD
hf  = h5py.File('GrainData.hdf5','r')

radius = hf['InitialSNII']['Radius'].value
GSD = hf['InitialSNII']['GSD'].value

#fill in the values that clearly look like a numerical problem
GSD_fill_in = interp(GSD[1:5],radius,GSD)
GSD[1:5] = GSD_fill_in

#fit the data into two powerlaws
popt1,pcov1 = curve_fit(fx,np.log10(radius[0:magic_number]),np.log10(GSD[0:magic_number]))
popt2,pcov2 = curve_fit(fx,np.log10(radius[magic_number:-1]),np.log10(GSD[magic_number:-1]))

#read in the grain sizes we want to use from Bruce's work
draine_sizes = np.loadtxt('draine_sizes.txt')

#find where the maximum draine size hits the mckinnon size array, and
#then wipe out all the stuff in mckinnon's array before that, and then
#merge them.
mask_idx = find_nearest(radius,np.max(draine_sizes))
if radius[mask_idx] < np.max(draine_sizes): mask_idx +=1 #in case the nearest index was just one less than the max value
final_compiled_radii = np.asarray(list(draine_sizes)+list(radius[mask_idx:-1]))

#now make GSD:
idx_to_break_at = find_nearest(final_compiled_radii,radius[magic_number])
gsd1 = 10.**(fx(np.log10(final_compiled_radii[0:idx_to_break_at]),popt1[0],popt1[1]))
gsd2 = 10.**(fx(np.log10(final_compiled_radii[idx_to_break_at::]),popt2[0],popt2[1]))
final_compiled_gsd = np.asarray(list(gsd1)+list(gsd2))

##make a new GSD array thats twice as long as original GSD.  The values that we prepend on will just look like GSD[0]
#GSD2 = np.zeros(len(GSD)*2)
#GSD2[0:len(radius)] = GSD[0]
#GSD2[len(radius)::] = GSD

##make new radiuses that go down to 1e-10 cm
#radius2 = np.zeros(len(radius)*2)
#radius2[0:len(radius)] = 10.**(np.linspace(-10,-4,len(radius)))
#radius2[len(radius)::] = radius

hf.close()

#write
outfile = 'GrainData_extrap.hdf5'

with h5py.File('GrainData.hdf5','r') as input_file, h5py.File(outfile,'w') as output_file:
#    output_file.create_group('CondEffAGB')
    output_file.copy(input_file['CondEffAGB'],'CondEffAGB')
    output_file.copy(input_file['CondEffSNII'],'CondEffSNII')
    output_file.copy(input_file['SNDestruction'],'SNDestruction')
    output_file.copy(input_file['VelocitiesCNM'],'VelocitiesCNM')
    output_file.copy(input_file['VelocitiesDC1'],'VelocitiesDC1')
    output_file.copy(input_file['VelocitiesDC2'],'VelocitiesDC2')
    output_file.copy(input_file['VelocitiesMC'],'VelocitiesMC')
    output_file.copy(input_file['VelocitiesWIM'],'VelocitiesWIM')
    output_file.copy(input_file['VelocitiesWNM'],'VelocitiesWNM')

    #now get the Initial SNII data
    #output_file.copy(input_file['InitialSNII'],'InitialSNII')
    output_file.create_group('InitialSNII')
    output_file['InitialSNII']['NPoints'] = len(final_compiled_radii)
    output_file['InitialSNII']['Radius'] = final_compiled_radii
    output_file['InitialSNII']['GSD'] = final_compiled_gsd

fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(radius[0:40],GSD[0:40])
ax.plot(radius[0:magic_number],10.**(fx(np.log10(radius[0:magic_number]),popt1[0],popt1[1])))
ax.plot(radius[magic_number:-1],10.**(fx(np.log10(radius[magic_number:-1]),popt2[0],popt2[1])))

#ax.loglog(radius,10.**fx(np.log10(np.log10(radius)),1,3,3))
#ax.set_ylim([1.e-5,1.e22])
fig.savefig('junk.png',dpi=300)
