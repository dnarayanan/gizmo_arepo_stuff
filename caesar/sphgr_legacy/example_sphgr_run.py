import yt
import yt.analysis_modules.sphgr.api as sphgr
import ipdb
import numpy as np

snapshot = '/Users/desika/Dropbox/zooms/bobby_dm/snap_N128L16_085.hdf5'

ds = yt.load(snapshot)

#this is a class that has attributes like galaxies, halos, and the yt
#dataset, as well as the standard sphgr methods. this is like
#importing the old picklefiles
obj = sphgr.SPHGR() 
obj.yt_dataset = ds  #note, alternatively i could have just said obj =
                     #sphgr.SPHGR(ds)

obj.member_search() #obviously this does member search

#now we need to write the output to a memberfile (except it'll be hdf5
#now instead of pickle)

obj.save_sphgr_data("halos.hdf5")

#example analysis

h = obj.halos[0]  #this will grab the most massive halo
#now do a dir(h)

#if you were to want a galaxy then you would do something like:

#g = obj.halos.mass['dm']

