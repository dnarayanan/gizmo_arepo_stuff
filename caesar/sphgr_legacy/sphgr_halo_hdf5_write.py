#function call: 

#>halo_hdf5_write_step1 snapshot_name output_hdf5_file_name

#to load halos to identify for higher-res resimulation in zooms

import yt
import yt.analysis_modules.sphgr.api as sphgr
import ipdb
import numpy as np
from sys import argv

scriptname,snapshot,fname_out = argv

ds = yt.load(snapshot)

#save the halos -- note you only have to do this once
obj = sphgr.SPHGR(ds)
obj.member_search()
obj.save_sphgr_data(fname_out)






