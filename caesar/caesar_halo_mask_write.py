import numpy as np
import caesar,yt
import sys
sys.path.insert(0,'inputfiles/')
from caesar_halo_mask_write_input import *

obj = caesar.load(caesarfile)
ic = icfile
ds = yt.load(snapshot)
ic_ds = yt.load(ic)
obj.yt_dataset = ds
obj.halos[halonum].write_IC_mask(ic_ds,outfile)
