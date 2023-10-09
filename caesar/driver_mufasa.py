import numpy as np
import yt
import caesar
import sys
sys.path.insert(0,'inputfiles/')
from driver_input_mufasa import *

caesar.drive(snapdir,snapname,snapnums,progen=True)
