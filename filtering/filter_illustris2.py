import h5py
import argparse
import sys,os
sys.path.insert(0,'paul_illustris_tools')
import paul_illustris_tools.readhaloHDF5 as readhalo
import numpy as np
from glob import glob

base_path = '/blue/narayanan/desika.narayanan/arepo_runs/m25n512b_09c6391_nolimits_filter_ready/output/'
snap_num = 16

DEBUG = False

print('halo reader')
h = readhalo.HaloReader(str(base_path), '0'+str(snap_num), int(snap_num))

print('getting offsets')


num_galaxies = h.halo_offset.shape[0]
num_galaxies = 1000


#first make sure the output directory is in place
output_dir = base_path+'/filtered_snaps/snap'+str(snap_num).zfill(3)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(output_dir+" does not exist: creating now")

for gal in range(num_galaxies):

    print("working on galaxy: ",gal)
    outfile = output_dir+'/galaxy_'+str(gal)+'.hdf5'
    print("outfile: ",outfile)
    
    gas_offset = h.halo_offset[gal][0]
    dust_offset = h.halo_offset[gal][3]
    star_offset = h.halo_offset[gal][4]
    gas_len = h.cat.SubhaloLenType[gal][0]
    star_len = h.cat.SubhaloLenType[gal][4]
    dust_len = h.cat.SubhaloLenType[gal][3]
    
    if DEBUG:
        print('star offset:',star_offset)
        print('star len:',star_len)
        print('gas offset:',gas_offset)
        print('gas len:',gas_len)
        print('dust offset:',dust_offset)
        print('dust len:',dust_len)

    
    files_dir = base_path+'/snapdir_'+'0'+str(snap_num)+'/'
    files = sorted(glob(files_dir+'*.hdf5'), key = lambda name: int(name.split('.')[-2]))
    this_file_start0 = 0
    this_file_start3 = 0
    this_file_start4 = 0
    

    if DEBUG: print('starting with gas particles')
    for snap_file in files:
        with h5py.File(snap_file, 'r') as input_file:
            
            print(snap_file)
            
            this_file_end0 = this_file_start0 + len(input_file['PartType0']['Masses'])
            
            #gas_checks
            if (gas_offset+gas_len) > this_file_end0:
                this_file_start0 += len(input_file['PartType0']['Masses'])
                continue
            
            elif this_file_start0 > (gas_offset+gas_len):
                
                break
            
            else:
                output_file = h5py.File(outfile,"a")
                output_file.copy(input_file['Header'], 'Header')
                output_file.copy(input_file['Config'], 'Config')
                output_file.create_group('PartType0')
                
                
                if DEBUG: print('writing from file',snap_file)
            
                for k in input_file['PartType0']:
                    if k in output_file['PartType0']:
                        if DEBUG: print('[gas] we are writing for the second time, i.e. galaxy data exists across two snapshots')
                        in_data = input_file['PartType0'][k][0 : leftover]
                        output_file['PartType0'][k] = np.hstack(output_file['PartType0'][k], in_data)
                        if DEBUG: print(len(in_data))
                    else: 
                        if DEBUG: print('[gas] writing for the first time')
                        in_data = input_file['PartType0'][k][gas_offset-this_file_start0 : min(this_file_end0,(gas_offset+gas_len)-this_file_start0)]
                        output_file['PartType0'][k] = in_data
                        leftover = (gas_offset+gas_len) - this_file_end0
                this_file_start0 += len(input_file['PartType0']['Masses'])
            
                output_file.close()



    if DEBUG: print('now doing dust particles')
    for snap_file in files:
        with h5py.File(snap_file, 'r') as input_file:
            
            this_file_end3 = this_file_start3 + len(input_file['PartType3']['Masses'])
            
            #dust_checks
            if (dust_offset+dust_len) > this_file_end3:
                if DEBUG: print('[dust] got to first check')
                if DEBUG: print('[dust] gal:', (dust_offset+dust_len))
                this_file_start3 += len(input_file['PartType3']['Masses'])
                continue
            
            elif this_file_start3 > (dust_offset+dust_len):
                if DEBUG: print('[dust] breaking due to second check')
                break
            
            else:
                output_file = h5py.File(outfile,"a")
                output_file.create_group('PartType3')
                
                for k in input_file['PartType3']:
                    if k in output_file['PartType3']:
                        if DEBUG: print('[dust] we are writing for the second time, i.e. galaxy data exists across two files')
                        in_data = input_file['PartType3'][k][0 : leftover]
                        output_file['PartType3'][k] = np.hstack(output_file['PartType3'][k], in_data)
                        
                    else:
                        if DEBUG: print('[dust] writing for the first time')
                        in_data = input_file['PartType3'][k][dust_offset-this_file_start3 : min(this_file_end3,(dust_offset+dust_len)-this_file_start3)]
                        output_file['PartType3'][k] = in_data
                        
                        leftover = (dust_offset+dust_len) - this_file_end3
                this_file_start3 += len(input_file['PartType3']['Masses'])
                output_file.close()
                        

                
    if DEBUG: print('now doing star particles')
    for snap_file in files:
        with h5py.File(snap_file, 'r') as input_file:
            
            this_file_end4 = this_file_start4 + len(input_file['PartType4']['Masses'])
            
            #star_checks                                                                           
            if (star_offset+star_len) > this_file_end4:
                if DEBUG: print('[stars] got to first check')
                if DEBUG: print('[stars] gal:', (star_offset+star_len))
                this_file_start4 += len(input_file['PartType4']['Masses'])
                continue
            
            elif this_file_start4 > (star_offset+star_len):
                if DEBUG: print('[stars] breaking due to second check')
                break
            
            else: 
                output_file = h5py.File(outfile,"a")
                output_file.create_group('PartType4')
                
                for k in input_file['PartType4']:
                    if k in output_file['PartType4']:
                        if DEBUG: print('[stars] we are writing for the second time, i.e. galaxy data exists across two files')
                        in_data = input_file['PartType4'][k][0 : leftover]
                        output_file['PartType4'][k] = np.hstack(output_file['PartType4'][k], in_data)
                        
                    else:
                        if DEBUG: print('[stars] writing for the first time')
                        in_data = input_file['PartType4'][k][star_offset-this_file_start4 : min(this_file_end4,(star_offset+star_len)-this_file_start4)]
                        output_file['PartType4'][k] = in_data
                        leftover = (star_offset+star_len) - this_file_end4
                this_file_start4 += len(input_file['PartType4']['Masses'])
                output_file.close()

    re_out = h5py.File(outfile,"a")
    
    re_out['Header'].attrs.modify('NumFilesPerSnapshot', 1)
    re_out['Header'].attrs.modify('NumPart_ThisFile', np.array([gas_len, 0, 0, dust_len, star_len, 0]))
    re_out['Header'].attrs.modify('NumPart_Total', np.array([gas_len, 0, 0, dust_len, star_len, 0]))
    
    re_out.close()

    
