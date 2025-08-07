import h5py
import numpy as np
from scipy.optimize import curve_fit
import os

# ===================================================================
# --- 1. Configuration ---
# ===================================================================
# Path to your complete, original, "fast" HDF5 file
SOURCE_FILE = '/blue/narayanan/desika.narayanan/arepo_runs/cosmic_sands/m100/z0/pzr_radical_sobol/gal9265_ml12_areverse_extendedgsd/GrainData_extrap.hdf5' 

# Set this to the new, smaller value you want to test (e.g., 0.001)
MIN_GRAIN_SIZE_UM = 0.0004

# ===================================================================
# --- Automatically generated filenames and parameters ---
# ===================================================================
OUTPUT_FILE = f'GrainData_mgs{MIN_GRAIN_SIZE_UM:.4f}_cloned.hdf5'
original_min_size_um = 0.0004
SCALE_FACTOR = MIN_GRAIN_SIZE_UM / original_min_size_um

NUM_NEW_POINTS = 40
BREAK_RADIUS_CM = 1.25e-4 

# ===================================================================
# --- Helper Functions to Calculate New Data ---
# ===================================================================
# (These are the same as before, they just calculate the new numpy arrays)
def get_new_gsd_data(hf_in):
    original_radii = hf_in['InitialSNII/Radius'][:]
    original_gsd = hf_in['InitialSNII/GSD'][:]
    break_idx = np.argmin(np.abs(original_radii - BREAK_RADIUS_CM))
    log_r = np.log10(original_radii)
    log_gsd = np.log10(original_gsd)
    popt, _ = curve_fit(lambda x, a, b: a*x + b, log_r[:break_idx], log_gsd[:break_idx])
    slope, intercept = popt[0], popt[1]
    new_min_log_r = np.log10(np.min(original_radii) * SCALE_FACTOR)
    original_min_log_r = np.min(log_r)
    prepended_log_r = np.linspace(new_min_log_r, original_min_log_r, NUM_NEW_POINTS, endpoint=False)
    final_log_r = np.concatenate([prepended_log_r, log_r])
    final_radii = 10.**np.unique(final_log_r)
    final_log_r = np.log10(final_radii)
    new_break_idx = np.argmin(np.abs(final_radii - BREAK_RADIUS_CM))
    log_gsd_extrap = slope * final_log_r[:new_break_idx] + intercept
    interp_log_gsd = np.interp(final_log_r[new_break_idx:], log_r, log_gsd)
    final_log_gsd = np.concatenate([log_gsd_extrap, interp_log_gsd])
    final_gsd = 10.**final_log_gsd
    return final_radii, final_gsd

def get_new_1d_table_data(hf_in, group_name, value_name):
    original_radii = hf_in[f'{group_name}/Radius'][:]
    original_values = hf_in[f'{group_name}/{value_name}'][:]
    new_min_log_r = np.log10(np.min(original_radii) * SCALE_FACTOR)
    original_min_log_r = np.log10(np.min(original_radii))
    prepended_log_r = np.linspace(new_min_log_r, original_min_log_r, NUM_NEW_POINTS, endpoint=False)
    final_log_r = np.concatenate([prepended_log_r, np.log10(original_radii)])
    final_radii = 10.**np.unique(final_log_r)
    final_values = np.interp(final_radii, original_radii, original_values, left=original_values[0], right=original_values[-1])
    return final_radii, final_values

# ===================================================================
# --- The "Perfect Cloner" Core Logic ---
# ===================================================================
def clone_hdf5_object(source_group, dest_group, modified_data):
    """
    Recursively clones groups and datasets, preserving properties.
    Substitutes data from the modified_data dictionary where specified.
    """
    for name, obj in source_group.items():
        full_path = obj.name
        
        if isinstance(obj, h5py.Dataset):
            # Check if this is a dataset we have pre-calculated new data for
            if full_path in modified_data:
                data = modified_data[full_path]
                # Use original properties but new data
                dest_group.create_dataset(name, data=data, 
                                          compression=obj.compression, 
                                          chunks=obj.chunks)
            else:
                # Perform a direct copy, preserving all properties
                source_group.copy(name, dest_group)
        
        elif isinstance(obj, h5py.Group):
            # Create the new group
            new_group = dest_group.create_group(name)
            # Recurse into the group
            clone_hdf5_object(obj, new_group, modified_data)

# ===================================================================
# --- Main Execution ---
# ===================================================================
if __name__ == "__main__":
    if not os.path.exists(SOURCE_FILE):
        raise FileNotFoundError(f"Source file '{SOURCE_FILE}' not found.")

    # --- Step 1: Pre-calculate all new data arrays in memory ---
    print("--- Calculating all modified data arrays in memory ---")
    modified_data = {}
    with h5py.File(SOURCE_FILE, 'r') as hf_in:
        new_radii_snii, new_gsd = get_new_gsd_data(hf_in)
        modified_data['/InitialSNII/Radius'] = new_radii_snii
        modified_data['/InitialSNII/GSD'] = new_gsd
        modified_data['/InitialSNII/NPoints'] = np.int64(len(new_radii_snii))

        tables_to_modify = {
            'VelocitiesCNM': 'Velocity', 'VelocitiesDC1': 'Velocity', 'VelocitiesDC2': 'Velocity',
            'VelocitiesMC': 'Velocity', 'VelocitiesWIM': 'Velocity', 'VelocitiesWNM': 'Velocity',
            'CoulombCNM_Gra': 'Coulomb', 'CoulombWNM_Gra': 'Coulomb', 'CoulombWIM_Gra': 'Coulomb',
            'CoulombCNM_Sil': 'Coulomb', 'CoulombWNM_Sil': 'Coulomb', 'CoulombWIM_Sil': 'Coulomb'
        }
        for group, value_name in tables_to_modify.items():
            if group in hf_in:
                new_radii, new_values = get_new_1d_table_data(hf_in, group, value_name)
                modified_data[f'/{group}/Radius'] = new_radii
                modified_data[f'/{group}/{value_name}'] = new_values
                if f'/{group}/NPoints' in hf_in:
                    modified_data[f'/{group}/NPoints'] = np.int64(len(new_radii))
        
        if 'SNDestruction' in hf_in:
            original_radii_sd = hf_in['/SNDestruction/Radius'][:]
            original_xifrac = hf_in['/SNDestruction/XiFrac'][:]
            new_radii_sd, _ = get_new_1d_table_data(hf_in, 'SNDestruction', 'Radius')
            
            n_new = len(new_radii_sd)
            final_xifrac = np.zeros((n_new, n_new))
            for i_new, r_new in enumerate(new_radii_sd):
                i_old = np.argmin(np.abs(original_radii_sd - r_new))
                for j_new, r2_new in enumerate(new_radii_sd):
                    j_old = np.argmin(np.abs(original_radii_sd - r2_new))
                    final_xifrac[i_new, j_new] = original_xifrac[i_old, j_old]
            
            modified_data['/SNDestruction/Radius'] = new_radii_sd
            modified_data['/SNDestruction/XiFrac'] = final_xifrac
            if '/SNDestruction/NRadius' in hf_in:
                 modified_data['/SNDestruction/NRadius'] = np.int64(n_new)
    
    # --- Step 2: Create the new file by cloning the original ---
    print(f"\n--- Cloning '{SOURCE_FILE}' to '{OUTPUT_FILE}' with substitutions ---")
    with h5py.File(SOURCE_FILE, 'r') as hf_in, h5py.File(OUTPUT_FILE, 'w') as hf_out:
        clone_hdf5_object(hf_in, hf_out, modified_data)

    print(f"\n✅ File '{OUTPUT_FILE}' created successfully.")
    print("Please use this file and the corresponding MinGrainSize in your parameter file.")
