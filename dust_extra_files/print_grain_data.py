import h5py
import numpy as np

def print_hdf5_structure(name, obj):
    """Recursively prints the structure of an HDF5 file."""
    print(name)
    # If it's a dataset, print its shape and dtype
    if isinstance(obj, h5py.Dataset):
        print("  - Shape:", obj.shape)
        print("  - Dtype:", obj.dtype)
    # Print attributes for both groups and datasets
    for key, val in obj.attrs.items():
        print(f"  - Attr: {key} = {val}")

# --- Main Script ---
filename = 'GrainData_extrap.hdf5'

try:
    with h5py.File(filename, 'r') as f:
        print(f"--- Structure of {filename} ---")
        f.visititems(print_hdf5_structure)
        print("---------------------------------")
except FileNotFoundError:
    print(f"Error: The file '{filename}' was not found in the current directory.")
except Exception as e:
    print(f"An error occurred: {e}")
