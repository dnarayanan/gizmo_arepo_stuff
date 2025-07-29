import numpy as np
import caesar
import yt
import multiprocessing
from functools import partial

# ===============================================================
#                  CONFIGURATION PARAMETERS
# ===============================================================
# Set to a number (e.g., 1000) for a quick test, or None for the full run.
LIMIT = None

# --- File Paths
SNAPSHOT_FILE = '/blue/narayanan/desika.narayanan/arepo_runs/cosmic_sands/m100/dm_boxes/m100n256/output/snapdir_019/snapshot_019.0.hdf5'
CAESAR_FILE = '/blue/narayanan/desika.narayanan/arepo_runs/cosmic_sands/m100/dm_boxes/m100n256/output/Groups/caesar_snapshot_019.hdf5'
OUTPUT_FILE = 'isolated_halos_test.txt' # Using a different name for the test file

# --- Selection Criteria
MASS_RANGE = [5e11, 5e12]
ISOLATION_MASS_RATIO = 10.0
ISOLATION_RADIUS_MULTIPLE = 3.0
N_CORES = multiprocessing.cpu_count()
# ===============================================================

def calculate_r200c_from_value(mass_msun, rho_crit_msun_kpccm3):
    r_200c_cubed = (3 * mass_msun) / (4 * np.pi * 200 * rho_crit_msun_kpccm3)
    return r_200c_cubed**(1/3)

def periodic_distance(pos1, pos2, box_size):
    diff = np.abs(pos1 - pos2)
    diff[diff > box_size / 2.0] = box_size - diff[diff > box_size / 2.0]
    return np.linalg.norm(diff)

# --- Worker for Stage 1: Find candidates
def check_halo_isolation(central_halo_data, **kwargs):
    central_mass = central_halo_data['mass']
    if not (kwargs['mass_range'][0] <= central_mass <= kwargs['mass_range'][1]):
        return None

    r_200c = calculate_r200c_from_value(central_mass, kwargs['rho_crit_value'])
    check_radius = kwargs['iso_radius_multiple'] * r_200c
    companion_mass_limit = central_mass / kwargs['iso_mass_ratio']

    for other_halo_data in kwargs['all_halos_data']:
        if central_halo_data['index'] == other_halo_data['index']: continue
        if other_halo_data['mass'] < companion_mass_limit: break
        
        distance = periodic_distance(central_halo_data['pos'], other_halo_data['pos'], kwargs['box_size_kpccm'])
        if distance < check_radius:
            return None
    return central_halo_data['index']

# --- Worker for Stage 2: Calculate priority metrics
def calculate_priority_metric(central_halo_data, **kwargs):
    min_dist_to_massive_neighbor = float('inf')
    companion_mass_limit = central_halo_data['mass'] / kwargs['iso_mass_ratio']
    r_200c = calculate_r200c_from_value(central_halo_data['mass'], kwargs['rho_crit_value'])

    for other_halo_data in kwargs['all_halos_data']:
        if central_halo_data['index'] == other_halo_data['index']: continue
        if other_halo_data['mass'] < companion_mass_limit: break
        
        dist = periodic_distance(central_halo_data['pos'], other_halo_data['pos'], kwargs['box_size_kpccm'])
        if dist < min_dist_to_massive_neighbor:
            min_dist_to_massive_neighbor = dist
            
    isolation_metric = min_dist_to_massive_neighbor / r_200c
    return {'index': central_halo_data['index'], 'metric': isolation_metric}

if __name__ == '__main__':
    print("Loading data...")
    ds = yt.load(SNAPSHOT_FILE)
    obj = caesar.load(CAESAR_FILE, ds)

    print("Extracting halo data into a pickle-safe format...")
    halo_data_list = [{'index': h._index, 'mass': h.masses['total'].in_units('Msun').value,
                       'pos': h.pos.in_units('kpccm').value} for h in obj.halos]
    halo_data_list.sort(key=lambda h: h['mass'], reverse=True)
    
    # Apply the limit for testing
    if LIMIT is not None:
        halos_to_process = halo_data_list[:LIMIT]
        print(f"--- RUNNING IN TEST MODE: PROCESSING THE {LIMIT} MOST MASSIVE HALOS ---")
    else:
        halos_to_process = halo_data_list

    worker_kwargs = {
        'all_halos_data': halo_data_list, # Workers always check against the FULL list
        'mass_range': MASS_RANGE,
        'iso_mass_ratio': ISOLATION_MASS_RATIO,
        'iso_radius_multiple': ISOLATION_RADIUS_MULTIPLE,
        'box_size_kpccm': ds.domain_width[0].in_units('kpccm').value,
        'rho_crit_value': ds.critical_density.in_units('Msun/kpccm**3').value
    }

    print(f"Stage 1: Finding isolated candidates from a sample of {len(halos_to_process)} halos using {N_CORES} cores...")
    stage1_worker = partial(check_halo_isolation, **worker_kwargs)
    with multiprocessing.Pool(processes=N_CORES) as pool:
        isolated_indices = [idx for idx in pool.map(stage1_worker, halos_to_process) if idx is not None]

    print(f"Found {len(isolated_indices)} candidates. Now calculating priority metrics.")
    isolated_candidates_data = [h for h in halo_data_list if h['index'] in isolated_indices]
    
    print(f"Stage 2: Calculating metrics in parallel using {N_CORES} cores...")
    stage2_worker = partial(calculate_priority_metric, **worker_kwargs)
    with multiprocessing.Pool(processes=N_CORES) as pool:
        final_candidates = pool.map(stage2_worker, isolated_candidates_data)

    final_candidates.sort(key=lambda c: c['metric'], reverse=True)
    
    header_title = f"Found {len(final_candidates)} halos meeting the fiducial criteria:"
    header_params = f"(Mass Range = [{MASS_RANGE[0]:.1e}, {MASS_RANGE[1]:.1e}] Msun, Isolation Ratio > {ISOLATION_MASS_RATIO:.0f}:1)"
    column_names = f"{'Index':<8}{'Mass (Msun)':<18}{'Isolation Metric':<20}{'Spin Param':<15}{'Position (cMpc/h)':<30}"
    separator = f"{'-'*7:<8}{'-'*16:<18}{'-'*18:<20}{'-'*13:<15}{'-'*28:<30}"

    with open(OUTPUT_FILE, 'w') as f:
        print("\n" + "="*80)
        f.write(f"# {header_title}\n# {header_params}\n")
        print(header_title); print(header_params); print("="*80)
        
        print(column_names); f.write(f"#{column_names}\n")
        print(separator); f.write(f"#{separator}\n")

        for cand in final_candidates:
            halo = obj.halos[cand['index']]
            mass_str = f"{halo.masses['total'].in_units('Msun').value:.3e}"
            pos_str = np.array_str(halo.pos.in_units('Mpccm').value, precision=2)
            metric_str = f"{cand['metric']:.2f}"
            spin_str = f"{halo.spin_parameter:.3f}" if hasattr(halo, 'spin_parameter') and halo.spin_parameter is not None else "N/A"
            line_out = f"{halo._index:<8}{mass_str:<18}{metric_str:<20}{spin_str:<15}{pos_str:<30}"
            print(line_out)
            f.write(line_out + "\n")

    print("\n" + "="*80)
    print(f"✅ Results sorted by Isolation Metric and saved to {OUTPUT_FILE}")
