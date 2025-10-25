#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import numpy as np

def parse_cpu_file(filepath):
    """
    Parses a single Arepo cpu.txt file to extract simulation time
    and cumulative wall clock time.
    
    This version dynamically handles both the old (5-column) and
    new (6-column) 'total' line formats.
    """
    if not os.path.exists(filepath):
        print(f"Error: File not found at '{filepath}'")
        return []

    data_points = []
    print(f"Reading {os.path.basename(filepath)}...")
    try:
        with open(filepath, 'r') as f:
            current_sim_time = None
            for line in f:
                # Look for the line containing the simulation time
                if line.startswith('Step') and 'Time:' in line:
                    try:
                        time_str = line.split('Time:')[1].split(',')[0].strip()
                        current_sim_time = float(time_str)
                    except (ValueError, IndexError):
                        current_sim_time = None
                        
                # Look for the line with total timings for the step
                elif line.strip().startswith('total') and current_sim_time is not None:
                    try:
                        parts = line.split()
                        
                        # --- MODIFIED LOGIC TO HANDLE BOTH FORMATS ---
                        cumulative_wall_time = None
                        
                        # Original format (Arepo1): 5 columns
                        # total, diff, diff%, cumulative, cumulative%
                        if len(parts) == 5:
                            cumulative_wall_time = float(parts[3])
                            
                        # New format (Arepo2): 6 columns
                        # total, diff, diff%, max, cumulative, cumulative%
                        elif len(parts) == 6:
                            cumulative_wall_time = float(parts[4])
                        
                        else:
                            # Unknown format, skip this line
                            print(f"Warning: Unrecognized 'total' line format in {filepath}: {line.strip()}")
                            continue
                        # --- END OF MODIFIED LOGIC ---

                        data_points.append((current_sim_time, cumulative_wall_time))
                        current_sim_time = None # Reset after use
                        
                    except (ValueError, IndexError, TypeError) as e:
                        print(f"Warning: Could not parse line: {line.strip()} - {e}")
                        continue
                        
    except IOError as e:
        print(f"Error reading file {filepath}: {e}")

    return data_points

def plot_direct_comparison(all_data):
    """
    Generates the first plot, directly comparing the wall clock times.
    """
    print("\n--- Generating Plot 1: Direct Performance Comparison ---")
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 7))

    for label, data_points in all_data.items():
        if not data_points:
            continue
        
        data_points.sort()
        simulation_times, wall_times = zip(*data_points)
        wall_times_hours = [t / 3600.0 for t in wall_times]
        ax.plot(simulation_times, wall_times_hours, marker='.', markersize=4, linestyle='-', label=label)

    ax.set_xlabel('Scale Factor (a)', fontsize=12)
    ax.set_ylabel('Cumulative Wall Clock Time (Hours)', fontsize=12)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=12)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    output_filename = "comparison_timing_plot.png"
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"✅ Plot 1 saved to '{output_filename}'")
    plt.show()

def plot_speedup_ratio(all_data, baseline_label, comparison_label):
    """
    Generates the second plot, showing the speed-up ratio.
    """
    print(f"\n--- Generating Plot 2: Speed-up Ratio ---")
    print(f"Baseline run: '{baseline_label}'")
    print(f"Comparison run: '{comparison_label}'")

    if baseline_label not in all_data or comparison_label not in all_data:
        print("Error: Labels for speed-up calculation not found in the processed data.")
        return

    # --- Prepare data and handle potential emptiness ---
    base_data = all_data[baseline_label]
    comp_data = all_data[comparison_label]
    if not base_data or not comp_data:
        print("Error: Cannot calculate ratio because one of the runs has no data.")
        return
        
    base_data.sort()
    comp_data.sort()
    
    base_sim_times, base_wall_times = zip(*base_data)
    comp_sim_times, comp_wall_times = zip(*comp_data)

    # --- Interpolation ---
    comp_wall_times_interp = np.interp(base_sim_times, comp_sim_times, comp_wall_times)

    # --- Calculate the Ratio ---
    with np.errstate(divide='ignore', invalid='ignore'):
        # Ratio = Time_comparison / Time_baseline
        speedup_ratio = np.divide(comp_wall_times_interp, base_wall_times)

    # Clean up non-finite values (inf, nan) that can result from division by zero at the start
    final_sim_times = np.array(base_sim_times)
    final_ratio = np.array(speedup_ratio)
    finite_mask = np.isfinite(final_ratio)
    
    # --- Create the Plot ---
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 7))

    ax.plot(final_sim_times[finite_mask], final_ratio[finite_mask], marker='o', markersize=3, linestyle='-', color='C2')
    
    ax.axhline(1.0, color='k', linestyle='--', linewidth=1.5, label='No Speed-up (Ratio = 1)')

    ax.set_xlabel('Scale Factor (a)', fontsize=12)
    ax.set_ylabel(f"Speed-up Factor ({comparison_label} Time / {baseline_label} Time)", fontsize=12)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=12)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    output_filename = "speedup_ratio_plot.png"
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"✅ Plot 2 saved to '{output_filename}'")
    plt.show()

def main():
    """
    Main function to define files, parse them, and generate comparison plots.
    """
    # --- ⚙️ EDIT THIS SECTION ---
    # List of the full paths to your cpu.txt files
    file_paths = [
        '/blue/narayanan/desika.narayanan/arepo_runs/cosmic_sands/m100/z0/arepo2_test/output/cpu.txt',
        '/blue/narayanan/desika.narayanan/arepo_runs/cosmic_sands/m100/z0/sobol_044_production/gal9265_sobol_044/output/cpu.txt'

    ]

    # Corresponding labels for the plot legend
    labels = ['A2', 'A1']

    # Define which labels to use for the speed-up calculation.
    # To see how much faster 'A2' is, the ratio should be Time(A1) / Time(A2)
    baseline_run_label = 'A2'
    comparison_run_label = 'A1'
    # -----------------------------

    if len(file_paths) != len(labels):
        print("Error: The number of file paths must match the number of labels.")
        return

    # --- Parse all files and store data in a dictionary ---
    all_data = {}
    for filepath, label in zip(file_paths, labels):
        all_data[label] = parse_cpu_file(filepath)

    # --- Generate Plots ---
    plot_direct_comparison(all_data)
    plot_speedup_ratio(all_data, baseline_label=baseline_run_label, comparison_label=comparison_run_label)


if __name__ == "__main__":
    # Ensure you have numpy and matplotlib installed
    # pip install numpy matplotlib
    main()
