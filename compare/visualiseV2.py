import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import os

# --- Configuration ---
JSON_FILE_PATH = "output/results.json"  # Adjust if your file is elsewhere
OUTPUT_DIR = "output/plots"
REFERENCE_SUBDIVISION_FACTOR = 0
# --- End Configuration ---


# --- Helper Functions ---
def _save_and_close_plot(figure_path):
    """Saves the current matplotlib figure and closes it."""
    plt.savefig(figure_path)
    plt.close()


def _setup_common_line_plot_aesthetics(title, xlabel, ylabel, x_ticks_values, add_margins=True):
    """Sets common aesthetic elements for line plots."""
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.xticks(x_ticks_values)
    if add_margins:
        plt.margins(x=0.05)


# --- Individual Plotting Functions ---


def plot_ram_usage_for_system(system_name, system_df, subdivision_factors_num, output_dir):
    """Plots RAM usage vs. subdivision factor."""
    print(f"  Generating RAM plot for {system_name}...")
    plt.figure(figsize=(10, 6))
    plt.plot(subdivision_factors_num, system_df["max_ram_mb"], marker="o", linestyle="-", color="blue")
    _setup_common_line_plot_aesthetics(
        title=f"System '{system_name}': RAM Usage vs. Subdivision Factor", xlabel="Subdivision Factor", ylabel="Max RAM Usage (MB)", x_ticks_values=subdivision_factors_num
    )
    plt.tight_layout()
    _save_and_close_plot(os.path.join(output_dir, f"{system_name}_ram_usage.png"))


def plot_runtime_breakdown_for_system(system_name, system_df, subdivision_factors_num, subdivision_factors_str, output_dir):
    """Plots stacked bar chart of runtime breakdown."""
    print(f"  Generating Runtime plot for {system_name}...")
    plt.figure(figsize=(10, 6))

    t_init_duration = system_df["time_init"].values
    # Assuming time_molorb is cumulative walltime after molorb step finishes
    t_molorb_cumulative = system_df["time_molorb"].values
    t_wall_duration = system_df["time_wall"].values

    # Durations of each component for the stacked bar
    init_component = t_init_duration
    molorb_calc_component = np.maximum(0, t_molorb_cumulative - t_init_duration)
    overhead_component = np.maximum(0, t_wall_duration - t_molorb_cumulative)

    bar_width = 0.6
    x_pos = np.arange(len(subdivision_factors_num))  # Numeric positions for bars

    plt.bar(x_pos, init_component, label="Initialization (time_init)", color="skyblue", width=bar_width)
    plt.bar(x_pos, molorb_calc_component, label="MolOrb Calculation (time_molorb - time_init)", bottom=init_component, color="lightcoral", width=bar_width)
    # The bottom for overhead is the sum of previous components' heights
    plt.bar(x_pos, overhead_component, label="Overhead (time_wall - time_molorb)", bottom=init_component + molorb_calc_component, color="lightgrey", width=bar_width)

    plt.xlabel("Subdivision Factor")
    plt.ylabel("Wall Time (seconds)")
    plt.title(f"System '{system_name}': Runtime Breakdown vs. Subdivision Factor")
    plt.xticks(x_pos, subdivision_factors_str)  # Use string labels for factors on x-axis
    plt.legend()
    plt.grid(True, axis="y", linestyle="--", linewidth=0.5)  # y-axis grid for bar chart clarity
    plt.tight_layout()
    _save_and_close_plot(os.path.join(output_dir, f"{system_name}_runtime.png"))


def plot_convergence_for_system(system_name, system_df, subdivision_factors_num, ref_abssum, ref_charge, output_dir):
    """Plots convergence of abssum and charge relative to reference."""
    print(f"  Generating Convergence plot for {system_name}...")

    # Check if reference values are usable for calculating relative metrics
    can_plot_abssum_rel = abs(ref_abssum) > 1e-9
    can_plot_charge_rel = abs(ref_charge) > 1e-9

    if not can_plot_abssum_rel:
        print(f"  Warning: Reference 'abssum' is near zero ({ref_abssum}) for system '{system_name}'. Skipping its line in convergence plot.")
    if not can_plot_charge_rel:
        print(f"  Warning: Reference 'charge' is near zero ({ref_charge}) for system '{system_name}'. Skipping its line in convergence plot.")

    # Skip plot if no metrics can be meaningfully shown
    if not can_plot_abssum_rel and not can_plot_charge_rel:
        print("  Skipping convergence plot entirely as reference values for tracked metrics are unsuitable.")
        return

    plt.figure(figsize=(10, 6))
    plt.axhline(100, color="black", linestyle="--", linewidth=1, label=f"Reference ({REFERENCE_SUBDIVISION_FACTOR}=100%)")

    if can_plot_abssum_rel:
        rel_abssum = (system_df["abssum"] / ref_abssum) * 100
        plt.plot(subdivision_factors_num, rel_abssum, marker="o", linestyle="-", label="Absolute Sum (%)")
    if can_plot_charge_rel:
        rel_charge = (system_df["charge"] / ref_charge) * 100
        plt.plot(subdivision_factors_num, rel_charge, marker="s", linestyle=":", label="Charge (%)")

    _setup_common_line_plot_aesthetics(
        title=f"System '{system_name}': Convergence vs. Subdivision Factor", xlabel="Subdivision Factor", ylabel="Value Relative to Reference (%)", x_ticks_values=subdivision_factors_num
    )
    plt.gca().yaxis.set_major_formatter(PercentFormatter())
    plt.legend()  # Ensure legend is added after all plot elements
    plt.tight_layout()
    _save_and_close_plot(os.path.join(output_dir, f"{system_name}_convergence.png"))


def plot_mean_relative_error_for_system(system_name, non_ref_df, all_subdivision_factors_num, output_dir):
    """Plots mean relative error vs. subdivision factor."""
    print(f"  Generating Mean Relative Error plot for {system_name}...")
    if non_ref_df.empty:
        print(f"  Skipping Mean Relative Error plot for {system_name}: no non-reference data points.")
        return

    plt.figure(figsize=(10, 6))
    plt.plot(non_ref_df["subdivision_factor"], non_ref_df["mean_rel_err"], marker="o", linestyle="-")
    _setup_common_line_plot_aesthetics(
        title=f"System '{system_name}': Mean Relative Error vs. Subdivision Factor",
        xlabel="Subdivision Factor",
        ylabel="Mean Relative Error",
        x_ticks_values=all_subdivision_factors_num,  # Use all factors for consistent x-axis ticks
    )
    # Optional: Log scale if errors span large range and are positive
    # if not non_ref_df.empty and all(non_ref_df['mean_rel_err'].dropna() > 1e-9): # Check for positive values
    #     try:
    #         plt.yscale('log')
    #     except ValueError: # Handle cases where values might become non-positive after processing
    #         print(f"  Could not set log scale for Mean Relative Error plot for {system_name}.")
    #         pass

    plt.tight_layout()
    _save_and_close_plot(os.path.join(output_dir, f"{system_name}_mean_rel_error.png"))


def plot_mean_squared_error_for_system(system_name, non_ref_df, all_subdivision_factors_num, output_dir):
    """Plots mean squared error vs. subdivision factor."""
    print(f"  Generating Mean Squared Error plot for {system_name}...")
    if non_ref_df.empty:
        print(f"  Skipping Mean Squared Error plot for {system_name}: no non-reference data points.")
        return

    plt.figure(figsize=(10, 6))
    plt.plot(non_ref_df["subdivision_factor"], non_ref_df["mean_sq_err"], marker="o", linestyle="-", label="Mean Squared Error")
    _setup_common_line_plot_aesthetics(
        title=f"System '{system_name}': Mean Squared Error vs. Subdivision Factor", xlabel="Subdivision Factor", ylabel="Mean Squared Error", x_ticks_values=all_subdivision_factors_num
    )
    plt.legend()
    # Optional: Log scale (similar logic as above for mean_rel_err)
    plt.tight_layout()
    _save_and_close_plot(os.path.join(output_dir, f"{system_name}_mean_sq_error.png"))


def plot_max_squared_error_for_system(system_name, non_ref_df, all_subdivision_factors_num, output_dir):
    """Plots max squared error vs. subdivision factor."""
    print(f"  Generating Max Squared Error plot for {system_name}...")
    if non_ref_df.empty:
        print(f"  Skipping Max Squared Error plot for {system_name}: no non-reference data points.")
        return

    plt.figure(figsize=(10, 6))
    plt.plot(non_ref_df["subdivision_factor"], non_ref_df["max_sq_err"], marker="s", linestyle=":", label="Max Squared Error")
    _setup_common_line_plot_aesthetics(
        title=f"System '{system_name}': Max Squared Error vs. Subdivision Factor", xlabel="Subdivision Factor", ylabel="Max Squared Error", x_ticks_values=all_subdivision_factors_num
    )
    plt.legend()
    # Optional: Log scale (similar logic as above for mean_rel_err)
    plt.tight_layout()
    _save_and_close_plot(os.path.join(output_dir, f"{system_name}_max_sq_error.png"))


# --- Main Orchestration Function for a Single System ---


def generate_plots_for_system(system_name, system_df, output_dir):
    """Generates all plots for a single system."""
    print(f"--- Processing System: {system_name} ---")

    # Ensure data is sorted by subdivision factor for consistent plotting
    system_df = system_df.sort_values("subdivision_factor").reset_index(drop=True)

    # Extract subdivision factors for x-axis
    subdivision_factors_num = system_df["subdivision_factor"].values
    subdivision_factors_str = system_df["subdivision_factor"].astype(str).values

    # Plot 1: RAM Usage
    plot_ram_usage_for_system(system_name, system_df, subdivision_factors_num, output_dir)

    # Plot 2: Runtime Breakdown
    plot_runtime_breakdown_for_system(system_name, system_df, subdivision_factors_num, subdivision_factors_str, output_dir)

    # Reference data (guaranteed to exist as per problem statement)
    # .iloc[0] is safe because sorting ensures consistent selection if multiple ref runs exist,
    # and problem statement implies at least one exists.
    ref_data_row_df = system_df[system_df["subdivision_factor"] == REFERENCE_SUBDIVISION_FACTOR]
    if ref_data_row_df.empty:  # Should not happen based on problem statement but defensive
        print(f"  FATAL: Reference run for factor {REFERENCE_SUBDIVISION_FACTOR} not found for system '{system_name}' despite assumption. Skipping dependent plots.")
        return
    ref_data_row = ref_data_row_df.iloc[0]

    ref_abssum = ref_data_row["abssum"]
    ref_charge = ref_data_row["charge"]

    # Plot 3: Convergence
    plot_convergence_for_system(system_name, system_df, subdivision_factors_num, ref_abssum, ref_charge, output_dir)

    # Data excluding the reference point for error plots (error against self is 0)
    non_ref_df = system_df[system_df["subdivision_factor"] != REFERENCE_SUBDIVISION_FACTOR].copy()

    # Plot 4: Mean Relative Error
    plot_mean_relative_error_for_system(system_name, non_ref_df, subdivision_factors_num, output_dir)

    # Plot 5: Mean Squared Error (Separated from Max)
    plot_mean_squared_error_for_system(system_name, non_ref_df, subdivision_factors_num, output_dir)

    # Plot 6: Max Squared Error (Separated from Mean)
    plot_max_squared_error_for_system(system_name, non_ref_df, subdivision_factors_num, output_dir)


# --- Main Program Logic ---
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    try:
        with open(JSON_FILE_PATH, "r") as f:
            raw_data = json.load(f)
        print(f"Successfully loaded data from {JSON_FILE_PATH}")
    except FileNotFoundError:
        print(f"Error: File not found at {JSON_FILE_PATH}")
        return 1  # Indicate error
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {JSON_FILE_PATH}. Check its format.")
        return 1
    except Exception as e:
        print(f"An unexpected error occurred loading the JSON file: {e}")
        return 1

    if not isinstance(raw_data, list) or not raw_data:
        print("Error: JSON data is not a non-empty list of objects.")
        return 1

    df = pd.DataFrame(raw_data)

    # Validate required columns
    required_cols = ["system", "subdivision_factor", "max_ram_mb", "time_wall", "time_init", "time_molorb", "charge", "abssum", "mean_sq_err", "max_sq_err", "mean_rel_err"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: The JSON data is missing required columns: {', '.join(missing_cols)}")
        return 1

    # Group data by system and generate plots for each
    grouped_data = df.groupby("system")
    for system_name, system_df_group in grouped_data:
        # Pass a copy to avoid SettingWithCopyWarning, though current functions are mostly read-only on the df.
        generate_plots_for_system(system_name, system_df_group.copy(), OUTPUT_DIR)

    print(f"\nAll plots saved to directory: {OUTPUT_DIR}")
    return 0  # Indicate success


if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
