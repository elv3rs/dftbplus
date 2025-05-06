import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import os

# --- Configuration ---
JSON_FILE_PATH = "output/results.json"
OUTPUT_DIR = "output/plots"
REFERENCE_SUBDIVISION_FACTOR = 0
# --- End Configuration ---


def plot_system_data(system_name, system_df, output_dir):
    """Generates all plots for a single system."""

    print(f"--- Processing System: {system_name} ---")

    # Ensure data is sorted by subdivision factor for plotting
    system_df = system_df.sort_values("subdivision_factor").reset_index(drop=True)

    # Use subdivision factor as string for categorical bar plots if preferred,
    # or as numbers for line plots. We'll use numbers for line plots
    # and potentially strings for bar x-labels.
    subdivision_factors_num = system_df["subdivision_factor"].values
    subdivision_factors_str = system_df["subdivision_factor"].astype(str).values

    # --- 1. RAM Usage Plot ---
    print("  Generating RAM plot...")
    fig_ram, ax_ram = plt.subplots(figsize=(10, 6))
    ax_ram.plot(subdivision_factors_num, system_df["max_ram_mb"], marker="o", linestyle="-", color="blue")
    ax_ram.set_xlabel("Subdivision Factor")
    ax_ram.set_ylabel("Max RAM Usage (MB)")
    ax_ram.set_title(f"System '{system_name}': RAM Usage vs. Subdivision Factor")
    ax_ram.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax_ram.set_xticks(subdivision_factors_num)  # Ensure ticks for all factors
    ax_ram.margins(x=0.05)  # Add a little space at the ends
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{system_name}_ram_usage.png"))
    plt.close(fig_ram)

    # --- 2. Runtime Plot (Stacked Bar) ---
    print("  Generating Runtime plot...")
    fig_time, ax_time = plt.subplots(figsize=(10, 6))

    # Calculate overhead, ensuring it's not negative due to float precision
    time_init = system_df["time_init"].values
    time_molorb = system_df["time_molorb"].values
    time_wall = system_df["time_wall"].values
    time_overhead = np.maximum(0, time_wall - time_init - time_molorb)  # Clip negative values

    bar_width = 0.6  # Adjust for visual preference

    # Use numeric factors for positioning but string labels
    x_pos = np.arange(len(subdivision_factors_num))

    ax_time.bar(x_pos, time_init, label="Initialization (time_init)", color="skyblue", width=bar_width)
    ax_time.bar(x_pos, time_molorb - time_init, label="MolOrb Calculation (molorb -init)", bottom=time_init, color="lightcoral", width=bar_width)
    ax_time.bar(x_pos, time_overhead, label="Overhead (wall - init - molorb)", bottom=time_init + time_molorb, color="lightgrey", width=bar_width)

    ax_time.set_xlabel("Subdivision Factor")
    ax_time.set_ylabel("Wall Time (seconds)")
    ax_time.set_title(f"System '{system_name}': Runtime Breakdown vs. Subdivision Factor")
    ax_time.set_xticks(x_pos)
    ax_time.set_xticklabels(subdivision_factors_str)  # Use string labels
    ax_time.legend()
    ax_time.grid(True, axis="y", linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{system_name}_runtime.png"))
    plt.close(fig_time)

    # --- Convergence and Error Plots (Require Reference Run) ---
    reference_run = system_df[system_df["subdivision_factor"] == REFERENCE_SUBDIVISION_FACTOR]

    if reference_run.empty:
        print(f"  Warning: No reference run (subdivision_factor={REFERENCE_SUBDIVISION_FACTOR}) found for system '{system_name}'. Skipping convergence/error plots.")
    else:
        # Ensure only one reference run is used if multiple exist (take first)
        ref_data = reference_run.iloc[0]
        ref_abssum = ref_data["abssum"]
        ref_charge = ref_data["charge"]

        # Check for zero reference values to avoid division by zero
        can_plot_abssum = abs(ref_abssum) > 1e-9
        can_plot_charge = abs(ref_charge) > 1e-9

        if not can_plot_abssum:
            print(f"  Warning: Reference 'abssum' is near zero ({ref_abssum}) for system '{system_name}'. Skipping its convergence plot.")
        if not can_plot_charge:
            print(f"  Warning: Reference 'charge' is near zero ({ref_charge}) for system '{system_name}'. Skipping its convergence plot.")

        # --- 3. Convergence Plot (abssum, charge) ---
        if can_plot_abssum or can_plot_charge:
            print("  Generating Convergence plot...")
            fig_conv, ax_conv = plt.subplots(figsize=(10, 6))
            ax_conv.axhline(100, color="black", linestyle="--", linewidth=1, label=f"Reference ({REFERENCE_SUBDIVISION_FACTOR}=100%)")

            if can_plot_abssum:
                rel_abssum = (system_df["abssum"] / ref_abssum) * 100
                ax_conv.plot(subdivision_factors_num, rel_abssum, marker="o", linestyle="-", label="Absolute Sum (%)")

            if can_plot_charge:
                rel_charge = (system_df["charge"] / ref_charge) * 100
                ax_conv.plot(subdivision_factors_num, rel_charge, marker="s", linestyle=":", label="Charge (%)")

            ax_conv.set_xlabel("Subdivision Factor")
            ax_conv.set_ylabel("Value Relative to Reference (%)")
            ax_conv.set_title(f"System '{system_name}': Convergence vs. Subdivision Factor")
            ax_conv.legend()
            ax_conv.grid(True, which="both", linestyle="--", linewidth=0.5)
            ax_conv.yaxis.set_major_formatter(PercentFormatter())
            ax_conv.set_xticks(subdivision_factors_num)
            ax_conv.margins(x=0.05)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{system_name}_convergence.png"))
            plt.close(fig_conv)
        else:
            print("  Skipping convergence plot as both reference values are near zero.")

        # --- 4. Relative Error Plot ---
        print("  Generating Mean Relative Error plot...")
        fig_rel_err, ax_rel_err = plt.subplots(figsize=(10, 6))
        # Exclude reference point (error is definitionally 0) if desired, or plot it
        non_ref_df = system_df[system_df["subdivision_factor"] != REFERENCE_SUBDIVISION_FACTOR]
        ax_rel_err.plot(non_ref_df["subdivision_factor"], non_ref_df["mean_rel_err"], marker="o", linestyle="-")
        # Optional: plot reference point at 0
        # ax_rel_err.plot(REFERENCE_SUBDIVISION_FACTOR, 0, marker='o', color='red', label="Reference")

        ax_rel_err.set_xlabel("Subdivision Factor")
        ax_rel_err.set_ylabel("Mean Relative Error")
        ax_rel_err.set_title(f"System '{system_name}': Mean Relative Error vs. Subdivision Factor")
        ax_rel_err.grid(True, which="both", linestyle="--", linewidth=0.5)
        # Optional: Log scale if errors span large range and are positive
        # try:
        #     if all(non_ref_df['mean_rel_err'] > 0):
        #          ax_rel_err.set_yscale('log')
        # except Exception: pass # Handle empty data if needed
        ax_rel_err.set_xticks(subdivision_factors_num)
        ax_rel_err.margins(x=0.05)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{system_name}_mean_rel_error.png"))
        plt.close(fig_rel_err)

        # --- 5. Squared Error Plot ---
        print("  Generating Squared Errors plot...")
        fig_sq_err, ax_sq_err = plt.subplots(figsize=(10, 6))
        # Exclude reference point if desired
        ax_sq_err.plot(non_ref_df["subdivision_factor"], non_ref_df["mean_sq_err"], marker="o", linestyle="-", label="Mean Squared Error")
        ax_sq_err.plot(non_ref_df["subdivision_factor"], non_ref_df["max_sq_err"], marker="s", linestyle=":", label="Max Squared Error")
        # Optional: plot reference point at 0
        # ax_sq_err.plot(REFERENCE_SUBDIVISION_FACTOR, 0, marker='o', color='red', label="Reference")

        ax_sq_err.set_xlabel("Subdivision Factor")
        ax_sq_err.set_ylabel("Squared Error")
        ax_sq_err.set_title(f"System '{system_name}': Squared Errors vs. Subdivision Factor")
        ax_sq_err.legend()
        ax_sq_err.grid(True, which="both", linestyle="--", linewidth=0.5)
        # Optional: Log scale if errors span large range and are positive
        # try:
        #      if all(non_ref_df['mean_sq_err'] > 0) and all(non_ref_df['max_sq_err'] > 0):
        #         ax_sq_err.set_yscale('log')
        # except Exception: pass # Handle empty data
        ax_sq_err.set_xticks(subdivision_factors_num)
        ax_sq_err.margins(x=0.05)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{system_name}_sq_errors.png"))
        plt.close(fig_sq_err)


def main():
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load data from JSON file
    try:
        with open(JSON_FILE_PATH, "r") as f:
            raw_data = json.load(f)
        print(f"Successfully loaded data from {JSON_FILE_PATH}")
    except FileNotFoundError:
        print(f"Error: File not found at {JSON_FILE_PATH}")
        exit(1)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {JSON_FILE_PATH}. Check its format.")
        exit(1)
    except Exception as e:
        print(f"An unexpected error occurred loading the JSON file: {e}")
        exit(1)

    if not isinstance(raw_data, list) or not raw_data:
        print("Error: JSON data is not a non-empty list of objects.")
        exit(1)

    # Convert to DataFrame
    df = pd.DataFrame(raw_data)

    # Basic validation of required columns
    required_cols = ["system", "subdivision_factor", "max_ram_mb", "time_wall", "time_init", "time_molorb", "charge", "abssum", "mean_sq_err", "max_sq_err", "mean_rel_err"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: The JSON data is missing required columns: {', '.join(missing_cols)}")
        exit(1)

    # Group data by system
    grouped_data = df.groupby("system")

    # Generate plots for each system
    for system_name, system_df in grouped_data:
        plot_system_data(system_name, system_df.copy(), OUTPUT_DIR)  # Pass a copy to avoid SettingWithCopyWarning

    print(f"\nAll plots saved to directory: {OUTPUT_DIR}")


# --- Main Execution ---
if __name__ == "__main__":
    main()
