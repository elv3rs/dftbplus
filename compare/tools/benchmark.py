import numpy as np
import os
import typing
import subprocess
import re
import json
import sys

from cubelib import load_cube, set_system, clean_output_files

# --- Configuration ---
WAVEPLOT_EXECUTABLE = "/home/merlin/tmpfs/dftbuild/app/waveplot/waveplot"
RESULTS_JSON_FILE = "output/results.json"
TIMEOUT = 600

# --- ANSI Colors ---
RESET = "\033[0m"
RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
BLUE = "\033[94m"
MAGENTA = "\033[95m"


# --- Helper Functions ---
def _print_status(message: str, color: str = BLUE):
    print(f"{color}{message}{RESET}")


def _print_error(message: str):
    print(f"{RED}ERROR: {message}{RESET}", file=sys.stderr)


def _print_warning(message: str):
    print(f"{YELLOW}WARNING: {message}{RESET}", file=sys.stderr)


# --- Comparison and Saving ---
def calculate_errors(approx: np.ndarray, ref: np.ndarray) -> dict:
    """Calculates various error metrics between two numpy arrays."""
    diff = approx - ref
    sqdiff = diff**2
    norm_ref = np.linalg.norm(ref)
    norm_diff = np.linalg.norm(diff)
    rel_err = norm_diff / norm_ref if norm_ref > 1e-15 else 0.0
    return {
        "mean_sq_err": np.mean(sqdiff),
        "max_sq_err": np.max(sqdiff),
        "mean_rel_err": rel_err,
        "abssum": np.sum(np.abs(approx)),
    }


def save_results(results_dict: dict, filename: str = RESULTS_JSON_FILE):
    """Appends a results dictionary to a JSON file."""
    results_list = []
    if os.path.exists(filename):
        try:
            with open(filename, "r") as f:
                content = f.read()
                if content.strip():
                    results_list = json.loads(content)
                if not isinstance(results_list, list):
                    _print_warning(f"Content of {filename} is not a list. Overwriting.")
                    results_list = []
        except json.JSONDecodeError:
            _print_warning(f"Could not decode JSON from {filename}. Overwriting.")
            results_list = []
        except Exception as e:
            _print_error(f"Reading {filename}: {e}. Starting fresh.")
            results_list = []

    results_list.append(results_dict)
    with open(filename, "w") as f:
        json.dump(results_list, f, indent=4)


# --- Program Execution Logic ---


def _parse_waveplot_output(output: str, subdivision_factor: int) -> dict:
    """Parses the combined stdout/stderr from waveplot and /usr/bin/time."""
    try:
        parsed_data = {
            "max_ram_mb": int(re.search(r"MaxMemoryKB:\s*(\d+)", output).group(1)) / 1000,
            "time_wall": float(re.search(r"WallTimeSecs:\s*([\d.]+)", output).group(1)),
            "charge": float(re.search(r"Total charge:\s*([\d.eE+-]+)", output).group(1)),
            "time_molorb": round(float(re.search(r"MolorbTime:\s*([\d.eE+-]+)", output).group(1)), 2),
            "time_init": 0.0,
        }
        # optional field
        init_match = re.search(r"InitTime:\s*([\d.eE+-]+)", output)
        if init_match:
            parsed_data["time_init"] = round(float(init_match.group(1)), 2)

        return parsed_data

    except AttributeError as e:
        raise ValueError("Missing required field in waveplot output.") from e


def run_waveplot(system: str, subdivision_factor: int) -> typing.Tuple[dict, np.ndarray]:
    """
    Runs waveplot, parses output, handles cube file, returns info and data.
    """
    run_info = {"system": system, "subdivision_factor": subdivision_factor}

    # --- Prepare ---
    _print_status(f"Running: System='{system}', Subdivision={subdivision_factor}", BLUE)
    set_system(system, subdivision_factor)
    clean_output_files()
    output_cube_file = "wp-abs2.cube"

    # --- Execute ---
    cmd = f'/usr/bin/time -f "MaxMemoryKB: %M\\nWallTimeSecs: %e" {WAVEPLOT_EXECUTABLE}'
    try:
        result = subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True, timeout=TIMEOUT)
    except FileNotFoundError:
        _print_error(f"Command '/usr/bin/time' or executable '{WAVEPLOT_EXECUTABLE}' not found.")
        raise
    except subprocess.TimeoutExpired:
        _print_error(f"Command timed out after {TIMEOUT} seconds.")
        raise ValueError(f"Waveplot timed out (>{TIMEOUT}s)")

    # --- Parse Output & Handle Cube ---
    output = result.stdout + "\n" + result.stderr
    try:
        parsed_metrics = _parse_waveplot_output(output, subdivision_factor)
        run_info.update(parsed_metrics)
    except ValueError as e:
        _print_error(f"Parsing failed: {e}")
        print(f"--- Waveplot Output ---\n{output}\n--- End Output ---")
        raise  # Re-raise after printing context

    if not os.path.exists(output_cube_file):
        _print_error(f"Expected output file '{output_cube_file}' not found.")
        print(f"--- Waveplot Output ---\n{output}\n--- End Output ---")
        raise ValueError(f"Output file {output_cube_file} not generated.")

    volume_data, _ = load_cube(output_cube_file)

    clean_output_files()

    return run_info, volume_data


# --- Main Execution ---
def main():
    systems = list(os.listdir("inp"))
    # systems = ["H"]

    _print_status(f"Found systems: {systems}", GREEN)

    # Create output directory
    os.makedirs("output", exist_ok=True)
    if os.path.exists(RESULTS_JSON_FILE):
        _print_warning(f"{RESULTS_JSON_FILE} exists, results will be appended.")

    for system in systems:
        print(f"\n{MAGENTA}===== Processing System: {system} ====={RESET}")
        reference_vol = None
        reference_info = None

        try:
            # Warmup run
            _print_status("Warmup Run (S=0)", YELLOW)
            run_waveplot(system, 0)

            # Get reference output (Subdivision 0)
            _print_status("Reference Run (S=0)", YELLOW)
            reference_info, reference_vol = run_waveplot(system, 0)

            # Calculate "errors" for reference vs reference (should be ~0)
            ref_errors = calculate_errors(reference_vol, reference_vol)
            reference_info.update(ref_errors)
            save_results(reference_info, RESULTS_JSON_FILE)
            print(f"{GREEN}Reference (S=0) saved.{RESET} Errors: mse={ref_errors['mean_sq_err']:.2e}, rel={ref_errors['mean_rel_err']:.2e}")

            # Run with increasing subdivision factors
            for i in [-1] + list(range(1, 100)):
                try:
                    run_info, approx_vol = run_waveplot(system, i)

                    errors = calculate_errors(approx_vol, reference_vol)
                    run_info.update(errors)
                    save_results(run_info, RESULTS_JSON_FILE)
                    print(f"{GREEN} S={i}: OK{RESET} | Time={run_info['time_wall']:.2f}s | RAM={run_info['max_ram_mb']:.1f}MB | MSE={errors['mean_sq_err']:.2e} | REL_ERR={errors['mean_rel_err']:.2e}")

                except Exception as e:
                    _print_error(f"Error during S={i} for {system}: {e.__class__.__name__}: {e}")
                    _print_warning(f"Stopping subdivision ramp for {system}.")
                    break

        except (ValueError, FileNotFoundError, subprocess.CalledProcessError) as e:
            _print_error(f"Critical error during reference/warmup for {system}: {e}")
            _print_warning(f"Skipping system {system}.")
            continue  # Skip to next system
        except Exception as e:  # Catch any other unexpected error during ref/warmup
            _print_error(f"Unexpected error during reference/warmup for {system}: {e.__class__.__name__}: {e}")
            _print_warning(f"Skipping system {system}.")
            continue

        print(f"{MAGENTA}===== Finished System: {system} ====={RESET}")

    print(f"\n{GREEN}Benchmarking complete. Results saved in: {RESULTS_JSON_FILE}{RESET}")


if __name__ == "__main__":
    main()
