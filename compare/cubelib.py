import numpy as np
import os
import typing
import re


DEBUG = False


def dprint(*args, **kwargs):
    """Debug print function to control debug output."""
    if DEBUG:
        print(*args, **kwargs)


def load_cube(filename: str) -> typing.Tuple[np.ndarray, dict]:
    """
    Loads a Gaussian cube file.

    Args:
        filename (str): The path to the cube file.

    Returns:
        tuple: A tuple containing:
            - data (np.ndarray): A 3D NumPy array of the volumetric data (NX, NY, NZ).
            - metadata (dict): A dictionary containing metadata.
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file format is inconsistent (e.g., wrong number of values).
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Cube file not found: {filename}")

    metadata = {}
    atoms = []
    data = None

    with open(filename, "r") as f:
        # --- Read Header ---
        metadata["comments"] = [f.readline().strip(), f.readline().strip()]
        line3 = f.readline().split()
        if len(line3) != 4:
            raise ValueError(f"Expected 4 values on line 3 (N_atoms, Ox, Oy, Oz), got {len(line3)} in {filename}: {line3}")
        metadata["n_atoms"] = int(line3[0])
        metadata["origin"] = np.array([float(x) for x in line3[1:]])

        units = "Bohr"
        axes_vecs = []
        dims = []
        for i, axis in enumerate(["X", "Y", "Z"]):
            line = f.readline().split()
            if len(line) != 4:
                raise ValueError(f"Expected 4 values on line {4+i} ({axis}-axis), got {len(line)} in {filename}: {line}")
            n_voxels = int(line[0])
            vec = np.array([float(x) for x in line[1:]])
            if i == 0 and n_voxels < 0:
                units = "Angstrom"
            dims.append(abs(n_voxels))
            axes_vecs.append(vec)
            metadata[f"N{axis}"] = abs(n_voxels)
            metadata[f"{axis}_vec"] = vec
        metadata["units"] = units
        NX, NY, NZ = dims

        for _ in range(metadata["n_atoms"]):
            line = f.readline().split()
            if len(line) != 5:
                raise ValueError(f"Expected 5 values for atom line, got {len(line)} in {filename}: {line}")
            atom_num = int(line[0])
            charge = float(line[1])
            coords = np.array([float(x) for x in line[2:]])
            atoms.append((atom_num, charge, coords))
        metadata["atoms"] = atoms

        # --- Read Volumetric Data ---
        data_string = f.read()
        try:
            # Attempt to handle potential extra whitespace/newlines robustly
            flat_data = np.fromstring(data_string, sep="\n")
            # If reading by newline fails significantly, try by space (more common)
            if flat_data.size < NX * NY * NZ * 0.9:  # Heuristic check
                print(f"Warning: Reading by newline seemed ineffective in {filename}, trying space separation.")
                flat_data = np.fromstring(data_string, sep=" ")

        except Exception as e:
            raise ValueError(f"Error parsing volumetric data in {filename}: {e}\nCheck for non-numeric values.")

        expected_size = NX * NY * NZ
        if flat_data.size != expected_size:
            # Try cleaning potential trailing empty strings from splitting
            flat_data = flat_data[flat_data != ""]
            if flat_data.size != expected_size:
                raise ValueError(f"Read {flat_data.size} data points in {filename}, expected {expected_size} ({NX}x{NY}x{NZ}). Check file format.")

        # Reshape according to Fortran ("F") order
        try:
            data = flat_data.reshape((NX, NY, NZ), order="F")
        except ValueError as e:
            raise ValueError(f"Reshaping error in {filename}: {e}. Read {flat_data.size}, expected {expected_size}.")

    return data, metadata


def set_system(system: str, n=1, input_hsd_file: str = "waveplot_in.hsd"):
    """
    Modifies the waveplot input HSD file to use the specified system's input files.
    Assumes paths are relative and contain 'inp/<system_name>/'.
    """
    dprint(f"Configuring for system: {system}")
    changed = False
    try:
        content = []
        with open(input_hsd_file, "r") as f:
            content = f.readlines()

        new_content = []
        pattern = re.compile(r"(.*\"inp\/)([^/]+)(\/.*\")")  # Matches '... "inp/<old_system>/..."'
        for line in content:
            new_content.append(line)

            match = pattern.search(line)
            if match and match.group(2) != system:
                new_line = pattern.sub(rf"\g<1>{system}\g<3>", line)
                new_content[-1] = new_line
                changed = True
                print(f"  Updated path: {new_line.strip()}")

            if "SubdivisionFactor" in line:
                if int(line.split("=")[-1]) != int(n):
                    changed = True
                    new_content[-1] = line.split("=")[0] + f"= {n}\n"

        if changed:
            with open(input_hsd_file, "w") as f:
                f.writelines(new_content)
            dprint(f"  '{input_hsd_file}' updated for system '{system}'.")
        else:
            dprint(f"  '{input_hsd_file}' already configured for system '{system}'.")

    except FileNotFoundError:
        print(f"Error: Input HSD file '{input_hsd_file}' not found.")
        raise
    except Exception as e:
        print(f"Error modifying '{input_hsd_file}': {e}")
        raise


def clean_output_files():
    """Removes files matching the given pattern in the current directory."""
    junkfiles = ("waveplot_pin.hsd", "wp-abs2.cube")
    [os.remove(junk) for junk in junkfiles if os.path.exists(junk)]
