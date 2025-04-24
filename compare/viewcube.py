import numpy as np
import os
import time
import typing
import matplotlib.pyplot as plt
import napari


def load_cube(filename: str) -> typing.Tuple[np.ndarray, dict]:
    """
    Loads a Gaussian cube file.

    Args:
        filename (str): The path to the cube file.

    Returns:
        tuple: A tuple containing:
            - data (np.ndarray): A 3D NumPy array of the volumetric data (NX, NY, NZ).
            - metadata (dict): A dictionary containing metadata:
                - 'comments' (list): The first two comment lines.
                - 'n_atoms' (int): Number of atoms.
                - 'origin' (np.ndarray): 1x3 array of the volume origin (X, Y, Z).
                - 'NX', 'NY', 'NZ' (int): Number of voxels along each axis.
                - 'X_vec', 'Y_vec', 'Z_vec' (np.ndarray): 1x3 arrays of the voxel vectors.
                                                        Length is voxel dimension, direction is axis direction.
                                                        Units are Bohr if N is positive, Angstrom if negative.
                - 'units' (str): 'Bohr' or 'Angstrom', based on the sign of NX (conventionally).
                - 'atoms' (list): A list of tuples, one for each atom:
                                  (atom_number, charge, [coord_x, coord_y, coord_z])
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file format is inconsistent (e.g., wrong number of values).
    """
    metadata = {}
    atoms = []
    data = None

    with open(filename, "r") as f:
        # --- Read Header ---
        # Comments
        metadata["comments"] = [f.readline().strip(), f.readline().strip()]

        # Number of atoms and origin
        line3 = f.readline().split()
        if len(line3) != 4:
            raise ValueError(f"Expected 4 values on line 3 (N_atoms, Ox, Oy, Oz), got {len(line3)}: {line3}")
        metadata["n_atoms"] = int(line3[0])
        metadata["origin"] = np.array([float(x) for x in line3[1:]])

        # Voxel grid dimensions and vectors
        units = "Bohr"  # Default assumption
        axes_vecs = []
        dims = []
        for i, axis in enumerate(["X", "Y", "Z"]):
            line = f.readline().split()
            if len(line) != 4:
                raise ValueError(f"Expected 4 values on line {4+i} ({axis}-axis), got {len(line)}: {line}")
            n_voxels = int(line[0])
            vec = np.array([float(x) for x in line[1:]])

            # Determine units from the sign (conventionally use X-axis sign)
            if i == 0 and n_voxels < 0:
                units = "Angstrom"

            dims.append(abs(n_voxels))
            axes_vecs.append(vec)
            metadata[f"N{axis}"] = abs(n_voxels)
            metadata[f"{axis}_vec"] = vec

        metadata["units"] = units
        NX, NY, NZ = dims

        # Atom information
        for _ in range(metadata["n_atoms"]):
            line = f.readline().split()
            if len(line) != 5:
                raise ValueError(f"Expected 5 values for atom line, got {len(line)}: {line}")
            atom_num = int(line[0])
            charge = float(line[1])
            coords = np.array([float(x) for x in line[2:]])
            atoms.append((atom_num, charge, coords))
        metadata["atoms"] = atoms

        # --- Read Volumetric Data ---
        # Read the rest of the file into a single string
        data_string = f.read()

        # Parse the flat data list
        try:
            flat_data = np.fromstring(data_string, sep=" ")
        except Exception as e:
            raise ValueError(f"Error parsing volumetric data: {e}\nCheck for non-numeric values.")

        expected_size = NX * NY * NZ
        if flat_data.size != expected_size:
            raise ValueError(f"Read {flat_data.size} data points, expected {expected_size} ({NX}x{NY}x{NZ})")

        # Reshape according to Fortran ("F") order (X outer, Y middle, Z inner loop)
        data = flat_data.reshape((NX, NY, NZ), order="F")

    return data, metadata


def draw_slice(volumetric_data, name):
    half = volumetric_data.shape[0] // 2
    slice_to_plot = volumetric_data[half, :, :]
    plt.imshow(slice_to_plot.T, cmap="viridis", origin="lower", aspect="auto")
    plt.colorbar()
    plt.title(f"Slice at ix = {volumetric_data.shape[0] // 2}")
    plt.xlabel("Y index")
    plt.ylabel("Z index")
    filename = f"{name}.png"
    plt.savefig(f"{name}.png")
    plt.close()
    print(f"Slice saved as {filename}")


def setSubdivisionFactor(n=1):
    with open("subdivision_factor.hsd", "w") as f:
        f.write(f"SubdivisionFactor = {n}")


def compare():
    runs = {}
    for i in range(1, 50):
        setSubdivisionFactor(i)
        cmd = "time /home/merlin/tmpfs/dftbuild/app/waveplot/waveplot"
        start = time.perf_counter()
        print(f"Running waveplot with subdivision factor {i}")
        print(cmd)
        os.system(cmd)
        elapsed = round(time.perf_counter() - start, 2)
        print(f"Elapsed time: {elapsed} seconds")
        if os.path.exists("wp-abs2.cube"):
            os.rename("wp-abs2.cube", f"wpMain-{i}.cube")
        else:
            break
        volMain, meta = load_cube(f"wpMain-{i}.cube")
        volChunked, meta = load_cube("wpChunked.cube")
        diffSq = get_sqdiff_sum("wpMain.cube", f"wpMain-{i}.cube")
        sqdiffsum = np.sum(diffSq)
        runs[i] = [elapsed, sqdiffsum]
        print(f"Sum of squared differences: {sqdiffsum}")
    print("Results:")
    print(f"{'Subdivision Factor':<20} {'Elapsed Time (s)':<20} {'Sum of Squared Differences':<30}")
    for i, (elapsed, sqdiffsum) in runs.items():
        print(f"{i:<20} {elapsed:<20} {sqdiffsum:<30}")


def get_sqdiff_sum(volMain, volChunked, view=False):
    volMain, meta = load_cube(volMain)
    volChunked, meta = load_cube(volChunked)
    diffSq = (volMain - volChunked) ** 2

    print(f"Sum of squared differences: {np.sum(diffSq)}")
    if view:
        draw_slice(diffSq, "diffSq")

        viewer = napari.Viewer()
        viewer.add_image(volMain, name="volMain", colormap="turbo", rendering="attenuated_mip")
        viewer.add_image(volChunked, name="volChunked", colormap="turbo", rendering="attenuated_mip")
        viewer.add_image(np.sqrt(diffSq), name="diffSq", colormap="turbo", rendering="attenuated_mip")
        napari.run()

    return diffSq


def view_percent_diff(volMain, volChunked):
    volMain, meta = load_cube(volMain)
    volChunked, meta = load_cube(volChunked)
    diff = np.abs((volMain - volChunked) / np.where(np.abs(volMain) > 1e-6, volMain, 1.0))
    avg_diff = np.mean(diff)
    print(f"Average percent difference: {avg_diff}")

    draw_slice(diff, "diff Percent")
    print("Starting napari")
    napari.view_image(volMain, name="volMain", colormap="turbo", rendering="attenuated_mip")
    napari.view_image(volChunked, name="volChunked", colormap="turbo", rendering="attenuated_mip")
    napari.view_image(diff, name="diff Percent", colormap="turbo", rendering="attenuated_mip")
    napari.run()


def main():
    # data = {}
    # for name in ["wpMain.cube", "wpChunked.cube"]:
    #    vol, meta = data[f"{name}.volume"], data[f"{name}.meta"] = load_cube(name)
    #    print(vol.shape)

    diffSq = get_sqdiff_sum("wpMain.cube", "wpChunked.cube", view=True)
    draw_slice(diffSq, "diffSq")


if __name__ == "__main__":
    view_percent_diff("wpMain.cube", "wpChunked.cube")
    # compare()
    # main()
