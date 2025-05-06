import cubelib
import napari
import sys
import os
import numpy as np

# Usage check
if len(sys.argv) < 2:
    print("Usage: python view.py <cube_file1> [cube_file2]")
    sys.exit(1)

cube_file1 = sys.argv[1]

if not os.path.exists(cube_file1):
    print(f"File not found: {cube_file1}")
    sys.exit(1)

vol1, _ = cubelib.load_cube(cube_file1)

nameA = cube_file1.split(".")[0]
# Start Napari viewer
viewer = napari.Viewer()

viewer.add_image(vol1, name=nameA, colormap="turbo", rendering="attenuated_mip")

# Optional second file
if len(sys.argv) > 2:
    cube_file2 = sys.argv[2]
    if not os.path.exists(cube_file2):
        print(f"File not found: {cube_file2}")
        sys.exit(1)

    nameB = cube_file2.split(".")[0]

    vol2, _ = cubelib.load_cube(cube_file2)
    viewer.add_image(vol2, name=nameB, colormap="turbo", rendering="attenuated_mip")

    # Compute and show difference
    diff = vol1.astype(np.float32) - vol2.astype(np.float32)
    viewer.add_image(diff, name="Diff (A - B)", colormap="turbo", rendering="attenuated_mip")

napari.run()
