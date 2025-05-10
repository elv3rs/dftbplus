import pathlib
import os
import shutil
import random

random.seed(0)

instance_target = 32
inpath = pathlib.Path("orig_inp")
outpath = pathlib.Path("inp")
dftb_executable = "/home/merlin/tmpfs/dftbuild/app/dftb+/dftb+"


def random_offset():
    scale = 1.0
    return [random.random() / scale for _ in range(3)]


def shift_geometry(geometry, offset):
    space = " " * 1
    distance = 4

    with open(geometry, "r") as f:
        lines = f.readlines()

    with open(geometry, "w") as f:
        for line in lines:
            parts = line.split()

            if len(line.split()) == 5:
                parts[2:5] = map(float, parts[2:5])
                parts[2] += offset[0]
                parts[3] += offset[1]
                parts[4] += offset[2]
                parts[2:5] = map(lambda x: f"{x:.10E}", parts[2:5])
            processed = space.join(p.ljust(distance) for p in parts)
            f.write(processed + "\n")


def run_dftb(target):
    cwd = os.getcwd()
    os.chdir(target)
    os.system(dftb_executable)
    os.chdir(cwd)


def main():

    offsets = [random_offset() for _ in range(instance_target)]

    for system in inpath.iterdir():
        for i, offset in enumerate(offsets):
            target = outpath / f"{system.stem}_{i+1:02d}"
            assert not target.exists()
            print("Creating", target)
            shutil.copytree(system, target)
            if i > 0:
                geometry = target / "geometry.gen"
                print("Shifting", geometry)
                shift_geometry(geometry, offset)

            print("Running DFTB+")
            run_dftb(target)


if __name__ == "__main__":
    main()
