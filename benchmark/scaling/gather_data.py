"""Verteilt <N> Kohlenstoffatome zuffällig auf einen Würfel mit Kantenlänge <cube> """

import subprocess
import random
import shutil
from pathlib import Path
import json
import time
import sys

dftb = "/home/merlin/tmpfs/dftbuild/app/dftb+/dftb+"
waveplot = "/home/merlin/tmpfs/dftbuild/app/waveplot/waveplot"
VERSION = "28.04-sparse-22x"

VERBOSE_CLEAN = False

main_dir = Path(__file__).parent
dftb_dir = main_dir / "dftb"
dftb_cache = dftb_dir / "cache"
result_file = main_dir / "results.json"

Plot_every = 120


def main():
    if len(sys.argv) == 2:
        N = int(sys.argv[1])
        print("Loading dftb output for", N)
        prepare_dftb(N)
        run_dftb(N)
        exit()

    isWarm = False
    clean_dftb()
    clean_waveplot()
    lastPlot = time.time()

    results = {}
    if result_file.exists():
        with open(result_file, "r") as f:
            results = json.load(f)

    if VERSION not in results:
        results[VERSION] = {}

    target_N = [i**2 for i in range(1, 43)]
    for N in target_N:
        if str(N) in results.get(VERSION, {}):
            print(f"Skipping {N}")
            continue

        if not isWarm:
            print("Warming up")
            timeOnce(10)
            isWarm = True

        runtime = timeOnce(N)
        print(" -> Timing Result:", runtime)
        results[VERSION][N] = runtime
        with open(result_file, "w") as f:
            json.dump(results, f, indent=2)

        if time.time() - lastPlot > Plot_every:
            lastPlot = time.time()
            plot_results(results, dpi=100)

    print("Done!")
    plot_results(results, show=True)


def timeOnce(N):
    print(f"{N} Atoms")
    try:
        prepare_dftb(N)
        run_dftb(N)
        time.sleep(2)
        runtime = run_waveplot()
    finally:
        clean_dftb()
        clean_waveplot()
    return runtime


def get_gen_file(N):
    cube = 400
    string = f"{N}  C\nC\n"
    for i in range(0, N):
        xyz = [f"{random.uniform(0, cube):.11E}" for _ in "123"]
        string += f"{i+1}    1    " + "   ".join(xyz) + "\n"
    return string


def clean_dftb():
    if VERBOSE_CLEAN:
        print("Cleaning dftb")
    files = ["detailed.xml", "eigenvec.bin", "input.gen"]
    for file in files:
        f = dftb_dir / file
        if f.exists():
            f.unlink()
            if VERBOSE_CLEAN:
                print(f" -> Removed {file}")


def prepare_dftb(N):
    # print("Preparing dftb")
    string = get_gen_file(N)
    gen_file = dftb_dir / "input.gen"
    with open(gen_file, "w") as f:
        f.write(string)
    print(" -> Wrote input.gen")

    assert (dftb_dir / "dftb_in.hsd").exists()


def run_dftb(N):
    # Check if output is cached
    detailed = dftb_cache / f"{N:03d}_detailed.xml"
    eigenvec = dftb_cache / f"{N:03d}_eigenvec.bin"
    if detailed.exists() and eigenvec.exists():
        print("Using cached dftb calculation")
        shutil.copy(detailed, dftb_dir / "detailed.xml")
        shutil.copy(eigenvec, dftb_dir / "eigenvec.bin")
        return

    print("Running dftb")
    start = time.time()
    subprocess.run([dftb], cwd=dftb_dir, check=True, capture_output=True)
    print(" -> Done after", round(time.time() - start, 1), "s")
    # Ensure output was produced
    assert (dftb_dir / "detailed.xml").exists()
    assert (dftb_dir / "eigenvec.bin").exists()
    # cp to cache
    dftb_cache.mkdir(exist_ok=True)
    for file in ["detailed.xml", "eigenvec.bin"]:
        src = dftb_dir / file
        dst = dftb_cache / f"{N:03d}_{file}"
        shutil.copy(src, dst)


def clean_waveplot():
    if VERBOSE_CLEAN:
        print("Cleaning waveplot")
    files = ["waveplot_pin.hsd", "wp-1-1-4-abs2.cube", "wp-abs2.cube", "wp-abs2diff.cube", "wp-1-1-4-real.cube"]
    for file in files:
        f = main_dir / file
        if f.exists():
            f.unlink()
            if VERBOSE_CLEAN:
                print(f" -> Removed {file}")


def run_waveplot():
    # time waveplot execution. Return median and stddev
    """
    real	0m0.003s
    user	0m0.000s
    sys	    0m0.003s
    """
    print("Timing waveplot")
    p = subprocess.run(["time", "-p", waveplot], cwd=main_dir, capture_output=True, check=True)
    # print(p.stderr.decode())
    real, runtime_s = p.stderr.decode().split("\n")[-4].split()
    assert real == "real"
    return float(runtime_s)


def plot_results(results, show=False, dpi=300):
    import matplotlib.pyplot as plt

    for version, data in results.items():
        x = [float(i) for i in data.keys()]
        y = [float(i) for i in data.values()]
        plt.plot(x, y, label=version, marker="o", linestyle="")

    plt.legend()
    plt.xlabel("Atoms [N]")
    plt.ylabel("Runtime [s]")
    plt.savefig("runtime.png", dpi=dpi)
    print("Saved plot to runtime.png")
    if show:
        plt.show()
    plt.close()


if __name__ == "__main__":
    main()
