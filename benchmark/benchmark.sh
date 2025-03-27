#!/bin/bash

WP="/home/merlin/tmpfs/dftbuild/app/waveplot/waveplot"

time_waveplot() {        
    cd "$1" || exit 1
    echo "Timing $1 system:"
    hyperfine "$WP" \
        --export-json "../$1-first-dev.json" \
        --cleanup "rm -f waveplot_pin.hsd wp-1-1-4-abs2.cube wp-abs2.cube wp-abs2diff.cube wp-1-1-4-real.cube"

    cd ..
}

time_waveplot small
time_waveplot medium
time_waveplot large

