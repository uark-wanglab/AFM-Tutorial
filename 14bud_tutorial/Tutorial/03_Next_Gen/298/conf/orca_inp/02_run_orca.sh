#!/bin/bash

CWD=$(pwd)

targets=$(ls -d final*)

for item in $targets; do
    cd $item && sbatch runme.pinnacle ${item}.orca.inp
    cd $CWD
done
