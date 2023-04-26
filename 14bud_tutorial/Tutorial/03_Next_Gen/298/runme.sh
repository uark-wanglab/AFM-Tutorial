#!/bin/bash

if [ $NP == 32 ]; then
    gmx_d grompp -f grompp.mdp -p BUU_water.top -c conf.gro -n index.ndx 2> gerr.txt 1> gout.txt
    gmx_d mdrun -ntmpi 16 -table BUU.xvg -tableb BUU_b0.xvg 2> rerr.txt 1> rout.txt
fi

if [ $NP == 40]; then
    gmx_d grompp -f grompp.mdp -p BUU_water.top -c conf.gro -n index.ndx 2> gerr.txt 1> gout.txt
    gmx_d mdrun -ntmpi 16 -table BUU.xvg -tableb BUU_b0.xvg 2> rerr.txt 1> rout.txt
fi
    

