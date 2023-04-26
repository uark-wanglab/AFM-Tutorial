#!/bin/bash

echo "2 0" | gmx_d trjconv -f traj_comp.xtc -s ../topol.tpr -n ../index.ndx -split 80 -dt 80 -o final.gro -nzero 3 -center -pbc mol

ls -1 final* | xargs ./BUU-master

mv final*gro grofiles
