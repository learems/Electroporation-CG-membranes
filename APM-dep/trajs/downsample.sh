#!/bin/bash

source /usr/local/gromacs-2019.4/bin/GMXRC

mem=mem1
gmx trjconv -f ../../../../APM2x2/$mem/equilibration/conf.8.xtc -s ../../../../APM2x2/$mem/equilibration/conf.8.tpr -o $mem.equil10nsBeforeEField.xtc -b 40000 -e 50000 -skip 5 <<< "0"

