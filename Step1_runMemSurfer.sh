#!/bin/bash

MEM_HOME=/home/lea/Software/memsurfer
export PYTHONPATH=$MEM_HOME/external/lib/python3.7/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=$MEM_HOME/external/lib:$MEM_HOME/external/lib64:$LD_LIBRARY_PATH

for PM in APM-dep ; do
  for mem in 1  ; do
    for fname in equil10nsBeforeEField ; do
      mkdir $PM/mem$mem/
      mkdir $PM/mem$mem/${fname}

      python memsurf.py $PM/trajs/  $PM mem$mem mem$mem.gro mem$mem.$fname.xtc $fname

    done
  done
done
