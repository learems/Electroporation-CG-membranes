#!/bin/bash

MEM_HOME=/home/lea/Software/memsurfer
export PYTHONPATH=$MEM_HOME/external/lib/python3.7/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=$MEM_HOME/external/lib:$MEM_HOME/external/lib64:$LD_LIBRARY_PATH

#for PM in APM2x2 BPM2x2 APM2x2-hyp BPM2x2-hyp APM2x2-40us ; do
for PM in APM2x2-40us ; do
      	for mem in 1 ; do
#    for fname in AN_EF0mV_nores_run1 AN_EF0mV_nores_run2 AN_EF0mV_nores_run3 AN_EF0mV_nores_run4 AN_EF0mV_nores_run5 ; do
    for fname in frame40eq ; do
#    for fname in frame.tp+0 ; do
        mkdir $PM
        mkdir $PM/mem$mem/
        mkdir $PM/mem$mem/${fname}

	if [[ $fname == *"EF0mV"* ]]; then
          python memsurf.py ../../$PM/mem$mem/${fname}   $PM mem$mem conf.8.gro conf.9.every20ps.xtc $fname
        elif [[ $fname == *"eq"* ]]; then
          python memsurf.py ../../$PM/extracted-frames/  $PM mem$mem mem$mem.gro mem$mem.$fname.xtc $fname
        elif [[ $fname == *"tp+0"* ]]; then
          python memsurf.poreframe.py ../../$PM/extracted-frames-new/  $PM mem$mem mem$mem.gro mem$mem.$fname.xtc $fname
	fi

    done
  done
done
