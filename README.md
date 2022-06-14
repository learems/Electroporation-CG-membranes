# Electroporation-CGmem-MemSurfer
This repository contains a collection of scripts used for the analysis of trajectories from molecular dynamics simulation of electroporation of coarse-grained Martini membranes with complex lipid composition. The results are published at [eLife](https://doi.org/10.7554/eLife.74773). 

The scripts were run on Ubuntu 20.04 with Anaconda-based python 3.7 suite and Matlab R2021a. 

The analysis proceeds in 5 steps: 
1. Extract local membrane properties from each frame of a trajectory using the [MemSurfer tool](https://github.com/LLNL/MemSurfer) by running Step1_runMemSurfer.sh. There is an example trajectory in the folder APM-dep/trajs/. This trajectory contains 20 frames saved every 500 ps. The original analysis was based on trajectories with 100 frames saved every 100 ps. The script writes the extracted properties into txt files in directory APM-dep/mem1/equil10nsBeforeEField/. The text files for this example can be found in APM-dep/mem1.zip
2. Load data into Matlab using Step2_LoadDataToMatlab.m
3. Extract local values from nonporated and porated locations using Step3_Extractvalues.m
4. Plot histograms and calculate the distance between histograms (KL divergence) using Step4_HistsAndKLD.m
5. Make a surface plot of the local membrane properties using Step5_MakeSurfacePlot.m (optional)
6. Export data into tables to be used for machine learning analysis using Step6_ExportTablesForML.m The scripts used for ML analysis are found in https://github.com/learems/Electroporation-CGmem-MachineLearning 
