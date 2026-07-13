This is a software that implements the fit of force fields using the linear model SNAP (Spectral Neighbour Analysis Potential) for both periodic and non-periodic systems. The software can also be used to perform molecular dynamics and geometrical optimization, 
and in the case in which a training set is not available yet, FieldsForge can be driven by an external program present in this repo to perform active learning for the optimal construction of rapid training set (1).
For isolated systems, the fit of electric dipoles can also be performed. Its compilation requires the software suite MolForge (https://github.com/LunghiGroup/MolForge) and LAMMPS (https://github.com/lammps/lammps/) tested on version 29 Oct 2020). FieldsForge is known to compile on Linux machines with GNU compilers (v12.2.0) and openmpi(4.1.1).

FieldsForge depends on the linear algebra routines Lapack (v3.8.0) and Blas(v3.8.0).

BASIC INSTRUCTIONS

The program can be compiled with the available Makefile.

The software is executed using:

$ ./fforge snap.input

"snap.input" is the text file containing the instructions for the program. In the repo is provided a template with explanations of the meaning of the single variables as comments.
To test the correct execution of the program, please go the "test_benzene" folder.
In the folder "test_benzene", the data for a dataset of 79 molecules of benzene are provided. In the README.txt of that directory the formats for the single input files are explained.

INSTRUCTIONS FOR ACTIVE LEARNING

Here the instructions to perform active learning on HPC facilities is given calling the quantum chemistry code Orca. The scripts can be easily adapted to perform also in local and with any
other quantum chemistry software.

1.  Change the lines in job.sh according to your needs (frames in your initial training set, path of training directories,...). Further explanations are given in the file job.sh.

2. Prepare a folder with the same name specified in the job.sh file. In this folder, put the training data you will put the training data (energies, forces, dipoles and geometries).

3. In the same folder at point 2, provide the following files: input, active_learning.sh, snap_tr and snap_run.
"inp" is the input file for the quantum chemistry softwaree program that has to be used to perform ab initio calculations. Finally, active_learning.sh is the core driver of the active learning process. It calls "fforge" for the fit, molecular dynamics and calls the external quantum chemistry software for ab initio calculations. This file has to be modified according to your specific calculation. The files snap_tr and snap_run are needed to perform the training and MD during the active learning. Necessary instructions on how to prepare the files are given in "test_active_learning/snap_run" and "test_active_learning/snap_tr".

An example on how to prepare a folder for active learning is given in "test_active_learning".

If you use this program, please cite:
Valerio Briganti and Alessandro Lunghi 2023 Mach. Learn.: Sci. Technol. 4 035005
