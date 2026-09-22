# **FieldsForge**

FieldsForge is a software that implements the fit of force fields using the linear model SNAP (Spectral Neighbour Analysis Potential) for both periodic and non-periodic systems. The software can also be used to perform molecular dynamics and geometrical optimization, 
and in the case in which a training set is not available yet, FieldsForge can be driven by an external program present in this repo to perform active learning for the optimal construction of training sets.
For isolated systems, the fit of electric dipoles can also be performed. Its compilation requires the software suite MolForge (https://github.com/LunghiGroup/MolForge) and LAMMPS (https://github.com/lammps/lammps/) tested on version 29 Oct 2020. FieldsForge is known to compile on Linux machines with GNU compilers (v12.2.0) and openmpi(4.1.1).

FieldsForge depends on the linear algebra routines Lapack (v3.8.0) and Blas(v3.8.0).

## **License and Copyright**

FieldsForge: a FORTRAN2003 program for the fit of tensorial machine learning force fields, molecular dynamics and active learning.

Copyright (c) 2026 Università degli Studi di Firenze

Developed by Valerio Briganti, Università degli Studi di Firenze, Dipartimento di Chimica "Ugo Schiff", Via della Lastruccia, 3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

## **Basic instructions**

The program can be compiled with the available Makefile. Two compilation options are available, whether you run it locally or on cluster.

(LOCAL) $ make 
(CLUSTER) $ make PLATFORM=hpc

Check Makefile to change the relevant paths.

The software is executed using:

$ ./fforge snap.input

"snap.input" is the text file containing the instructions for the program. In the repo is provided a template with explanations of the meaning of the single variables as comments.

To train and validate a force field, please go the "test_benzene" folder.
To perform an active learning cycle, please go the "test_active_learning" folder.

## **How to cite FieldsForge**

When you use this program, please cite:
- Valerio Briganti and Alessandro Lunghi 2023 Mach. Learn.: Sci. Technol. 4 035005
- all the relevant papers for MolForge and LAMMPS
