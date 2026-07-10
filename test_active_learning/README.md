Here an example on how to run active learning for benzene starting from one configuration in the training set.

1. In the folder, where you will sbatch job.sh, prepare a subfolder "test_active_learning" where you will put the following files: energies, geometries, forces and dipoles, inp, snap_run and snap_tr, and active_learning.sh.

inp: input file for ORCA
snap_tr: training instructions for "snap"
snap_run: instructions for MD and active learning for "snap".
active_learning.sh: file containing loading of modules relevant for the cluster and parsing of the outpout files of ORCA.

2. Prepare job.sh with the relevant information for the system. They are all contained in the "USER SECTION".

3. Launch job.sh. All the results of the fit and MD are constantly updated in the "test_active_learning_out" folder. 

No names of the subdirectories or training data files are hardocoded, so they can be changed accordingly in job.sh.
