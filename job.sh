#!/bin/bash

#SBATCH -N 1      # Request 1 node
#SBATCH --ntasks-per-node=6
#SBATCH -t 0-48:00:00      # Request 12 hours and 0 minutes
#SBATCH --nodelist=node8
#SBATCH -p bergamo         # Use the "compute" partition
#SBATCH -J MD     # Job name

#modules to upload
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/postdoc6/lammps/examples/COUPLE/fortran2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/postdoc6/lammps/src
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/postdoc6/MolForge/libs
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/postdoc6/dftd3-lib-0.9/lib

export OMP_NUM_THREADS=2
export OMPI_MCA_btl='^uct,ofi'
export OMPI_MCA_pml='ucx'
export OMPI_MCA_mtl='^ofi'

#########################################	 USER SECTION (MODIFY ONLY THIS PART)

export WORK_DIR=/home/postdoc3/SNAP_PBC		### directory where you put the subdir related to the specific compound under study
export SNAP_EXE_DIR=/home/postdoc3/SNAP_PBC	### directory where the executable of snap is located (exe_name: snap)
compound="test_active_learning"					### subdir in WORK_DIR where all the input data is located

natoms=12							### number of atoms in the system
ABINIT_EXE_DIR=/opt/orca/6.0.0					### directory containing the ab initio  executable
from_scratch="yes"						### if set to "yes", deletes the common Record.txt file, restarting the AL session from scratch
consecutive_runs=5					### number of consecutive successful runs that will lead to the conclusion of the AL session
max_loop=500							### maximum number of MD runs to perform in the AL session
nconfig=1							###initial configurations in the training set
init_AL=1							### start numbering of AL directories on compute node (important if you already have AL running on nodes)
tot_AL=0							### how many parallel AL sessions to run

FILE_ENE_TO_COPY=ener_tr				###these files contain your starting data
FILE_DIP_TO_COPY=dipoles_tr
FILE_GEO_TO_COPY=geo_tr
FILE_FORCE_TO_COPY=forces_tr


##########################################

(( num_AL=${init_AL}+${tot_AL} ))

nodes=($(scontrol show hostnames $SLURM_NODELIST))
# Loop to create directories and prepare input data on each node
for i in $(seq ${init_AL} ${num_AL} ); do
(( k=${i}-1 ))
node_name="${nodes[$k]}"
   # Dynamically assign TMP_DIR for each node 
   TMP_DIR_NODE="/tmp/${SLURM_JOB_ID}/AL_${i}"  # Unique TMP_DIR per node/tas
   echo "Sending command to $node_name"

   mkdir -p ${WORK_DIR}/${compound}_out/AL_${i}
   # Use srun or ssh to run a command on the node
   srun -N1 -n1 --nodelist=$node_name bash -c "
   mkdir -p ${TMP_DIR_NODE}
   rm -f ${TMP_DIR_NODE}/*
   
    cp ${WORK_DIR}/${compound}/snap_tr   ${TMP_DIR_NODE}     
    cp ${WORK_DIR}/${compound}/snap_run	 ${TMP_DIR_NODE}
    cp ${SNAP_EXE_DIR}/fforge ${TMP_DIR_NODE}

    rm -f ${TMP_DIR_NODE}/traj_MD*
    rm -f ${TMP_DIR_NODE}/etotal*

    cp ${WORK_DIR}/${compound}/inp ${TMP_DIR_NODE}
    cp ${WORK_DIR}/${compound}/*gbw ${TMP_DIR_NODE}
    cp ${WORK_DIR}/${compound}/active_learning.sh ${TMP_DIR_NODE}
    "

    FILE_ENE=ener_tr_AL_++
    FILE_DIP=dipoles_tr_AL_++
    FILE_GEO=geo_tr_AL_++
    FILE_FORCE=forces_tr_AL_++

    # Copy the initial data to the common directory for multiple trajectories
    cp ${WORK_DIR}/${compound}/${FILE_ENE_TO_COPY} ${WORK_DIR}/${compound}/${FILE_ENE}
    cp ${WORK_DIR}/${compound}/${FILE_DIP_TO_COPY} ${WORK_DIR}/${compound}/${FILE_DIP}
    cp ${WORK_DIR}/${compound}/${FILE_GEO_TO_COPY} ${WORK_DIR}/${compound}/${FILE_GEO}
    cp ${WORK_DIR}/${compound}/${FILE_FORCE_TO_COPY} ${WORK_DIR}/${compound}/${FILE_FORCE}

done


# Loop 2: Execute tasks in the background
for i in $(seq ${init_AL} ${num_AL} ); do
    TMP_DIR_NODE="/tmp/${SLURM_JOB_ID}/AL_${i}"

    # Run in a bash subshell in the background (replacing srun)
    (
        cd ${TMP_DIR_NODE}
    sed -i 's/natoms=.*/natoms='"${natoms}"'/g' active_learning.sh
    sed -i 's#WORK_DIR=.*#WORK_DIR='"${WORK_DIR}"'#g' active_learning.sh
    sed -i 's#ABINIT_EXE_DIR=.*#ABINIT_EXE_DIR='"${ABINIT_EXE_DIR}"'#g' active_learning.sh
    sed -i 's/compound=.*/compound='"${compound}"'/g' active_learning.sh
    sed -i 's/project_name_cp2k=.*/project_name_cp2k='"${project_name_cp2k}"'/g' active_learning.sh
    sed -i 's/from_scratch=.*/from_scratch=\"'"${from_scratch}"'\"/g' active_learning.sh
    sed -i 's/consecutive_runs=.*/consecutive_runs='"${consecutive_runs}"'/g' active_learning.sh
    sed -i 's/number_AL=.*/number_AL='"${i}"'/g' active_learning.sh
    sed -i 's/nconfig_in=.*/nconfig_in='"${nconfig}"'/g' active_learning.sh
    sed -i 's/max_loop=.*/max_loop='"${max_loop}"'/g' active_learning.sh
    sed -i 's/procs=.*/procs='"${procs}"'/g' active_learning.sh

    # Fill in the Training data block of the namelist input files with the paths of the common files
# Fill in the Training data block of the namelist input files with the paths of the common files
    sed -i "s#^\( *nconfig *= *\).*#\1${nconfig}#" snap_tr
    sed -i "s#^\( *geometry_file *= *\).*#\1\"${WORK_DIR}/${compound}/${FILE_GEO}\"#" snap_tr
    sed -i "s#^\( *energy_file *= *\).*#\1\"${WORK_DIR}/${compound}/${FILE_ENE}\"#" snap_tr
    sed -i "s#^\( *forces_file *= *\).*#\1\"${WORK_DIR}/${compound}/${FILE_FORCE}\"#" snap_tr
    sed -i "s#^\( *dipoles_file *= *\).*#\1\"${WORK_DIR}/${compound}/${FILE_DIP}\"#" snap_tr
    sed -i "s#^\( *record_file *= *\).*#\1\"${WORK_DIR}/${compound}/Record.txt\"#" snap_tr

    sed -i "s#^\( *nconfig *= *\).*#\1${nconfig}#" snap_run
    sed -i "s#^\( *geometry_file *= *\).*#\1\"${WORK_DIR}/${compound}/${FILE_GEO}\"#" snap_run
    sed -i "s#^\( *energy_file *= *\).*#\1\"${WORK_DIR}/${compound}/${FILE_ENE}\"#" snap_run
    sed -i "s#^\( *forces_file *= *\).*#\1\"${WORK_DIR}/${compound}/${FILE_FORCE}\"#" snap_run
    sed -i "s#^\( *dipoles_file *= *\).*#\1\"${WORK_DIR}/${compound}/${FILE_DIP}\"#" snap_run
    sed -i "s#^\( *record_file *= *\).*#\1\"${WORK_DIR}/${compound}/Record.txt\"#" snap_run
        # Run the active learning script natively
        ./active_learning.sh
    ) &
done

wait

