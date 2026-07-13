#! /bin/bash

compound="Dy_b"
WORK_DIR="/home/postdoc3/new/ML_proj/SNAP_PBC"

natoms=115

ABINIT_EXE_DIR=/opt/orca/6.0.0/
compound=Dy_b


from_scratch="yes"
#number of consecutive runs finished that will lead to the conclusion of the job
consecutive_runs=5
number_AL=1
nconfig_in=1
max_loop=100
procs=2

#common directory where you have your common files, be careful to check whether you want to destroy or not the "Record.txt." file
common_dir="${WORK_DIR}/${compound}"

file_record="${common_dir}/Record.txt"

############################################################################

#delete files in common directory according to whether you're restarting everything from scratch or not
if [ ${from_scratch} == "yes" ]
then

yes | rm ${file_record}
yes | rm ${common_dir}/*check

fi

nconfig=$(( ${nconfig_in} ))

###we are assuming that the following files already exist in the common directory, create them in the job file

FILE_ENE=${common_dir}/ener_tr_AL_++
FILE_DIP=${common_dir}/dipoles_tr_AL_++
FILE_GEO=${common_dir}/geo_tr_AL_++
FILE_FORCE=${common_dir}/forces_tr_AL_++
###################


for loop_1 in $(seq 1 ${max_loop} )
do

if [ -e ${file_record} ]
then

check=$( grep "OVER" ${file_record} )

if [ $check ]
then
echo "################## The program is about to be terminated. Active learning completed. ###########################"
exit

fi
fi

echo " Training the force field ${loop_1} "
module load mpi/openmpi-4.1.1
mpiexec -q -n 1 ./fforge snap_tr
#./exe

var1=$(( $RANDOM % 4095 ))
var2=$(( $RANDOM % 4095 ))
var3=$(( $RANDOM % 4095 ))
var4=$(( $RANDOM % 4095 ))

echo ${var1} ${var2} ${var3} ${var4}
var41=$(( ${var4}%2 ))
echo ${var41}

if [ ${var41} == 0 ]
then
var4=$(( ${var4}+1 ))
echo ${var4}
fi

sed -i '/iseed/c\iseed='${var1}','${var2}','${var3}','${var4}'' snap_run


echo "Running MD..."

mpiexec -q -n 1 ./fforge snap_run
echo "Finished MD"

#if the Record file exists,check whether the last consecutive runs have achieved success. If so, end AL and cancel the job tout court.

if [ -e ${file_record} ] 
then

last_10=$( tail -n ${consecutive_runs} ${file_record} | grep '1' | wc -l )

if [ ${last_10} == ${consecutive_runs} ]
then

cp Execution_times ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp AL_stats.txt ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp traj* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp etotal* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp forces_rms* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp snapcoeff* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp cumulative_geo_AL.xyz ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp energy_rms* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp out_AL* ${WORK_DIR}/${compound}_out/AL_${number_AL}

echo "OVER" >> ${file_record}

fi
fi

########################################
# backup of all files at login node for tracking, also here we're assuming that directory AL_${number_AL} exists
cp * ${WORK_DIR}/${compound}_out/AL_${number_AL}

# if the AL has found a new structure, signaled by "GO" in the AL_stats.txt file, launch DFT

last_line=$( tail -n 1 AL_stats.txt | awk '{print $1}' )

if [ $last_line == "GO" ] 
then
echo "Launching DFT..."

module load orca/6.0.0

${ABINIT_EXE_DIR}/orca inp >  out_AL_${loop_1}

echo "DFT terminated"
module unload mpi/openmpi-4.1.1
# adding new structure to cumulative geo AL

cat new_geo_AL.xyz >> cumulative_geo_AL.xyz

echo '0' >> ${file_record}

exec 9>"${FILE_GEO}.lock"
exec 10>"${common_dir}/geo_check.lock"
exec 11>"${FILE_ENE}.lock"
exec 12>"${common_dir}/ene_check.lock"
exec 13>"${FILE_DIP}.lock"
exec 14>"${common_dir}/dip_check.lock"
exec 15>"${FILE_FORCE}.lock"
exec 16>"${common_dir}/force_check.lock"

flock 9 || { echo "Failed to acquire lock."; exit 1; }
flock 10 || { echo "Failed to acquire lock."; exit 1; }
flock 11 || { echo "Failed to acquire lock."; exit 1; }
flock 12 || { echo "Failed to acquire lock."; exit 1; }
flock 13 || { echo "Failed to acquire lock."; exit 1; }
flock 14 || { echo "Failed to acquire lock."; exit 1; }
flock 15 || { echo "Failed to acquire lock."; exit 1; }
flock 16 || { echo "Failed to acquire lock."; exit 1; }


# orca parsing



# write energy 
grep 'FINAL SINGLE*' out_AL_${loop_1} | awk '{print $5}' | awk '{$1=sprintf("%5f",$1*627.509474063)}1'  >> ${FILE_ENE}
grep 'FINAL SINGLE*' out_AL_${loop_1} | awk '{print $5}' | awk '{$1=sprintf("%5f",$1*627.509474063)}1'  >> new_ener.txt


# write dipoles, to change at a later stage
grep 'Total Dipole*' out_AL_${loop_1} | awk '{print $5}' >> ${FILE_DIP}
grep 'Total Dipole*' out_AL_${loop_1} | awk '{print $6}' >> ${FILE_DIP}
grep 'Total Dipole*' out_AL_${loop_1} | awk '{print $7}' >> ${FILE_DIP}

grep 'Total Dipole*' out_AL_${loop_1} | awk '{print $5}' >> new_dipoles.txt
grep 'Total Dipole*' out_AL_${loop_1} | awk '{print $6}' >> new_dipoles.txt
grep 'Total Dipole*' out_AL_${loop_1} | awk '{print $7}' >> new_dipoles.txt

linea=$( grep -n 'gradient' inp.engrad | awk -F : '{print $1}' )

(( q=${linea} + ${natoms}*3 +1 ))
(( k=${natoms}*3 ))
head -n ${q} inp.engrad | tail -n ${k} >> new_gradients.txt
head -n ${q} inp.engrad | tail -n ${k} >> ${FILE_FORCE}

grep 'CENTER*' out_AL_${loop_1} | awk '{print $12,$13,$14,$15}'|sed 's/,//g'|sed 's/[)(]//g' >> ${FILE_SHIFT}

## end orca parsing

echo ${number_AL} >> ${common_dir}/geo_check
echo ${number_AL} >> ${common_dir}/force_check
echo ${number_AL} >> ${common_dir}/ene_check
echo ${number_AL} >> ${common_dir}/dip_check

if [ ${loop_1} == ${max_loop} ]
then
echo "${max_loop} configurations have been added"
fi

exec 9>&-
exec 10>&-
exec 11>&-
exec 12>&-
exec 13>&-
exec 14>&-
exec 15>&-
exec 16>&-

fi

#no matter whether we have found a new structure before or not, we update the number of configurations by giving a look at the "Record.txt"

new_structures=$( grep '0' ${file_record} | wc -l )
(( nconfig=${nconfig_in}+${new_structures} ))

sed -i 's/nconfig.*/nconfig='${nconfig}'/g' snap_tr
sed -i 's/nconfig_AL.*/nconfig_AL='${nconfig}'/g' snap_tr
sed -i 's/nconfig.*/nconfig='${nconfig}'/g' snap_run
sed -i 's/nconfig_AL.*/nconfig_AL='${nconfig}'/g' snap_run

cp Execution_times ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp AL_stats.txt ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp traj* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp etotal* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp forces_rms* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp snapcoeff* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp cumulative_geo_AL.xyz ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp energy_rms* ${WORK_DIR}/${compound}_out/AL_${number_AL}

done
