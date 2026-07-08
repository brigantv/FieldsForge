#! /bin/bash

natoms=256

WORK_DIR=/home/users/brigantv/active_learning 		
ABINIT_EXE_DIR=/home/users/brigantv/cp2k-2024.3/exe/local

compound=Mol9_ICMM
project_name_cp2k=mol_9 ### flag specific to CP2K

from_scratch="yes"      ### if set to "yes", it will delete the common record_file where the story of the already finish MD runs is store, "no" if otherwise
consecutive_runs=10	### #number of consecutive runs finished that will lead to the conclusion of this AL session

number_AL=1		### a tracker for parallel sessions of AL
nconfig_in=1		### starting number of configurations in the training set

max_loop=500		### maximum number of MD runs to perform in the AL session
procs=14		### number of processors to run the ab initio calculation


#common directory where you have your common files, be careful to check whether you want to destroy or not the "Record.txt." file
common_dir="${WORK_DIR}/${compound}"

# file_record is the file where you store the last run results (yes/no structure found)
file_record="${common_dir}/Record.txt"

##############################################################################################################################################################

#delete files in common directory according to whether you're restarting everything from scratch or not
if [ ${from_scratch} == "yes" ]
then

rm ${file_record}
rm ${common_dir}/*check

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

mpiexec -n 1 -quiet ./snap snap_tr

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

sed -i '/iseed/c\iseed=['${var1}','${var2}','${var3}','${var4}']' snap_run


echo "Running MD..."

mpiexec -n 1 -quiet ./snap snap_run

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
cp ${common_dir}/*++ ${WORK_DIR}/${compound}_out
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

mpirun -n $procs ${EXE_DIR}/cp2k.popt -i inp -o out_AL_${loop_1}

echo "DFT terminated"

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


# cp2k parsing


force_file="${project_name_cp2k}-force-1_0.xyz"

# write energy 
grep 'ENERGY| Total FORCE_E' out_AL_${loop_1} | tail -n 1 | awk '{printf "%1.12f \n",627.503*$9 }' >> ${FILE_ENE}  #a.u. -> kcal/mol
grep 'ENERGY| Total FORCE_E' out_AL_${loop_1} | tail -n 1 | awk '{printf "%1.12f \n",627.503*$9 }' >> new_ener.txt

# write forces
s=$(( natoms + 5 ))
x=$(awk 'NR>4&&NR<='${s}'{printf "%.17f\n", -$4}' ${force_file} | head -n -1)   # this file needs to be deleted or the next dft calc will add data
y=$(awk 'NR>4&&NR<='${s}'{printf "%.17f\n", -$5}' ${force_file} | head -n -1)
z=$(awk 'NR>4&&NR<='${s}'{printf "%.17f\n", -$6}' ${force_file} | head -n -1)

x=($x)
y=($y)
z=($z)

length_x=${#x[@]}

for (( i=0; i<$length_x; i++ )); do
	echo "${x[i]}" >> new_gradients.txt 
	echo "${x[i]}" >> ${FILE_FORCE}
	echo "${y[i]}" >> new_gradients.txt 
	echo "${y[i]}" >> ${FILE_FORCE}
	echo "${z[i]}" >> new_gradients.txt 
	echo "${z[i]}" >> ${FILE_FORCE}
done
mv ${force_file} ${force_file}_${loop_1}

# write dipoles, to change at a later stage

echo 0.0 >> ${FILE_DIP}
echo 0.0 >> ${FILE_DIP}
echo 0.0 >> ${FILE_DIP}

## end cp2k parsing

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

sed -i 's/nconfig=.*/nconfig='${nconfig}'/g' snap_tr
sed -i 's/nconfig_AL=.*/nconfig_AL=0/g' snap_tr

cp Execution_times ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp AL_stats.txt ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp traj* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp etotal* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp forces_rms* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp ${common_dir}/*++ ${WORK_DIR}/${compound}_out
cp snapcoeff* ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp cumulative_geo_AL.xyz ${WORK_DIR}/${compound}_out/AL_${number_AL}
cp energy_rms* ${WORK_DIR}/${compound}_out/AL_${number_AL}

done
