#!/bin/bash

# Automates profiling EIRENE
# This script launches jobs, it may have to be run many times
# Dependencies: jq

auto_prof_dir=$eir_dir/automated-profiling
source $auto_prof_dir/automation_script_header.sh

launch_profile() {
	N=$1
	n=$2
	c=$3
sbatch 2>&1 << EOF
#!/bin/bash
#SBATCH --job-name=eirprof_${N}-${n}-${c}
#SBATCH --account=project_2004753
#SBATCH --time=00:30:00
#SBATCH --partition=medium
#SBATCH -N $N
#SBATCH -n $n
#SBATCH -c $c

export OMP_NUM_THREADS=$c
./set_links.sh
srun ./eirene 1>eirene-2d.reference.out 2>eirene-2d.reference.err
cp eirene-2d.reference.out eirene-2d.reference_${N}-${n}-${c}.out
cp fort.200 eirene-2d.reference_${N}-${n}-${c}.out
cp fort.7 eirene-2d.reference_${N}-${n}-${c}.out
./rm_links.sh
EOF
	return $?
}


which jq > /dev/null || (echo "OOPS: Missing dependency: jq"; exit 1;)
#Iterate through cases
cases=$(jq '.cases | keys[]' $report_file | tr -d '"') 
for case_name in $cases
do
	echo
	echo "---$case_name---"
	#Check if it launched
	launch_success=$(read_report cases.${case_name}.launch_success)
	if [ "$launch_success" == "true" ];then
		echo "Already launched "
		continue
	fi
	
	#Check if configuration is bad
	bad_config=$(read_report cases.${case_name}.bad_config)
	if [ "$launch_success" == "true" ];then
		echo "Bad node configuration"
		continue
	fi
	
	####################
	# Launch SLURM job #
	####################
	
	N=$(read_report cases.${case_name}.n_nodes)
	n=$(read_report cases.${case_name}.n_mpi_ranks)
	c=$(read_report cases.${case_name}.n_omp_threads)
	printf "Attempting to launch ${case_name}/${sample}... "
	cd $top_dir/$case_name/$sample
	sbatch_return=$(launch_profile $N $n $c)
	if [[ $? -eq 0 ]];then
		report cases.${case_name}.launch_success "true"
		job_number=$(cut -d" " -f4 <<< $sbatch_return)
		report cases.${case_name}.job_number \"$job_number\"
		echo "Successfully launched sbatch job ${job_number}"
		echo $sbatch_return
	else
		echo "Failed to launch :(, exit code $? "
		report cases.${case_name}.launch_success "false"
		bad_config_msg="Requested node configuration is not available"
		if [[ $sbatch_return =~ $bad_config_msg ]]; then
			report cases.${case_name}.bad_config "true"
			echo "Bad node configuration"
		else
			echo $sbatch_return
			echo "Exiting"
			exit 100
		fi
	fi
	echo
done
