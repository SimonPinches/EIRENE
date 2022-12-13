#!/bin/bash

# Automates profiling EIRENE
# This script collects data from finished EIRENE profile runs
# Dependencies: jq

auto_prof_dir=$eir_dir/automated-profiling
source $auto_prof_dir/automation_script_header.sh

which jq > /dev/null || (echo "OOPS: Missing dependency: jq"; exit 1;)

cases=$(jq '.cases | keys[]' $report_file | tr -d '"') 
for case_name in $cases
do
	echo 
        echo "---$case_name---"
	cd $top_dir

	profile_generated=$(read_report cases.${case_name}.profile_success)
	if [ "$profile_copied" == "true" ];then
		echo "Already profiled"
		continue
	fi	

        #Check if we've already parsed this case
        profile_parsed=$(read_report cases.${case_name}.profile_parsed)
        if [ "$profile_parsed" == "true" ];then
                echo "Already parsed"
                continue
        fi      
	   
	#Check to see if the slurm job has finished
	job_number=$(read_report cases.${case_name}.job_number |tr -d '"')
	job_state=$(sacct -j $job_number --format jobid,state -n |sed -n -r "/^$job_number /s/($job_number| +)//gp")
	while [ "$job_state" != "COMPLETED" ]; do	        
		job_state=$(sacct -j $job_number --format jobid,state -n |sed -n -r "/^$job_number /s/($job_number| +)//gp")
		echo "Job ${job_number} not complete yet, job state is ${job_state} ... Waiting ..."
		sleep 10
	done
	report cases.${case_name}.job_state \"$job_state\"

	#Make a directory in the report dir for the case
	case_report_dir=$scalability_report_dir/${sample}_profiles/$case_name
	mkdir -p $case_report_dir

	####################
	# Generate profile #
	####################

	N=$(read_report cases.${case_name}.n_nodes)
	n=$(read_report cases.${case_name}.n_mpi_ranks)
	c=$(read_report cases.${case_name}.n_omp_threads)
	cd $case_name/$sample
	echo "Generating profile for slurm job ${job_number}"
	echo "Copying job artifacts"
	cp eirene-2d.reference_${N}-${n}-${c}.out slurm-${job_number}.out $case_report_dir/
	echo "Copied output ${sample}_profiles/$case_name/profile-${N}-${n}-${c}"
	report cases.${case_name}.output_copied "true"

	echo "Parsing output"
	report cases.${case_name}.timing {}
	cpu_time=$(grep CPU_TIME ${case_report_dir}/eirene-2d.reference_${N}-${n}-${c}.out | sed 's/[[:blank:]]*$//; s/.*[[:blank:]]//')
	wall_time=$(echo ${cpu_time:0:8})
	report cases.${case_name}.output_parsed "true"
	report cases.${case_name}.timing.wall_time \"${wall_time}\"
done
