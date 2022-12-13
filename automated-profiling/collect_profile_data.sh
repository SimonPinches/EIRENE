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

	#Check if we've already processed this case
	profile_generated=$(read_report cases.${case_name}.profile_success)
	if [ "$profile_generated" == "true" ];then
		echo "Already profiled"
		continue
	fi	

	#Check to see if the slurm job has finished
	job_number=$(read_report cases.${case_name}.job_number |tr -d '"')
	job_state=$(sacct -j $job_number --format jobid,state -n |sed -n -r "/^$job_number /s/($job_number| +)//gp")
	report cases.${case_name}.job_state \"$job_state\"

	while [ "$job_state" != "COMPLETED" ]; do	        
		echo "Job ${job_number} not complete yet, job state is ${job_state} ... Waiting ..."
		sleep 10
		job_number=$(read_report cases.${case_name}.job_number |tr -d '"')
		job_state=$(sacct -j $job_number --format jobid,state -n |sed -n -r "/^$job_number /s/($job_number| +)//gp")
		report cases.${case_name}.job_state \"$job_state\"
	done

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
#	printf "Parsing gprof..."
#	gprof ./eirene gmon_${N}-${n}-${c}.out > $case_report_dir/profile-${N}-${n}-${c}

#	if [[ $? -eq 0 ]];then
#		echo "Successfully generated profile ${sample}_profiles/$case_name/profile-${N}-${n}-${c}"
#		report cases.${case_name}.profile_success "true"
#	else
#		echo "Failed to generate profile :("
#		report cases.${case_name}.profile_success "false"
#	fi
done

#TODO: this loop is grafter on here because I've ran the upper loop separately, they could be combined
for case_name in $cases
do
	echo 
        echo "---$case_name---"
	cd $top_dir

	#Check if we've already processed this case
	profile_parsed=$(read_report cases.${case_name}.profile_parsed)
	if [ "$profile_parsed" == "true" ];then
		echo "Already parsed"
		continue
	fi	

	#Check if we've generated the profile
	profile_generated=$(read_report cases.${case_name}.profile_success)
	if [ "$profile_generated" != "true" ];then
		echo "Not profiled"
		continue
	fi	

	#Make a directory in the report dir for the case
	case_report_dir=$scalability_report_dir/${sample}_profiles/$case_name
	mkdir -p $case_report_dir
	cd $case_report_dir
	
	############################
	# Parse gprof flat profile #
	############################

	report cases.${case_name}.flat_profile {}

	N=$(read_report cases.${case_name}.n_nodes)
	n=$(read_report cases.${case_name}.n_mpi_ranks)
	c=$(read_report cases.${case_name}.n_omp_threads)

	echo "Parsing output"
#	sed -n -r '/^ *[0-9]/{p; :loop n; /^$/q; p; b loop}' ${case_report_dir}/profile-${N}-${n}-${c} | \
#	while read pcent cumul self n_calls self_ms total_ms name
#	do
#		#Sometimes the last three numerical values are missing, so the name field is shifted "left"
#		[ -z "$name" ] && name=$total_ms
#		[ -z "$name" ] && name=$self_ms
#		[ -z "$name" ] && name=$n_calls;n_calls=""
#		if [[ "$self" != "0.00" ]]; then
#			report cases.${case_name}.flat_profile.${name} {}
#			report cases.${case_name}.flat_profile.${name}.time \"${self}\"
#			report cases.${case_name}.flat_profile.${name}.percentage \"${pcent}\"
#			report cases.${case_name}.flat_profile.${name}.n_calls \"${n_calls}\"
#		fi
#	done 

done
