#!/bin/bash

# Automates profiling EIRENE
# This script collects data from finished EIRENE profile runs
# Dependencies: jq

if [ ! -v eir_dir ]
then
    echo Error: eir_dir not set. Exiting
    exit -1
fi

auto_prof_dir=$eir_dir/automated-profiling
source $auto_prof_dir/automation_script_header.sh

if [ ! -d $auto_prof_dir ]
then
    echo Error: $auto_prof_dir not found. Exiting
    exit -1
fi

which jq > /dev/null || (echo "OOPS: Missing dependency: jq"; exit 1;)

echo Collecting profile data

cases=$(jq '.cases | keys[]' $report_file | tr -d '"') 
for case_name in $cases
do
	echo 
        echo "---$case_name---"
	cd $top_dir

	# Check launch success
	launch_success=$(read_report cases.${case_name}.launch_success)
	if [ "$launch_success" == "false" ];then
		echo "Case failed to launch"
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
		if [ "$job_state" == "FAILED" ]
		then
		        echo "Job ${job_number} failed, job state is ${job_state}"
			break
		fi
		if [ "$job_state" == "null" ]
		then
		        echo "Invalid job number"
			break
		fi
		echo "Job ${job_number} not complete yet, job state is ${job_state} ... Waiting ..."
		sleep 10
		job_state=$(sacct -j $job_number --format jobid,state -n |sed -n -r "/^$job_number /s/($job_number| +)//gp")
	done
	report cases.${case_name}.job_state \"$job_state\"
	
	#Make a directory in the report dir for the case
	case_report_dir=$scalability_report_dir/${sample}_profiles/$case_name
	echo "creating $case_report_dir"
	mkdir -p $case_report_dir

	####################
	# Copy output      #
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

	if [ "$job_state" == "FAILED" ]
	then
	        echo "Job ${job_number} failed, skipping parsing"
	else
    	        echo "Parsing output"
		if [ $N -gt 1 ] || [ $n -gt 1 ] || [ $c -gt 1 ]; then
   		        output=output.0000
		else
		        output=eirene-2d.reference_${N}-${n}-${c}.out
		fi
		report cases.${case_name}.timing {}
		cpu_time=$(grep CPU_TIME ${case_report_dir}/eirene-2d.reference_${N}-${n}-${c}.out | sed 's/[[:blank:]]*$//; s/.*[[:blank:]]//')
		wall_time=$(echo ${cpu_time:0:8})
		atom_time=$(grep 'TOTAL=' ${output} | sed '1q;d' |  sed 's/[[:blank:]]*$//; s/.*[[:blank:]]//')
		mol_time=$(grep 'TOTAL=' ${output} | sed '2q;d' |  sed 's/[[:blank:]]*$//; s/.*[[:blank:]]//')
		ion_time=$(grep 'TOTAL=' ${output} | sed '3q;d' |  sed 's/[[:blank:]]*$//; s/.*[[:blank:]]//')
		report cases.${case_name}.output_parsed "true"
		echo "Adding wall time: ${wall_time}"
		report cases.${case_name}.timing.wall_time ${wall_time}
		report cases.${case_name}.timing.atom_time ${atom_time}
		report cases.${case_name}.timing.mol_time ${mol_time}
		report cases.${case_name}.timing.ion_time ${ion_time}
	fi
done
