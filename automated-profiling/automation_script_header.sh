# Loccation of the eirene repository
# This isn't necessary as the automated scripts are in the eirene repository,
# however it may be desirable to build and use a different version so this is left in
# in case that functionality is desired. possible TODO: alter scripts to allow this
# Juelich repo
#eirene_repo=git@jugit.fz-juelich.de:eirene/eirene.git

# Choose either the eirene samples juelich repository or local cases
# Juelich
eirene_samples_repo=git@jugit.fz-juelich.de:eirene/EIRENE-sample-cases.git
eirene_samples_branch=reduced_profiling_cases

# Uncomment the following line and replace with the full path to a local samples directory to avoid
# cloning from the Juelich repository which can take a long time and may require a password depending
# on your setup. If you want to use the automatically cloned version but stop it cloning set the line
# below to its location and if using run_profiling comment out the rm -rf EIRENE_SAMPLES_LOCAL line.
#export eirene_samples_dir=$HOME/f-eirene/EIRENE-sample-cases

local_samples_repo=EIRENE_SAMPLES_LOCAL

# Set the sample case to be used, at present this only works for the 2D-D_slab case
sample=2D-D_slab

# Submission script parameters
# Variables for the job submission system
# You can also alter the submission script directly
# Change this to a valid account name (no spaces)
export account_name=project_2004753
# and partition
export partition_name=medium
# If it is necessary to set a QoS value uncomment this line and set it
#qos_name=normal


top_dir=$PWD
scalability_report_dir=$top_dir/SCALABILITY_REPORT
report_filename=report_${sample}.json
report_file=$scalability_report_dir/$report_filename

ensure_report_exists() {
	if [ ! -f $report_file ]
	then
		touch $report_file
	fi
	if [ -z "$(jq "." < "$report_file")" ]
	then
		echo "{}" > $report_file
	fi
}

#Report writes to a json file that acts as a data store
report(){
	ensure_report_exists
	if [ $# -ne 2 ]
	then 
		echo "OOPS: expected 2 arguments to report(), got $#: $@"
		exit 10
	fi
	#Rewrite the key into a form that jq understands
	#This way we don't need to worry about quoting while calling report
	raw_key=$1
	IFS=. read -a key_path <<< $raw_key
	key=""
	for key_part in ${key_path[*]}
	do
		key+=.\"$key_part\"
	done
	tmpfile=$(mktemp)
	data=$(<$report_file)
	val=$2
	jq "$key=$val" <<< $data > $tmpfile || (echo "OOPS: jq can't write to file"; exit 11;)
	mv $report_file $report_file.bak
	mv $tmpfile $report_file
}

read_report(){
	ensure_report_exists
	if [ $# -ne 1 ]
	then
		echo "OOPS: expected 1 argument to read_report(), got $#: $@"
		exit 20
	fi
	raw_key=$1
	IFS=. read -a key_path <<< $raw_key
	key=""
	for key_part in ${key_path[*]}
	do
		key+=.\"$key_part\"
	done

	data=$(<$report_file)
	jq "$key" <<< $data || (echo "OOPS: jq can't read from file"; exit 21;)
}
