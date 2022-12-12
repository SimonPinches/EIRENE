eirene_repo=git@jugit.fz-juelich.de:eirene/eirene.git
#eirene_samples_repo=git@jugit.fz-juelich.de:eirene/EIRENE-sample-cases.git
eirene_samples_repo=$HOME/f-eirene/EIRENE-sample-cases

eir_dir=eirene
eirene_branch=develop_openmp

local_samples_repo=EIRENE_SAMPLES_LOCAL

sample=2D-D_slab

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
