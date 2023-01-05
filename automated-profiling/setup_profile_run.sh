#!/bin/bash

# Automates profiling EIRENE
# This script clones EIRENE and sets up the profile cases in source
# Dependencies: jq

auto_prof_dir=$eir_dir/automated-profiling

source $auto_prof_dir/automation_script_header.sh


which jq > /dev/null || (echo "OOPS: Missing dependency: jq"; exit 1;)

#############################################
# Setup top-level directories and git repos #
#############################################

##eir_dir=$EIRDIR

mkdir -p $scalability_report_dir


# Clone or copy EIRENE and build it
# This is not necessary as it stands as this script is within the EIRENE repository.
# However this may change so it will  be left in as it is for now
if [ ! -d $eir_dir ]
then
	echo Cloning eirene into $eir_dir
	git clone $eirene_repo $eir_dir
	cd $eir_dir
	git checkout $eirene_branch
else
	echo Using existing eirene repo in $eir_dir
	cd $eir_dir
fi

git_ref=$(git rev-parse HEAD)

report timestamp "\"`date`\""
report repo \"$eirene_repo\"
report branch \"$eirene_branch\"
report eirene_git_ref \"$git_ref\"
report sample \"$git_ref\"

################
# Build EIRENE #
################

mkdir -p buildRelease
cd buildRelease
echo Building EIRENE in $PWD
module load cmake
FC=mpif90 cmake -DJSON_MODULES=$HOME/lib/json-fortran_gcc/jsonfortran-gnu-8.3.0/lib/ -DLibJSON=$HOME/lib/json-fortran_gcc/jsonfortran-gnu-8.3.0/lib/libjsonfortran.a -DOPENMP=on ../src/
make -j EIRENE

cd $top_dir

# Either clone eirene samples into a local reference repo or copy from existing directory
if [ -d $eirene_samples_dir ]
then
        echo Copying existing eirene samples repo at $eirene_samples_repo into $local_samples_repo
        cp -r $eirene_samples_dir $local_samples_repo
else
	echo Cloning eirene samples into $local_samples_repo
	git clone $eirene_samples_repo $local_samples_repo
fi

cd $local_samples_repo
samples_git_ref=$(git rev-parse HEAD)
report eirene_samples_git_ref \"$samples_git_ref\"

########################
## Generate the cases ##
########################

#Mahti has 256 Cores
max_num_node_threads=256

#This limit is due to the way EIRENE uses output files
max_num_threads=9999

#node_range=( 1 2 4 )
#rank_per_node_range=( 1 2 4 8 16 32 )
#thread_per_rank_range=( 1 2 4 8 16 32 64)
#testing
node_range=( 1 )
rank_per_node_range=( 1 2 4 8 )
thread_per_rank_range=( 1 2 4 8 )

report cases {}
for N in ${node_range[*]}
do
for n_per_N in ${rank_per_node_range[*]}
do
for c in ${thread_per_rank_range[*]}
do
	cd $top_dir
	node_threads=$(($n_per_N * $c))
	total_threads=$(($node_threads * $N))
	[[ $node_threads -gt max_num_node_threads ]] && continue
	[[ $total_threads -gt max_num_threads ]] && continue
	n=$(($n_per_N * $N))
	case_name="Profile_${N}-${n}-${c}"
	echo "--- $case_name ---"
	report cases.${case_name} {}
	report cases.${case_name}.n_nodes $N
	report cases.${case_name}.n_mpi_ranks $n
	report cases.${case_name}.n_omp_threads $c
	if [ ! -d $case_name ]
	then    # Temporary hack for a local repo
	        mkdir $case_name
	        cp -r $local_samples_repo/$sample $case_name/$sample
		#git clone $top_dir/$local_samples_repo $case_name
	        #git clone $local_samples_repo $case_name
	else
		echo $case_name directory already exists...
	fi
	cd $case_name
	cd $sample
	echo Building $case_name/$sample
	make -j
	cd ../..
	echo
done
done
done
