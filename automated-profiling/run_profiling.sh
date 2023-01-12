#!/bin/bash

##############################################
# Setup
##############################################

# The location of an existing eirene installation
eir_dir=$HOME/f-eirene/eirene
export eir_dir

# Remove existing profile data
rm -rf SCALABILITY_REPORT
rm -rf Profile*
# Comment this line if you want to keep an already cloned local version
rm -rf EIRENE_SAMPLES_LOCAL

# And Uncomment the following line and replace with the full path to a local samples directory to avoid cloning from the Juelich repository
# which can take a long time and may require a password depending on your setup
eirene_samples_dir=$HOME/f-eirene/EIRENE-sample-cases
#eirene_samples_dir=$HOME/f-eirene/scratch/eirene-automated-profiling/EIRENE_SAMPLES_LOCAL
export eirene_samples_dir

###############################################
# Calling profiling scripts
###############################################

# Set up profile directories
$eir_dir/automated-profiling/setup_profile_run.sh
if [ $? -ne 0 ]; then
        exit -1
fi

# Launch jobs
$eir_dir/automated-profiling/launch_jobs.sh
if [ $? -ne 0 ]; then
        exit -1
fi

# Collect the profile data
$eir_dir/automated-profiling/collect_profile_data.sh
if [ $? -ne 0 ]; then
        exit -1
fi

# Plot the profile data using python
$eir_dir/automated-profiling/plot_scaling.py
if [ $? -ne 0 ]; then
        exit -1
fi
