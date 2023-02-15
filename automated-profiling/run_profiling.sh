#!/bin/bash

##############################################
# Setup
##############################################

# The location of an existing eirene installation
export eir_dir=$HOME/f-eirene/eirene

# Remove existing profile data
rm -rf SCALABILITY_REPORT
rm -rf Profile*
# Comment this line if you want to keep an already cloned local version in the rundir
rm -rf EIRENE_SAMPLES_LOCAL
# You will also need need to set the location of the samples repository in automation_script.sh

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
