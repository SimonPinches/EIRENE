#!/bin/bash
# The location of an existing eirene installation
eir_dir=$HOME/f-eirene/eirene

export eir_dir

# Remove existing profile data
rm -rf SCALABILITY_REPORT EIRENE_SAMPLES_LOCAL Profile*

# Set up profile directories
$eir_dir/automated-profiling/setup_profile_run.sh

# Launch jobs
$eir_dir/automated-profiling/launch_jobs.sh

# Collect the profile data
$eir_dir/automated-profiling/collect_profile_data.sh

# Plot the profile data using python
$eir_dir/automated-profiling/plot_scaling.py
