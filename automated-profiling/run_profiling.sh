#!/bin/bash
rm -rf SCALABILITY_REPORT EIRENE_SAMPLES_LOCAL Profile*

export eir_dir=$HOME/f-eirene/eirene
$eir_dir/automated-profiling/setup_profile_run.sh
