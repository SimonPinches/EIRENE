The shell scripts were originally taken from Oskar Lappi at
https://version.helsinki.fi/lapposka/eirene-omp-profiling-2021
and altered for this project.

Requirements:
python
jq

The set consists of 5 scripts. They can be run individually or collectively using run_profile.sh

1) setup_profiling_run.sh - Build Eirene, clone sample cases, create a directory for each combination of nodes, tasks and threads

2) launch_jobs.sh - Launch a slurm job for each profile case

3) collect_profile_data.sh - Collect timing data from each run and add to a json report file

4) plot_scaling.py - Plot png figures of strong scaling, efficiency and speedup

5) automation_script_header - Contains parameter settings and reporting functions

6) run_profiling.sh - Runs all steps required to produce scaling plots


The following parameters in the must be set for this to work
eir_dir - the location of the eirene repository, either in run_profiling or the runtime environment

The automation_script_header.sh script contains parameters for slurm job submission
account_name
partition_name
It may also be necessary to set a QoS value

A location of an existing eirene samples repository can also be set in automation_script_header.sh

The range of scaling can also be set in setup_profile_run.sh