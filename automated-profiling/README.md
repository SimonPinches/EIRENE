This profiling setup worked in July/August of 2021, if EIRENE's build system has been stable since it should still work

The profiling was done on Mahti

Step 0. Set parameters for profiling run by editing `setup_profile_run.sh`

Step 1. `./setup_profile_run.sh`, clone & build EIRENE, generate cases

Step 2. `./launch_jobs.sh`,  launch SLURM jobs for all the cases

Step 3. `./collect_profile_data.sh`, collect the generated gprof data
