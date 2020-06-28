#!/bin/bash

# See here for more info on GIZMO: https://teams.fhcrc.org/sites/citwiki/SciComp/Pages/Gizmo%20Cluster%20Quickstart.aspx

# 1. Transfer files to cluster
#     a. Should be a directory (e.g. 'z.stepped.wedge/'), containing an R/ directory, containing a file called MAIN.R or MAIN.Rmd

# 2. Call this bash script (run_r.sh) from the console via one of the following:
#   BOTH
#     Variable `cluster` can be either 'gizmo' or 'bayes'
#     Variable `type` can be either 'R' or 'Rmd'
#     Variable `project` should be the name of the directory (e.g. 'z.stepped.wedge')
#     Variable `par_type` ("parallelization type") can be "inner", "outer", or "none" (!!!!! expand documentation of these options)
#   GIZMO
#     sbatch --export=cluster='gizmo',type='R',project='z.stepped.wedge',add_to_tid=0 -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out run_r.sh
#       add -t 7-0 (days-hours) for jobs that will take longer than 72 hours
#       add --partition=largenode to use large nodes (up to 24 cores)
#       add --array=1-10 to run run_r.R ten times (in a job array)
#       add -c 4 to request four cores
#   BAYES
#     cd Desktop
#     qsub -v cluster='bayes',type='R',project='z.stepped.wedge' -cwd -e ./io/ -o ./io/ run_r.sh
#       add -t 1-10 to run run_r.R ten times (in a job array)

# 3. Check on jobs, delete jobs, see old jobs
#   GIZMO
#     a. squeue -u akenny (checks on all jobs)
#     b. scancel 123 (deletes job 123)
#     c. sacct -u akenny --format=JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,nnodes,ncpus
#   BAYES
#     a. qstat (checks on all jobs)
#     b. qdel 123 (deletes job 123)



# Gizmo code
if [ "$cluster" == "gizmo" ]; then

  # Load R module
  # !!!!! Wrap this in an if block that only loads modules if they are not already loaded
  module use /app/easybuild/modules/all
  module load R/3.6.2-foss-2019b-fh1
#  module load JAGS/4.2.0-foss-2016b # !!!!! Make this a command line argument
  
fi



# If needed, add some integer to .tid (to deal with job array size limits)
if [ -z "$add_to_tid" ]; then
  add_to_tid=0
fi



# Run script
Rscript run_r.R -cwd --cluster $cluster --project $project --type $type --add_to_tid $add_to_tid
