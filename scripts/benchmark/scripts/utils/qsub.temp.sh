#!/bin/bash
# declare a name for this job
#PBS -N hxj5_job
# request the queue for this job
#PBS -q cgsd
# request a total of x processors for this job (y nodes and z processors per node)
#PBS -l nodes=1:ppn=1
# request memory
#PBS -l mem=1gb
# specify walltime
#PBS -l walltime=01:00:00
# out log file
#PBS -o run.out.log
# err log file
#PBS -e run.err.log

# if conda environment is needed
#source /home/xianjie/.bashrc       # this cmd will activate conda base env.
#eval "$(conda shell.bash hook)"    # this will allow run conda inside shell scripts. from https://github.com/conda/conda/issues/7980
#bin_conda=conda
#$bin_conda activate <env>
#$bin_conda deactivate

#change to the directory where you submitted the job
cd $PBS_O_WORKDIR

#include the full path to the name of your program
work_dir=`cd $PBS_O_WORKDIR; pwd`
$work_dir/run.sh

exit 0

# Frequently used PBS commands: (copied from http://physics.princeton.edu/it/faq/pbs-maui-cmds.html) 
#qsub              #submit a job, see man qsub
#qdel -p jobid     #will force purge the job if it is not killed by qdel 
#qstat             #list information about queues and jobs
#showq             #calculated guess which job will run next
#xpbs              #GUI to PBS commands
#qstat -q          #list all queues on system
#qstat -Q          #list queue limits for all queues
#qstat -a          #list all jobs on system
#qstat -s          #list all jobs with status comments
#qstat -r          #list all running jobs
#qstat -f jobid    #list full information known about jobid
#qstat -Qf queueid #list all information known about queueid
#qstat -B          #list summary information about the PBS server
#qstat -iu userid  #get info for queued jobs of userid
#qstat -u userid   #get info for all the jobs of userid
#qstat -n -1 jobid #will list nodes on which jobid is running in one line
#checkjob jobid    #will list job details
