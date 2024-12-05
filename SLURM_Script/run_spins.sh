#!/bin/bash
#SBATCH --time=47:59:00 # Time limit, just under the 48hr queue limit for default on rocket
#SBATCH --array=1-4%2 # List of jobs to run % Number of concurrent jobs (1-4%2 runs jobs 1-4, 2 at a time)
#SBATCH --ntasks=16 # Number of tasks
#SBATCH --nodes=1  # Force rocket to not split nodes for 2D things
#SBATCH --mail-type=ALL
module load intel/2021a Autoconf/2.71-GCCcore-10.3.0 Automake/1.16.3-GCCcore-10.3.0 libtool/2.4.6-GCCcore-10.3.0

#Set Directory
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" case_list)
cd $DIR

#run the script
srun ./yourcase.x;
# Uncomment the below line if you also want to do a derivatives (requires the derivatives.x and spins.conf_derivs to be present)
# un ./derivatives.x --config=spins.conf_deriv;
#
# Optional steps to move the outputs you want to transfer to local usage into an easy to transfer "transfer" folder:
#mkdir transfer
#mv rho.* ./transfer/;
#mv u.* ./transfer/;
#mv w.* ./transfer/;
#mv vorty.* ./transfer/;
#mv diagnostics.txt ./transfer/;
#mv spins.conf ./transfer/;
#mv stresses_bottom* ./transfer/;
#mv *grid* ./transfer/;
#mv plot* ./transfer/;

cd /nobackup/$USERNAME/SPINS;
