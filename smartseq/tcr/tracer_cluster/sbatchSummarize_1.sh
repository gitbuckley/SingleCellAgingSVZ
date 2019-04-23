#!/bin/bash
#SBATCH --account=abrunet1
#SBATCH --job-name="tracerConda"
#SBATCH --nodes 1 # Number of tasks or threads
#SBATCH -t 1-00:01 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=16G # see also --mem 
#SBATCH -o out.%j # File to which STDOUT will be written
#SBATCH -e err.%j # File to which STDERR will be written
#SBATCH --mail-user=buckley7@stanford.edu
#SBATCH --mail-type=END,FAIL # What type of mail to send

# Purpose: Summarize "tracer assemble" job results across several cells.
date

module load anaconda
source activate tracer_teichlab_2018-03-08

cd /srv/gsfs0/projects/brunet/Buckley/5.BenNSCProject/2.Tcell/

tracer summarise -c 5.Tracer/CONFIG \
				--resource_dir 5.Tracer/resources \
				5.Tracer/Output_Conda_Mar15

date
echo "Finished"
