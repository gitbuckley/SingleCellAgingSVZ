# Purpose: Submit Tracer jobs to cluster. 1 job per cell.

# Create list of unique truncated file names from a given directory
from os import listdir
from os.path import isfile, join
from subprocess import call
import re

mypath = "/srv/gsfs0/projects/brunet/Buckley/5.BenNSCProject/2.Tcell/0.Raw"
print(mypath)

# Get all files in current directory (exlude directories)
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

# Filter for file names that contain a certain string
onlyfqfiles = [f for f in onlyfiles if re.search("Young|Old", f)]

# Trim off the last few characters of each files name.
trimmedfiles = list(sorted(set([name[0:-13] for name in onlyfqfiles])))
print("Number of files:  " + str(len(trimmedfiles)))
print("Example file:  " + trimmedfiles[1])

# Submit JOBS.
for file in trimmedfiles:
    cmd = "sbatch --export=F='{0}' --job-name='{0}' sbatchCondaTracer.sh".format(file)
    print(cmd)
#    exit_status = call(cmd, shell=True)
 #   if exit_status is 1:
  #      print("Job:  {} -----> failed to submit.").format(cmd)
