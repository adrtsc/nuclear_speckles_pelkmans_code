import sys
import os


SLURM_COMMAND = """#! /bin/sh
#SBATCH --array=0-{0}
#SBATCH -o /home/atschan/PhD/slurm_reports/slurm-%A_%a.out
#SBATCH -e /home/atschan/PhD/slurm_reports/slurmerror-%A_%a.out
#SBATCH --mem-per-cpu=8000m
#SBATCH --cpus-per-task=1
#SBATCH --time=240

n="$SLURM_ARRAY_TASK_ID"

exec python get_transcript_features_20230724.py {0} $n
"""

# load settings

n_jobs = 80

with open("temp.sh", "w") as f:
    f.write(SLURM_COMMAND.format(n_jobs))
os.system("sbatch temp.sh")
os.unlink("temp.sh")
