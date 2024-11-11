import sys
import os


SLURM_COMMAND = """#! /bin/sh
#SBATCH --array=0-{0}
#SBATCH -o /home/atschan/PhD/slurm_reports/slurm-%A_%a.out
#SBATCH -e /home/atschan/PhD/slurm_reports/slurmerror-%A_%a.out
#SBATCH --mem-per-cpu=60000m
#SBATCH --cpus-per-task=1
#SBATCH --time=240

n="$SLURM_ARRAY_TASK_ID"

exec python 02_MLR_combinatorial.py {0} $n
"""

# load settings

n_jobs = 0

with open("temp.sh", "w") as f:
    f.write(SLURM_COMMAND.format(n_jobs))
os.system("sbatch temp.sh")
os.unlink("temp.sh")
