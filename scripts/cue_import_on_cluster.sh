#!/bin/bash
#SBATCH --job-name=muscle_tissue # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=32       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --time=00:60:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email when job fails
#SBATCH --mail-user=cm7897@princeton.edu

module purge
module load anaconda3/2023.9
conda activate env_maldi


# CHOOSE THE SCIRPT YOU NEED TO RUN BELOW
# =======================================

python import_mat_substract_background.py
#python import_mat.py

# =======================================
