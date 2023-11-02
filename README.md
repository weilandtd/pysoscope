# pysoscope

## Convert maldi files to numpy


# Running the scripts on the cluster 
To run on the cluster first make the conda environment:

```
conda create --name env_maldi --file reqrequirements.txt
```

Modify the SLURM script we recommend using `nano` or `vim`. 
```
cd /path/to/pysoscope
cd ./scripts
nano cue_import_on_cluster.sh
```
Make sure you change your email adress, the wall time and the scirpt you want to run:
```
...
#SBATCH --time=00:30:00          # total run time limit (HH:MM:SS)
...

#SBATCH --mail-user=USERNAME@princeton.edu
...

# CHOOSE THE SCIRPT YOU NEED TO RUN BELOW
# =======================================

python import_mat_substract_background.py
# python import_mat.py

# =======================================
```

to que the job on the cluster use:
```
sbatch cue_import_on_cluster.py
88791 # Returns job id
```

Cancel the job
```
scancel 88791 #Job id
```