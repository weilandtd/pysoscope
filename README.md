# pysocope

# Convert maldi files to npz

To run on the cluster first make the conda environment:

```
conda create --name env_maldi --file reqrequirements.txt
```

Modify the SLURM script with your username

```
#SBATCH --mail-user=USERNAME@princeton.edu
```

Que the job
```
sbatch convery_maldi_script.py
88791 # Returns job id
```

Cancel the job
```
scancel 88791 #Job id
```