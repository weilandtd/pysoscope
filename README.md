
 <img src="./images/pysoscope_logo_v1-01.png" width="150"/>
 
# pysoscope 

## Getting started
You can install this package using `pip`
```
pip install -e /path/to/pysoscope
```

## Running the scripts on the cluster 
To run on the cluster first make the conda environment:

```
conda create --name env_maldi --file requirements.txt
conda acivate env_maldi
pip install - e /path/to/pysoscope
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
