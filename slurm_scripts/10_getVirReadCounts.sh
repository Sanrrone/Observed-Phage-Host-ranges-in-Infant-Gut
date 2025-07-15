#!/bin/bash
#
#SBATCH --job-name=s10_rema2
#SBATCH --output=log/s10_%a.log
#SBATCH --error=err/s10_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --array=46-55 #dont perform all at the same time, disk quota. de a 15 ta bien
#SBATCH --time=5:00:00 # set to 3 hours when all samples are run
#SBATCH --cpus-per-task=10
#SBATCH --mem=11G
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

new="/scratch/project_2007362"
n=$SLURM_ARRAY_TASK_ID
source global_env.sh
 
sname_mag=`sed -n "${n} p" $sepsafile`
sep=$(echo $sname_mag | awk -F"_" '{print $1}')
sa=$(echo $sname_mag | sed "s/${sep}_//g" )

set -ex 
awk -F"\t" -v h=$sep -v sa=$sa -v p="/scratch/project_2007362/sandro/HeP_samples/1_hostremoval/" '{if($1==h && $4==sa){print p$2".fq.gz"}}' $sampletable | 
	xargs -P 10 -n 1 -I {} bash ./getvircovs.sh {} ${sname_mag}

wait
#P 10 parallel processes
#n 1 command per independent seq 1,2,3,4...
#I {} the value of seq command, 1,2,3,4...

