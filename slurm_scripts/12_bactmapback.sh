#!/bin/bash
#
#SBATCH --job-name=s12_bactmap
#SBATCH --output=log/s12_%a.log
#SBATCH --error=err/s12_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --array=1-55
#SBATCH --time=2:00:00 # set to 3 hours when all samples are run # 8h
#SBATCH --cpus-per-task=20
#SBATCH --mem=15G #150
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
	xargs -P 6 -n 1 -I {} bash ./getbactcovs.sh {} ${sname_mag}

wait
#P 10 parallel processes
#n 1 command per independent seq 1,2,3,4...
#I {} the value of seq command, 1,2,3,4...

