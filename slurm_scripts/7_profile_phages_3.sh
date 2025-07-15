#!/bin/bash
#
#SBATCH --job-name=s7_3_virustaxo
#SBATCH --output=log/s7_vtaxo_%a.log
#SBATCH --error=err/s7_vtaxo_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=00:15:00
#SBATCH --array=1-55
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"
c=$(nproc)

source global_env.sh

n=$SLURM_ARRAY_TASK_ID
sname_assembly=`sed -n "${n} p" $sepsafile`
vircontigs="${sname_assembly}_virmining/minedviruses.fna"
PDIR="$new/software/PhaBOX-main"

export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate $VIRUSTAXO_ENV
##############################################################################################
cd ${home_project}/3_phage_annotation
SDIR="/scratch/project_2007362/software/VirusTaxo"

set -e
python $SDIR/predict.py --model_path $SDIR/uhgv_genus_mq.pkl --seq ${home_project}/2_assembly/$vircontigs > ${sname_assembly}_virustaxo.tsv


