#!/bin/bash
#
#SBATCH --job-name=s7_1_phaxxx
#SBATCH --output=log/s7_pha_%a.log
#SBATCH --error=err/s7_pha_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=00:40:00
#SBATCH --array=7-9,38-41,30-32,36,25-29,12,17,18
#SBATCH --mem=5G
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
PDIR="$new/software/PhaBOX"

export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate $PHABOX_ENV
##############################################################################################

cd $new/sandro/HeP_samples/3_phage_annotation
virfasta=${home_project}/2_assembly/$vircontigs

rm -rf ${sname_assembly}_tmpgcn ${sname_assembly}_phagcn ${sname_assembly}_phatyp ${sname_assembly}_tmptyp

set -ex
phabox2 --task phagcn --dbdir $PDIR/phabox_db_v2 --contigs $virfasta --threads $c --midfolder ${sname_assembly}_tmpgcn --outpth ${sname_assembly}_phagcn --reject 0.5 --len $lfilt
mv ${sname_assembly}_phagcn/final_prediction/phagcn_prediction.tsv ${sname_assembly}_phagcn.tsv

phabox2 --task phatyp --dbdir $PDIR/phabox_db_v2 --contigs $virfasta --threads $c --midfolder ${sname_assembly}_tmptyp --outpth ${sname_assembly}_phatyp --reject 0.5 --len $lfilt
mv ${sname_assembly}_phatyp/final_prediction/phatyp_prediction.tsv ${sname_assembly}_phatyp.tsv

rm -rf ${sname_assembly}_phagcn ${sname_assembly}_phatyp

