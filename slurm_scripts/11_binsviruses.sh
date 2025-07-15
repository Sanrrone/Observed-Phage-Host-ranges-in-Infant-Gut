#!/bin/bash
#
#SBATCH --job-name=s11_vrhyme
#SBATCH --output=log/s11_vrhyme_%a.log
#SBATCH --error=err/s11_vrhyme_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=00:15:00 # set to 3 hours when all samples are run
#SBATCH --array=1-55 #1-55
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

n=$SLURM_ARRAY_TASK_ID
c=$(nproc)

source global_env.sh

set -e

sname_assembly=`sed -n "${n} p" $sepsafile`
vircontigs="${sname_assembly}_virmining/minedviruses.fna"

export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate $VRHYME_ENV

cd ${home_project}/2_assembly
bamf=${sname_assembly}_virbams

sed "s/ /____/g" $vircontigs > tmp_${sname_assembly}.fna
rm -rf ${sname_assembly}_vrhyme $bamf/*.bai
vRhyme -i tmp_${sname_assembly}.fna -b $bamf/*_map.bam -t $c -o ${sname_assembly}_vrhyme
tsv=$(ls ${sname_assembly}_vrhyme/vRhyme_best_bins.*.membership.tsv | tail -n 1 )
sed "s/____/ /g" $tsv > ${sname_assembly}_vrhyme/${sname_assembly}_vrhyme.tsv
rm -rf ${sname_assembly}_vrhyme/vRhyme_alternate_bins ${sname_assembly}_vrhyme/vRhyme_best_bins_fasta $bamf/*.bai

