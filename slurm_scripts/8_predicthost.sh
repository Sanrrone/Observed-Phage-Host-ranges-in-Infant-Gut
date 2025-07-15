#!/bin/bash
#
#SBATCH --job-name=s8_iphop
#SBATCH --output=log/s8_iphop_%a.log
#SBATCH --error=err/s8_iphop_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=16:00:00
#SBATCH --array=47
#SBATCH --mem=150G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"
cpu=$(nproc)
n=$SLURM_ARRAY_TASK_ID

source global_env.sh

sname_assembly=`sed -n "${n} p" $sepsafile`
vircontigs="${home_project}/2_assembly/${sname_assembly}_virmining/minedviruses.fna"


export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate $IPHOP_ENV
####################################################################

cd ${home_project}/4_hostpred

sed "s/ /____/g" $vircontigs | sed "s:/:___:g" > ${sname_assembly}_tmp.fasta
set -ex

rm -rf ${sname_assembly}_splits ${sname_assembly}_sub
python $new/software/split_multifasta.py ${sname_assembly}_tmp.fasta ${sname_assembly}_sub 3

i=0

for ff in $(ls ${sname_assembly}_sub/*fasta)
do
	tmpname=$(basename $ff | sed "s/.fasta//g")
	iphop predict -d /scratch/project_2007362/software/iphop_db/Aug_2023_pub_rw -t $cpu -o ${sname_assembly}_${tmpname}_iphop -f ${ff} --no_qc &
        i=$((i+1))
        if [ $i -ge 3 ];then
                wait
                i=0
        fi
done
wait

echo "Virus,AAI to closest RaFAH reference,Host genus,Confidence score,List of methods" > ${sname_assembly}_iphop.csv
for ff in $(ls ${sname_assembly}_p*_iphop/Host_prediction_to_genus_m90.csv)
do
	tail -n +2 ${ff} | sed "s/____/ /g" | sed "s:___:/:g"

done >> ${sname_assembly}_iphop.csv

rm -rf ${sname_assembly}_tmp.fasta ${sname_assembly}_sub ${sname_assembly}_*_iphop


