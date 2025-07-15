#!/bin/bash
#
#SBATCH --job-name=s2_snowball
#SBATCH --output=log/s2_%a.log
#SBATCH --error=err/s2_%a.err
#SBATCH --partition=hugemem
#SBATCH --array=32,51
#SBATCH --account=Project_2007362
#SBATCH --time=2-10:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=400G
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

source global_env.sh

n=$SLURM_ARRAY_TASK_ID
c=$(nproc)
 
sname_mag=`sed -n "${n} p" $sepsafile`
#samples=$(awk -F"\t" -v s="$sname_mag" '{if($1"_"$4==s){print $2}}' $sampletable)
samples=$(awk -F"\t" -v s=$sname_mag -v p="/scratch/project_2007362/sandro/HeP_samples/1_hostremoval/" '{if($1"_"$4==s){printf p$2".fq.gz "}}' $sampletable)

cd ${home_project}/2_assembly/
rm -rf ${sname_mag}_tmpassembly && mkdir ${sname_mag}_tmpassembly
cd ${sname_mag}_tmpassembly

cat $samples > all.fq.gz
rm -rf final
megahit --presets meta-sensitive -r all.fq.gz -t $c -m 400e9 -o final
perl $SOFT_HOME/removesmalls.pl $lfilt final/final.contigs.fa | sed "s/>\([^ ]*\) .*/>\1/g" > ../${sname_mag}.fna
rm all.fq.gz
cd ..
rm -rf ${sname_mag}_tmpassembly
