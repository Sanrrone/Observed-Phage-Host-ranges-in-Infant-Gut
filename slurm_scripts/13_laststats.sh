#!/bin/bash
#
#SBATCH --job-name=summarize
#SBATCH --output=log/s14.log
#SBATCH --error=err/s14.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_207362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

source global_env.sh

cat $sepsafile | while read h
do 
	awk -v sep=$h 'BEGIN{FS="[>]"} /^>/{val=$2;next}  {print sep"\t"val"\t"length($0);val=""} END{if(val!=""){print val}}' ${home_project}/2_assembly/${h}_checkv/minedviruses.fna
done > supp_files/summary_contigs.tsv
