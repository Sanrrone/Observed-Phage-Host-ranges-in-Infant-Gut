#!/bin/bash
#
#SBATCH --job-name=hepB_stats
#SBATCH --output=log/run_%a.log
#SBATCH --error=err/run_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-55
#SBATCH --mem=30G
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

n=$SLURM_ARRAY_TASK_ID
c=$(nproc)

new="/scratch/project_2007362"
sname=`sed -n "${n} p" hep_sa.txt`
h=$sname
bdb="$new/software/hep_bactDB/hep_passbact.fna"
ml checkm2/1.1.0
ml seqtk

#rm -rf $sname
mkdir -p $sname
cd $sname

mkdir -p genomes
grep ">${sname}__" $bdb | sed "s/>//g"> cids.txt
awk -F"__" '{print $2}' cids.txt | sort -u | while read bid
do
	grep "__${bid}__" cids.txt > bids.txt
	seqtk subseq $bdb bids.txt > genomes/${bid}.fna
done
rm -rf ${h}_checkm ${h}_checkm2.tsv
checkm2 predict --threads $c --input genomes --output-directory ${h}_checkm --tmpdir $new/tmp --remove_intermediates
cp ${h}_checkm/quality_report.tsv ${h}_checkm2.tsv
rm -rf ${h}_checkm bids.txt  cids.txt
