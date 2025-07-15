#!/bin/bash
#
#SBATCH --job-name=s5_pullphages
#SBATCH --output=log/s5_pullphages_%a.log
#SBATCH --error=err/s5_pullphages_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=06:00:00
#SBATCH --array=18,25,26,27,28,34,37-40 # 1-55
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

source global_env.sh

n=$SLURM_ARRAY_TASK_ID
sname=`sed -n "${n} p" $sepsafile`
inputF=after_bactfilt.fna
SDIR=$new/software
PDIR="$new/software/PhaBOX"
c=$(nproc)

set -e

cd ${home_project}/2_assembly/${sname}_virmining

###############################
export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate $PHABOX_ENV

rm -rf ${sname}_tmpphamer
phabox2 --task phamer --dbdir $PDIR/phabox_db_v2 --contigs $inputF --threads $c --midfolder ${sname}_tmpphamer --outpth ${sname}_phamer --len $lfilt
 
set -e
 
rm -rf ${sname}_phamer/${sname}_tmpphamer
 
conda deactivate
awk -F'\t' -v fscore=$fscore '{if($5>=fscore)print $1}' ${sname}_phamer/final_prediction/phamer_prediction.tsv > ${sname}_phamer_exclusive.txt
awk -v elist="${sname}_phamer_exclusive.txt" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $inputF > after_phamer.fna
seqkit grep -w 0 -nf ${sname}_phamer_exclusive.txt $inputF > phages_phamer.fna

##############################

conda activate $MPP_ENV

rm -f chunk* *_mpp.tsv
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < after_phamer.fna > oneline.fna
awk 'BEGIN {n=0;} /^>/ {if(n%30==0){file=sprintf("chunk%d.fa",n);} print >> file; n++; next;} { print >> file; }' < oneline.fna
i=0
for ff in $(ls chunk*)
do
	python $mpp_home/predict.py -i $ff -o ${ff}_mpp.tsv &
	i=$((i+1))
	if [ $i -ge 5 ];then
		wait
		i=0

	fi
done
wait
echo -e "name\tclen\tscore" > ${sname}_mpp.tsv
cat chunk*_mpp.tsv >> ${sname}_mpp.tsv
rm -f oneline.fna chunk*

conda deactivate
awk -F'\t' -v fscore=$fscore '{if($3>=fscore)print $1}' ${sname}_mpp.tsv > ${sname}_mpp_exclusive.txt && rm -f ${sname}_mpp.tsv
awk -v elist="${sname}_mpp_exclusive.txt" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' after_phamer.fna > after_mpp.fna
seqkit grep -w 0 -nf ${sname}_mpp_exclusive.txt after_phamer.fna > phages_mpp.fna

#######################################

cat ${sname}_checkv/step2.fna phages_phamer.fna phages_mpp.fna > step3.fna


echo "DONE job $n step 5"



