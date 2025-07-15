#!/bin/bash
#
#SBATCH --job-name=s7_2_blast
#SBATCH --output=log/s7_blast_%a.log
#SBATCH --error=err/s7_blast_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=00:50:00
#SBATCH --array=1-55
#SBATCH --mem=6G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

source global_env.sh

n=$SLURM_ARRAY_TASK_ID

sname_assembly=`sed -n "${n} p" $sepsafile`
vircontigs="${home_project}/2_assembly/${sname_assembly}_virmining/minedviruses.fna"
ml blast

##############################################################################################
cd ${home_project}/3_phage_annotation

set -e

######### UHGV_mqplus
blastn -query $vircontigs -db $UHGVDB -num_threads 2 -evalue 0.001 -perc_identity 70 -out ${sname_assembly}.tsv -qcov_hsp_perc 70 -outfmt "6 qseqid sseqid pident qlen slen length"
awk -F"\t" '{if($3>=70 && $6/$4 >= 0.7 && $6/$5 >= 0.7)print $1"\t"$2"\t"$3"\t"($6/$4)*100"\t"($6/$5)*100}' ${sname_assembly}.tsv |
        awk -F"\t" 'BEGIN{getline;n[$1]=$3*$4;v[$1]=$0}{if($1 in n){if(n[$1]<$3*$4){n[$1]=$3*$4;v[$1]=$0}}else{n[$1]=$3*$4;v[$1]=$0}}END{for(k in n)print v[k]}'> ${sname_assembly}_blast.tsv

rm -f ${sname}.tsv

#####################
#TMP=emg_${sname_assembly}_tmp
# rm -rf $TMP && mkdir $TMP
#nextflow run EBI-Metagenomics/emg-viral-pipeline -r master --fasta "$new/sandro/HeP_samples/2_assembly/$vircontigs" --cores 1 --max-cores 1 -profile "local,singularity" --memory 8 --output testoutput --workdir $TMP

#kraken2 --db $UHGVDB --confidence 0.1 --report /dev/null --output ${sname_assembly}.kreads $vircontigs

#awk -F'\t' '{if($1=="C"){print}}' ${sname_assembly}.kreads |
#        awk -F"\t" 'function sum(array){sumn=0;for(k in array){sumn+=array[k]}return(sumn)};{split($NF,kpairs," ");delete n;for(i in kpairs){split(kpairs[i],kp,":");n[kp[1]]+=kp[2]};asort(n, sN);print $2"\t"$3"\t"sN[length(sN)]/sum(sN)}' > ${sname_assembly}.ktsv
#
#rm -f ${sname_assembly}.kreads

