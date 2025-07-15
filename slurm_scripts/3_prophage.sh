#!/bin/bash
#
#SBATCH --job-name=s3_prophage
#SBATCH --output=log/s3_%a.log
#SBATCH --error=err/s3_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=08:00:00
#SBATCH --array=1-55
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

source global_env.sh

n=$SLURM_ARRAY_TASK_ID
c=$SLURM_CPUS_PER_TASK
###############################################################################
################## ACTIVAR CONDA ENVIROMENT "prophages" #######################
###############################################################################
export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate $checkv_phageboost_env
ml blast
ml seqkit
##############################################################################

sname=`sed -n "${n} p" $sepsafile`
fscore=0.9

cd ${home_project}/2_assembly
echo "working on $sname"

rm -rf ${sname}_virmining
mkdir -p ${sname}_virmining
cd ${sname}_virmining

######### UHGV_mqplus
blastn -query ../${sname}.fna -db $UHGVDB -num_threads $c -evalue 0.001 -perc_identity 70 -out ${sname}.tsv -qcov_hsp_perc 70 -outfmt "6 qseqid sseqid pident qlen slen length"
awk -F"\t" '{if($3>=70 && $6/$4 >= 0.7)print $1"\t"$2"\t"$3"\t"($6/$4)*100}' ${sname}.tsv | 
	awk -F"\t" 'BEGIN{getline;n[$1]=$3*$4;v[$1]=$0}{if($1 in n){if(n[$1]<$3*$4){n[$1]=$3*$4;v[$1]=$0}}else{n[$1]=$3*$4;v[$1]=$0}}END{for(k in n)print v[k]}'> ${sname}_blast.tsv
rm -f ${sname}.tsv
awk -F"\t" '{print $1}' ${sname}_blast.tsv > ${sname}_blast_exclusive.txt

awk -v elist="${sname}_blast_exclusive.txt" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' ../${sname}.fna > ${sname}_after_blast.fna
seqkit grep -w 0 -nf ${sname}_blast_exclusive.txt ../${sname}.fna > blast.fna  


checkv end_to_end ${sname}_after_blast.fna ${sname}_checkv -t 4 --remove_tmp
cd ${sname}_checkv
#contig_id       contig_length   total_genes     viral_genes     host_genes
touch proviruses.fna
grep ">" proviruses.fna | sed "s/>//g" | sort -u > proids.txt

awk -F'\t' -v lfilt=$lfilt '{
if(NR==FNR){
        n[$1]=0
}else{
if($1 in n == 0){
	if($6=="Yes" && $2>lfilt){
       		 if($4==0 && $5==1){print;next}
       		 if($4>$5){
       		         print $1
       		 }else{
       		 	split($9,a,",");
      	 	 	n[a[1]]=$7;
       	 		n[a[2]]=$8;
       		 	if(n["viral"]>n["host"]){print $1}
       	 	}
	}
}
}}' proids.txt contamination.tsv > nocont.txt

cat nocont.txt proids.txt | sed "s/Scaffold_\([0-9]\+\).*/Scaffold_\1/g" > ${sname}_checkv_valids.txt
#rm -f nocont.txt tmp.fna proviruses.fna proids.txt

awk -v elist="${sname}_checkv_valids.txt" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' ../${sname}_after_blast.fna > ${sname}_after_checkv.fna
seqkit grep -w 0 -nf nocont.txt viruses.fna > tmp.fna
cat proviruses.fna tmp.fna ../blast.fna > step1.fna
rm -f tmp.fna
#####################
#phigaro -f tmp.fna -t 16 -e bed -o ${sname}_phigaro -d & # -d automatically filter 20k contigs
PhageBoost -f ${sname}_after_checkv.fna -j $c -o ${sname}_phageboost -meta 0 -t $fscore
#wait $!

#parse phageboost
find ${sname}_phageboost -name "*.fasta" -exec cat {} \; > ${sname}_prophages.fasta
perl $new/software/removesmalls.pl $lfilt ${sname}_prophages.fasta > tmp.fna
rm ${sname}_prophages.fasta && mv tmp.fna ${sname}_prophages.fasta

grep ">" ${sname}_prophages.fasta | sed "s/>//g" > ${sname}_phageboost_valids.txt


#parse phigaro
#bedtools getfasta -fo ${sname}_prophages_tmp.fasta -fi ${sname}.fna -bed ${sname}_phigaro.bed
#grep ">" ${sname}_prophages_tmp.fasta | sed "s/>//g" >  ${sname}_phigaro_exclusive.txt
#cat ${sname}_prophages_tmp.fasta >> ${sname}_prophages.fasta
#rm -f ${sname}_prophages_tmp.fasta


#cut -f1 ${sname}_phigaro.bed | sort -u > ${sname}_exclude.txt
awk -F"_phage" '{print $1}' ${sname}_phageboost_valids.txt | sort -u > ${sname}_phageboost_exclusive.txt
awk -v elist="${sname}_phageboost_exclusive.txt" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' ${sname}_after_checkv.fna > ${sname}_after_checkv_phageboost.fna 
rm -rf ${sname}_phageboost ${sname}.fna.fai

#########################
cat step1.fna ${sname}_prophages.fasta  > step2.fna
rm -f step1.fna ${sname}_prophages.fasta ${sname}_after_checkv.fna


##### getting host prediction
cat ${sname}_phageboost_exclusive.txt ${sname}_checkv_valids.txt > ${sname}_exclude.txt
seqkit grep -w 0 -nf ${sname}_exclude.txt ../${sname}_after_blast.fna > ${sname}_excluded.fna
export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate k2

kraken2 --confidence 0.3 --db $humgutdb --report /dev/null --output ${sname}.kreads ${sname}_excluded.fna
awk -F'\t' '{if($1=="C"){print}}' ${sname}.kreads |
        awk -F'\t' 'function sum(array){sumn=0;for(k in array){sumn+=array[k]}return(sumn)};{split($NF,kpairs," ");delete n;for(i in kpairs){split(kpairs[i],kp,":");n[kp[1]]+=kp[2]};asort(n, sN);print $2"\t"$3"\t"sN[length(sN)]/sum(sN)}' > ${home_project}/4_hostpred/${sname}_kr.tsv


rm -f ${sname_assembly}.kreads ${sname}_excluded.fna ${sname}_exclude.txt viruses.fna proviruses.fna nocont.txt proids.txt ${sname}_phageboost_exclusive.txt



echo "Done job $n"
