#!/bin/bash
#
#SBATCH --job-name=s4_filter
#SBATCH --output=log/s4_%a.log
#SBATCH --error=err/s4_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=00:15:00
#SBATCH --array=1-55 # 1-41
#SBATCH --cpus-per-task=1
#SBATCH --mem=22G
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

source global_env.sh

n=$SLURM_ARRAY_TASK_ID
sname=`sed -n "${n} p" $sepsafile`

cd ${home_project}/2_assembly

########### check_v
cd ${sname}_virmining
inputF=${sname}_checkv/${sname}_after_checkv_phageboost.fna

perl $new/software/removesmalls.pl $lfilt $inputF > tmp.fna
rm $inputF && mv tmp.fna $inputF


### finding false positives
export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate k2

cs=0.5
kraken2 --db $humgutdb --confidence $cs --report /dev/null --output ${sname}.kreads $inputF
awk -F'\t' '{if($1=="C"){print}}' ${sname}.kreads |
        awk -F"\t" -v cs=$cs 'function sum(array){sumn=0;for(k in array){sumn+=array[k]}return(sumn)};{split($NF,kpairs," ");delete n;for(i in kpairs){split(kpairs[i],kp,":");n[kp[1]]+=kp[2]};asort(n, sN);score=sN[length(sN)]/sum(sN);if(score>=cs)print $2"\t"$3"\t"score}' > potential_bact_contigs.txt

##tax
tidarr=$(awk -F'\t' '{print $2}' potential_bact_contigs.txt | sort -u | paste -s -d, -)

if [ "$tidarr" != "" ];then

echo "
library(dplyr)
library(tibble)

kr2df<-read.table('$humgutdb/inspect.txt', sep = '\t', quote = '')
kr2df\$V6<-trimws(kr2df\$V6)
colnames(kr2df)<-c('relabu','rawabu','rawabulvl','taxlvl','tid','name')
kr2df<-kr2df[,c('tid','taxlvl','name')]
tlvl<-c('D'=1,'P'=2,'C'=3,'O'=4,'F'=5,'G'=6,'S'=7)

get_idx<-function(x, lvl){
  tmpidx<-which(kr2df[1:x,'taxlvl']==lvl)
  ifelse(length(tmpidx)>0,max(tmpidx),1000000)
}

t_lineage<-function(tid){
  idx<-which(kr2df\$tid==tid)
  ldf<-data.frame(tid=tid,D=NA,P=NA,C=NA,O=NA,F=NA,G=NA,S=NA)
  oriname<-kr2df[idx,'name']

  lowest_idx<-get_idx(idx,'D')
  for(lvl in tlvl){
    curr_idx<-get_idx(idx,names(tlvl[lvl]))
    if(lowest_idx<=curr_idx){
      ldf[1,names(tlvl[lvl])]<-kr2df[curr_idx,'name']
      lowest_idx<-curr_idx
    }
  }
  
  colnames(ldf)<-c('tid','domain','phylum','class','order','family','genus','species')
  ldf\$species<-ifelse(is.na(ldf\$family),oriname,ldf\$species)
  return(ldf)
}


taxdf<-lapply(c($tidarr), function(x)t_lineage(x)) %>% bind_rows()
write.table(taxdf[,c('tid','domain')],'kingdoms.tsv', quote = F, row.names = F, col.names = F, sep = '\t')

" > taxr.R

module load r-env/430
srun apptainer_wrapper exec Rscript --vanilla taxr.R
awk -F'\t' '{if(NR==FNR){n[$1]=$2}else{if(n[$2]=="Eukaryota" || n[$2]=="Bacteria"){print $1}}}' kingdoms.tsv potential_bact_contigs.txt > bact_contigs.txt
else
	echo "" > bact_contigs.txt
fi

#####
awk -v elist="bact_contigs.txt" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $inputF > after_bactfilt.fna
cp potential_bact_contigs.txt ${sname}_kraken_bacttag.txt
rm -f bact_contigs.txt kingdoms.tsv resolved.json potential_bact_contigs.txt ${sname}.kreads potential_bact_contigs.txt taxr.R

echo "DONE job $n step 4"
