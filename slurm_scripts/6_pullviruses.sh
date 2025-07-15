#!/bin/bash
#
#SBATCH --job-name=s6_dvppha
#SBATCH --output=log/s6_dvppha_%a.log
#SBATCH --error=err/s6_dvppha_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=03:30:00 #40min
#SBATCH --array=1-55 # 1-55
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

source global_env.sh

n=$SLURM_ARRAY_TASK_ID
sname=`sed -n "${n} p" $sepsafile`

inputF=after_mpp.fna #input
virfasta=minedviruses.fna #output
SDIR=$new/software
cpu=$(nproc)

cd ${home_project}/2_assembly/${sname}_virmining

###############################################################################
export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate $DVP_ENV

python $SDIR/DeepVirFinder/dvf.py -i $inputF -m $SDIR/DeepVirFinder/models  -o . -l $lfilt -c $cpu
mv ${inputF}_gt${lfilt}bp_dvfpred.txt ${sname}_dvp.tsv

conda deactivate

awk -F'\t' -v fscore=$fscore '{if($3>=fscore && $4<0.01)print $1}' ${sname}_dvp.tsv | sort -u > ${sname}_dvp_exclusive.txt
awk -v elist="${sname}_dvp_exclusive.txt" 'BEGIN{while((getline<elist)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $inputF > after_dvp.fna
seqkit grep -w 0 -nf ${sname}_dvp_exclusive.txt $inputF > dvpviruses.fna

###########################################
conda activate $VS2_ENV

rm -rf sub
nf=3

rm -f chunk* vs2_*
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < after_dvp.fna > oneline.fna
awk 'BEGIN {n=0;} /^>/ {if(n%500==0){file=sprintf("chunk%d.fa",n);} print >> file; n++; next;} { print >> file; }' < oneline.fna
i=0
for ff in $(ls chunk*)
do
        virsorter run -w vs2_${ff} -i $ff  --min-length 2000 -j 4 --min-score 0.9 --high-confidence-only --provirus-off --rm-tmpdir all &
        i=$((i+1))
        if [ $i -ge 10 ];then
                wait
                i=0
 
        fi
done
wait

cp vs2_*/final-viral-combined.fa > vs2viruses.fna
rm -rf chunk* vs2_*

conda deactivate
############################################


cat step3.fna dvpviruses.fna vs2viruses.fna > $virfasta
sed -i "s/ /_checkv_/g" $virfasta

#rm -f after_*

awk 'BEGIN{FS="[>]"} /^>/{val=$2;next}  {print val"\t"length($0);val=""} END{if(val!=""){print val}}' $virfasta > tmplen.tsv

ml blast
makeblastdb -in $virfasta -out blastdb -dbtype nucl
blastn -query $virfasta -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -num_threads $cpu -max_target_seqs 100000 -perc_identity 50
python $new/software/aaicluster/blastani.py -i blast.tsv -o ani.tsv
python $new/software/aaicluster/cluster.py --fna $virfasta --ani ani.tsv --out ${sname}_ani_s.tsv --min_ani 98 --min_qcov 0 --min_tcov 85 
rm -f blast.tsv ani.tsv blastdb*
rm -f representatives.txt
awk -F"\t" '{if(NR==FNR){n[$1]=$2}else{if($1==$2){print $1 >> "representatives.txt"};if($1!=$2){split($2,a,",");mlen=0;id="";for(c in a){if(mlen<n[c]){mlen=n[c];id=c};print c>>"representatives.txt"}}}}' tmplen.tsv ${sname}_ani_s.tsv
rm -f tmplen.tsv

seqkit grep -w 0 -nf representatives.txt $virfasta > c98.fna
perl $new/software/removesmalls.pl $lfilt c98.fna > minedviruses.fna

awk -v sep="$sname" 'BEGIN{FS="[>]"} /^>/{val=$2;next}  {print sep"\t"val"\t"length($0);val=""} END{if(val!=""){print val}}' minedviruses.fna >> ~/phage_approach/supp_files/summary_contigs.tsv
grep "phage" minedviruses.fna | awk -F"_" -v h=$sname '{gsub(">","",$0);print h"\t"$0"\t"$1"_"$2"\t"$4"\t"$5}' >> ~/phage_approach/supp_files/prophages_coordinates.tsv
grep "checkv" minedviruses.fna | awk -F"_" '{gsub(">","",$0);print $0"\t"$1"_"$2"\t"$5}' | awk -F"\t" -v h=$sname '{split($3,a,"-");split(a[2],b,"/");print h"\t"$1"\t"$2"\t"a[1]"\t"b[1]}' >> ~/phage_approach/supp_files/prophages_coordinates.tsv

rm c98.fna representatives.txt
echo "DONE job $n step 6"


