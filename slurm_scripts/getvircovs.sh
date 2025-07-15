#!/bin/bash
#
new="/scratch/project_2007362"

ml minimap2
ml samtools
c=4 # number of cores
r3=$1 # read fq.gz file
sname_assembly=$2 # assembly (contigs)

source global_env.sh
set -e

sname=$(basename $r3 | awk -F".fq.gz" '{print $1}')

vircontigs="${sname_assembly}_virmining/minedviruses.fna"
if [ $(grep -c "${sname}_" blacklist.txt) -gt 0 ];then
        echo "sample $sname in blacklist"
        exit
fi

cd $new/sandro/HeP_samples/2_assembly

mkdir -p ${sname_assembly}_virbams

sed "s/ /____/g" $vircontigs | sed "s:/:__:g" > tmp_${sname}.fna
samtools faidx tmp_${sname}.fna
awk 'BEGIN{FS="[>]"} /^>/{val=$2;next}  {print val"\t"length($0);val=""} END{if(val!=""){print val}}' tmp_${sname}.fna > tmpsizes_${sname}.tsv
minimap2 -t $c -ax sr tmp_${sname}.fna $r3 > ${sname}.sam
samtools view -h -F 2308 ${sname}.sam | samtools sort | samtools view -bh > ${sname}_map.bam
samtools view -h -f 4 ${sname}.sam | samtools sort | samtools view -bh > ${sname}_nomap.bam

bfile=${sname}_map.bam

samtools view $bfile | awk -F'\t' 'BEGIN{print "flag\toccurrences"} {a[$2]++} END{for(i in a)print i"\t"a[i]}' > ${sname}_map.flags
samtools view $bfile | cut -f1,3 | uniq | awk '{print $1"\t"$2}' > ${sname}_map.massign
awk -F'\t' '{n[$2]+=1}END{for(k in n)print k"\t"n[k]}' ${sname}_map.massign > ${sname}_map.rawcounts
samtools depth $bfile | awk -F'\t' '{if($3>0){n[$1]+=1;depth[$1]+=$3};if($3>=3){n3[$1]+=1;depth3[$1]+=$3};if($3>=5){n5[$1]+=1;depth5[$1]+=$3};if($3>=10){n10[$1]+=1;depth10[$1]+=$3}}END{for(contig in n){printf("%s\t%d\t%d\t%d\t%d\t",contig,n[contig],n3[contig],n5[contig],n10[contig]);if(n3[contig]==0){n3[contig]=1};if(n5[contig]==0){n5[contig]=1};if(n10[contig]==0){n10[contig]=1}print depth[contig]/n[contig]"\t"depth3[contig]/n3[contig]"\t"depth5[contig]/n5[contig]"\t"depth10[contig]/n10[contig] }}' > ${sname}_map.basescov
awk -F'\t' '{if(NR==FNR){n[$1]=$2}else{if($1 in n){print $1"\t"n[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}}' ${sname}_map.rawcounts ${sname}_map.basescov > ${sname}.depths
awk -F'\t' 'BEGIN{print "ID\tlength\traw_abundance\tbpcov1\tbpcov3\tbpcov5\tbpcov10\tmdepth1\tmdepth3\tmdepth5\tmdepth10"}{if(NR==FNR){n[$1]=$2}else{if($1 in n)print $1"\t"n[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}' tmpsizes_${sname}.tsv ${sname}.depths > ${sname}_map.cov

sed -i "s/____/ /g" ${sname}_map.cov
sed -i "s:__:/:g" ${sname}_map.cov
rm -f tmp_${sname}.fna tmpsizes_${sname}.tsv tmp_${sname}.fna.fai tmpsizes_${sname}.tsv
#for debugging
rm -f ${sname}.depths ${sname}_map.flags

samtools fastq -n $bfile > ${sname}_map.fq
samtools fastq -n ${sname}_nomap.bam > ${sname}_nomap.fq

gzip -c ${sname}_map.fq > ../5_phages_readmapback/map/${sname}_virmap.fq.gz
gzip -c ${sname}_nomap.fq > ../5_phages_readmapback/nomap/${sname}_virnomap.fq.gz
mv -f ${sname}_map.cov ../5_phages_readmapback/map/.
mv -f ${sname}_map.bam ${sname_assembly}_virbams/.

rm -f ${sname}.sam ${sname}_nomap*.bam ${sname}_map.fq ${sname}_nomap.fq ${sname}_map.basescov ${sname}_map.massign ${sname}_map.rawcounts

echo "Done job $sname"
