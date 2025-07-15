#!/bin/bash
#
new="/scratch/project_2007362"

source global_env.sh

r3=$1
sname=$(basename $r3 | sed "s/.fq.gz//g")
c=4

#ml kraken
#ml bracken
ml minimap2
##############################################################################################
cd ${home_project}/6_tax_profile/bact
set -ex

#kraken2 --db $humgutdb --confidence 0.1 --threads $cpu --use-names --report ${sname}.kreport --output /dev/null ${home_project}/5_phages_readmapback/nomap/${sname}_virnomap.fq.gz

#bracken -d $new/software/krakenDB/uhggv2 -i ${sname}.kreport -o ${sname}_s.bracken -r 250 -l S -t 5
#bracken -d $new/software/krakenDB/uhggv2 -i ${sname}.kreport -o ${sname}_g.bracken -r 250 -l G -t 5
#bracken -d $new/software/krakenDB/uhggv2 -i ${sname}.kreport -o ${sname}_f.bracken -r 250 -l F -t 5

#rm -f ${sname}_bracken_*

#$new/software/centrifuge-1.0.4/centrifuge -x $new/software/centrifugeDB/nt/nt -p $cpu -q -U $new/sandro/HeP_samples/5_phages_readmapback/nomap/${sname}_virnomap.fq.gz -S ${sname}_centrifuge_report.tsv --report-file /dev/null
#awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8}' ${sname}_centrifuge_report.tsv > ${sname}.tmp
#rm -f ${sname}_centrifuge_report.tsv && mv ${sname}.tmp ${sname}_centrifuge_report.tsv

#$new/software/centrifuge-1.0.4/centrifuge-kreport -x $new/software/centrifugeDB/nt/nt ${sname}_centrifuge_report.tsv > ${sname}_centrifuge.kreport
#rm -f ${sname}_centrifuge_report.tsv

minimap2 -t $c --split-prefix=${sname}_tmp -a $new/software/sep_bactDB/m2/sepB.mmi ${home_project}/5_phages_readmapback/nomap/${sname}_virnomap.fq.gz | samtools view -h -F 2308 | samtools sort | samtools view -bh > ${sname}_map.bam
#minimap2 -t $c --split-prefix=${sname}_tmp -a $new/software/HumGutDB/m2/hg.mmi ${home_project}/5_phages_readmapback/nomap/${sname}_virnomap.fq.gz | samtools view -h -F 4 | samtools view -bh > ${sname}.tmp.bam

#for flag in 99 147 83 163 67 131 115 179
#do
#        samtools -bh -f $flag ${sname}.tmp.bam > ${sname}.flag_$flag.bam
#done	
#samtools merge -o ${sname}_map.bam ${sname}.flag_*
#rm ${sname}.flag_* ${sname}.tmp.bam

bfile=${sname}_map.bam

samtools sort $bfile > ${sname}_tmp.bam
rm $bfile && mv ${sname}_tmp.bam $bfile
samtools view $bfile | awk -F'\t' 'BEGIN{print "flag\toccurrences"} {a[$2]++} END{for(i in a)print i"\t"a[i]}' > ${sname}_map.flags
samtools view $bfile | cut -f1,3 | uniq | awk '{print $1"\t"$2}' > ${sname}_map.massign
awk -F'\t' '{n[$2]+=1}END{for(k in n)print k"\t"n[k]}' ${sname}_map.massign > ${sname}_map.rawcounts
samtools depth $bfile | awk -F'\t' '{if($3>0){n[$1]+=1;depth[$1]+=$3};if($3>=3){n3[$1]+=1;depth3[$1]+=$3};if($3>=5){n5[$1]+=1;depth5[$1]+=$3};if($3>=10){n10[$1]+=1;depth10[$1]+=$3}}END{for(contig in n){printf("%s\t%d\t%d\t%d\t%d\t",contig,n[contig],n3[contig],n5[contig],n10[contig]);if(n3[contig]==0){n3[contig]=1};if(n5[contig]==0){n5[contig]=1};if(n10[contig]==0){n10[contig]=1}print depth[contig]/n[contig]"\t"depth3[contig]/n3[contig]"\t"depth5[contig]/n5[contig]"\t"depth10[contig]/n10[contig] }}' > ${sname}_map.basescov
awk -F'\t' '{if(NR==FNR){n[$1]=$2}else{if($1 in n){print $1"\t"n[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}}' ${sname}_map.rawcounts ${sname}_map.basescov > ${sname}.depths
#kraken:taxid|3018277|HumGut_18277_1
#for humgutdb
#awk -F'\t' 'BEGIN{print "ID\tlength\traw_abundance\tbpcov1\tbpcov3\tbpcov5\tbpcov10\tmdepth1\tmdepth3\tmdepth5\tmdepth10"}{if(NR==FNR){n[$1]=$16}else{split($1,a,"|");$1=a[3];split(a[3],b,"_");c=b[1]"_"b[2];if(c in n)print $1"\t"n[c]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}' $HUMGUTTSV ${sname}.depths > ${sname}.bactstats

#for sepbactdb
awk -F'\t' 'BEGIN{print "ID\tlength\traw_abundance\tbpcov1\tbpcov3\tbpcov5\tbpcov10\tmdepth1\tmdepth3\tmdepth5\tmdepth10"}{if(NR==FNR){n[$1]=$2}else{if($1 in n){print $1"\t"n[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}}' $new/software/sep_bactDB/sep_bact.tsv ${sname}.depths > ${sname}.bactstats
rm -f ${sname}_map.massign ${sname}_map.rawcounts ${sname}_map.basescov ${sname}.depths ${sname}_map.flags ${sname}_map.bam
