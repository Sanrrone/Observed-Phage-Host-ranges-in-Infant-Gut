#!/bin/bash
#
new="/scratch/project_2007362"

source global_env.sh

r3=$1
sname=$(basename $r3 | sed "s/.fq.gz//g")
c=2

ml kraken
##############################################################################################
cd ${home_project}/6_tax_profile/bact
set -ex

kraken2 --db $new/software/sep_bactDB/k2/sepbact --confidence 0.5 --threads $c --use-names --report ${sname}.k2c05.kreport --output /dev/null ${home_project}/5_phages_readmapback/nomap/${sname}_virnomap.fq.gz

#bracken -d $new/software/krakenDB/uhggv2 -i ${sname}.kreport -o ${sname}_s.bracken -r 250 -l S -t 5
#bracken -d $new/software/krakenDB/uhggv2 -i ${sname}.kreport -o ${sname}_g.bracken -r 250 -l G -t 5
#bracken -d $new/software/krakenDB/uhggv2 -i ${sname}.kreport -o ${sname}_f.bracken -r 250 -l F -t 5

