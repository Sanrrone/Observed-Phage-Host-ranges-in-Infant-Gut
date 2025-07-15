### general path files
home_project=$new/sandro/S_samples # main path as starting to point to build further paths
sampletable="supp_files/supp_table1.tsv" #tsv file with the samples
onlysep="supp_files/sep.txt"
sepsafile="supp_files/sep_sa.txt"
sfile="supp_files/samples.txt"

#### software parameters
fscore=0.9 # confidence score for software
lfilt=2000 # minimum contig length for viruses (or in general)


## db
humgutdb=$new/software/unphage_humgutdb/k2/hg
humgutdb_m2=$new/software/unphage_humgutdb/m2/hg.mmi
UHGVDB=$new/software/UHGV_mqplus/blast/uhgv
HUMGUTTSV=$new/software/HumGutDB/HumGut.tsv

### enviroments & software paths
SOFT_HOME=$new/software
checkv_phageboost_env="prophages" # conda activate $checkv_phageboost_env
MPP_ENV="metaphapred" # conda activate $MPP_ENV
mpp_home=/scratch/project_2007362/software/MetaPhaPred
PHABOX_ENV="phabox2"
VIRUSTAXO_ENV=virustaxo
IPHOP_ENV=iphop
VRHYME_ENV=vRhyme
DVP_ENV=dvp
VS2_ENV=vs2
PROP_ENV=propagate
