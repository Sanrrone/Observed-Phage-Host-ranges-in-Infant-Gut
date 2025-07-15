#!/bin/bash
#
#SBATCH --job-name=s9_cluster
#SBATCH --output=log/s9_met2_%a.log
#SBATCH --error=err/s9_met2_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=04:00:00
#SBATCH --array=1-55
#SBATCH --mem=60G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

n=$SLURM_ARRAY_TASK_ID

source global_env.sh

pair=`sed -n "${n} p" $sepsafile`
sname=$pair
cpu=$(nproc)
method=2

export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"

##############################################################################################
cd ${home_project}/2_assembly/${sname}_virmining
virfasta="minedviruses.fna"


####just identity calc with identity meshclust
if [ $method -eq 100 ];then
sed "s/ /____/g" $virfasta  > tmp.fna
/scratch/project_2007362/software/Identity-master/bin/identity -d tmp.fna -o ${sname}.idt -t 0 -a y -c $cpu

sed -i "s/____/ /g" ${sname}.idt
sed -i "s/>//g" ${sname}.idt


fi

############## pyani | not recommended, worth mushrooms
if [ $method -eq 0 ];then
conda activate pyani
	
sed "s/ /__/g" $virfasta | sed "s:/:____:g" > tmp.fna
#cat ~/phage_approach/ICTV_all.fasta >> tmp.fna

rm -rf tmppyani pyaniout && mkdir tmppyani
mv tmp.fna tmppyani
cd tmppyani
cat tmp.fna | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($1,2) ".fasta")}
        print $0 >> filename
        close(filename)
}'
rm -f tmp.fna
cd ..

ani="ANIb"
average_nucleotide_identity.py -i tmppyani -o pyaniout -m $ani --seed 1

sed -i "s:____:/:g" pyaniout/${ani}_percentage_identity.tab pyaniout/${ani}_alignment_coverage.tab
sed -i "s/__/ /g" pyaniout/${ani}_percentage_identity.tab pyaniout/${ani}_alignment_coverage.tab

mv pyaniout/${ani}_percentage_identity.tab pyaniout/${sname}_percentage_identity.tab
mv pyaniout/${ani}_alignment_coverage.tab pyaniout/${sname}_alignment_coverage.tab

rm -rf tmppyani


fi


############### MESHCLUST
if [ $method -eq 1 ];then
sed "s/ /__/g" $virfasta > tmp.fna

/scratch/project_2007362/software/Identity-master/bin/meshclust -d $virfasta -c $cpu -t 0.98 -o ${sname}_s.clstr
#awk -F"\t" '$4=="C"{gsub(">","");print $2}' ${sname}_s.clstr > representatives.txt
#seqkit grep -w 0 -nf representatives.txt $virfasta > sp95.fna

/scratch/project_2007362/software/Identity-master/bin/meshclust -d $virfasta -c $cpu -t 0.7 -o ${sname}_f.clstr
#awk -F"\t" '$4=="C"{gsub(">","");print $2}' ${sname}_f.clstr > representatives.txt
#seqkit grep -w 0 -nf representatives.txt $virfasta > f70.fna


fi

############# ALFATCLUST
if [ $method -eq 2 ];then
conda activate atclust   

sed "s/ /____/g" $virfasta > tmp.fna                                                 
#sed "s/ .*//g" ~/phage_approach/ICTV_all.fasta >> tmp.fna
virfasta=tmp.fna

$new/software/ALFATClust/main/alfatclust.py -i $virfasta -o ${sname}_s.atclust -e ${sname}_s.repatclust -t $cpu -b dna --seed 1 -l 0.94 -d 0.01 #max value 94
$new/software/ALFATClust/main/alfatclust.py -i $virfasta -o ${sname}_f.atclust -e ${sname}_f.repatclust -t $cpu -b dna --seed 1 -l 0.7
$new/software/ALFATClust/main/alfatclust.py -i $virfasta -o ${sname}_norank.atclust -e ${sname}_norank.repatclust -t $cpu --seed 1 -b dna -l 0.1

#awk 'BEGIN{RS="#Cluster____"}{for(i=2;i<NF;i++){print $1"\t"$i};print $1"\t"$NF}' ${sname}_s.tmp > ${sname}_s.atclust.clstr
#awk 'BEGIN{RS="#Cluster____"}{for(i=2;i<NF;i++){print $1"\t"$i};print $1"\t"$NF}' ${sname}_f.tmp > ${sname}_f.atclust.clstr

sed -i "s/____/ /g" ${sname}_s.atclust ${sname}_f.atclust ${sname}_s.repatclust ${sname}_f.repatclust

rm -f tmp.fna

fi
############## ANI NUCL
if [ $method -eq 3 ];then
sed "s/ /____/g" $virfasta > tmp.fna
virfasta=tmp.fna
ml blast

makeblastdb -in $virfasta -out blastdb -dbtype nucl 
blastn -query $virfasta -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -num_threads $cpu -max_target_seqs 50000 -perc_identity 10

python $new/software/aaicluster/blastani.py -i blast.tsv -o ani.tsv

python $new/software/aaicluster/cluster.py --fna $virfasta --ani ani.tsv --out ${sname}_ani_s.tsv --min_ani 98 --min_qcov 0 --min_tcov 85
python $new/software/aaicluster/cluster.py --fna $virfasta --ani ani.tsv --out ${sname}_ani_f.tsv --min_ani 70 --min_qcov 0 --min_tcov 85
python $new/software/aaicluster/cluster.py --fna $virfasta --ani ani.tsv --out ${sname}_ani_norank.tsv --min_ani 1 --min_qcov 0 --min_tcov 70


sed -i "s/____/ /g" ${sname}_ani_s.tsv ${sname}_ani_f.tsv ${sname}_ani_norank.tsv
#rm -f blast.tsv ani.tsv blastdb*

#exit

fi
############# VIRCLUST (PROT + PHYLO)
if [ $method -eq 4 ];then
BASE=$new/sandro/HeP_samples/2_assembly/${sname}_checkv
TMP=$BASE/virclust_${sname}_tmp

#rm -rf $TMP && mkdir $TMP
#export SINGULARITY_TMPDIR=$TMP
#export APPTAINER_TMPDIR=$TMP
#export NXF_SINGULARITY_CACHEDIR=$new/software/singularity/
#export SINGULARITY_CACHEDIR=$new/software/singularity/

#sed -i "s/ /__/g" clustered/sp95.fna
#cd clustered

#echo "working on $new/sandro/HeP_samples/2_assembly/${sname}_checkv"
#bash /scratch/project_2007362/software/virclustv2_singularity/virclust.bash projdir=$BASE/clustered infile=$(pwd)/$virfasta step1A=T step2A=T step3A=T step4A=T step5A=T cpu=4 boot_pv_a=yes bootstrap_no_a=100 clust_dist_a=0.925
#sed "s/____/ /g" clustered/04a-06a_genome_clustering_PC/05/virDF.tsv > ${sname}_virclust_family.tsv
#sed -i "s/___/-/g" ${sname}_virclust_family.tsv
#sed -i "s:__:/:g" ${sname}_virclust_family.tsv

#rm -rf $TMP tmp.fna
#exit

fi
############# ANI PROT (AAI)
if [ $method -eq 5 ];then
conda activate c5 # just an empty enviroment that contains numpy

AAI=/scratch/project_2007362/software/aaicluster
sed "s/ /____/g" $virfasta > tmp.fna            
virfasta=tmp.fna    

prodigal -a ${sname}.faa -f gff -g 11 -i $virfasta -m -n -o ${sname}_prodigal.tsv -p single
diamond makedb --in ${sname}.faa --db viral_proteins --threads $cpu
diamond blastp --query ${sname}.faa --db viral_proteins --out blastp.tsv --outfmt 6 --evalue 1e-5 --max-target-seqs 10000 --query-cover 50 --subject-cover 50
rm viral_proteins.dmnd

python $AAI/amino_acid_identity.py --in_faa ${sname}.faa --in_blast blastp.tsv --out_tsv aai.tsv
#Amino acid identity is computed based on the average BLAST percent identity between all genes shared between each pair of genomes (E-value <1e-5)

#Filter edges and prepare MCL input
python $AAI/filter_aai.py --in_aai aai.tsv --min_percent_shared 20 --min_num_shared 16 --min_aai 40 --out_tsv ${sname}_genus_edges.tsv
python $AAI/filter_aai.py --in_aai aai.tsv --min_percent_shared 10 --min_num_shared 8 --min_aai 20 --out_tsv ${sname}_family_edges.tsv

#Here we're keeping edges between genomes with >=20% AAI and genomes with either 8 shared genes or at least 20% of shared genes (relative to both genomes)

#Perform MCL-based clustering
mcl ${sname}_genus_edges.tsv -te $cpu -I 2.0 --abc -o ${sname}_genus_clusters.txt
mcl ${sname}_family_edges.tsv -te $cpu -I 1.2 --abc -o ${sname}_family_clusters.txt

#sed -i "s/____/ /g" ${sname}_genus_clusters.txt
sed -i "s/____/ /g" ${sname}_family_clusters.txt
rm -f blastp.tsv aai.tsv  tmp.fna

fi

##CDHIT
if [ $method -eq 6 ];then
sed "s:/:____:g" $virfasta > tmp.fna
sed -i "s/ /__/g" tmp.fna
virfasta=tmp.fna
ml cdhit

rm -rf tmpcdhit && mkdir tmpcdhit
cd-hit-est -i $virfasta -o tmpcdhit/cluster -c 0.95 -aS 0.85 -n 10 -M 5000 -d 0 -sc 1
sed "s:____:/:g" tmpcdhit/cluster.clstr > ${sname}_cdhit95.clstr
sed -i "s/__/ /g" ${sname}_cdhit95.clstr
rm -rf tmpcdhit

fi

#MANIAC
if [ $method -eq 7 ];then
	conda activate maniac
	sed "s:/:____:g" $virfasta > tmp.fna
	sed -i "s/ /__/g" tmp.fna
	virfasta=tmp.fna

	sed "s/IFILE/$virfasta/g" $new/software/MANIAC/template.yml > ${sname}.yml
	sed -i "s/OFILE/maniac/g" ${sname}.yml
	sed -i "s/MSIZE/5/g" ${sname}.yml
	rm -rf .snakemake maniac
	snakemake --cores $cpu --quiet --snakefile $new/software/MANIAC/MANIAC --configfile ${sname}.yml
	mv maniac/genome-alignment.csv maniac/${sname}_maniac.csv
	rm -rf .snakemake
fi
