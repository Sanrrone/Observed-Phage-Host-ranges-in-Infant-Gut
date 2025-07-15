#!/bin/bash
#
#SBATCH --job-name=testparallel
#SBATCH --output=test.txt
#SBATCH --error=test.err
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=5G
#SBATCH --time=03:00:00
#SBATCH --account=Project_2003252
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_207362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"


seq 1 1366 | xargs -P 10 -n 1 -I {} bash ./0_countreads.sh {} supp_files/samples.txt

exit
# creating output folders
mkdir -p 2_assembly
mkdir -p 3_phage_annotation
mkdir -p 4_hostpred
mkdir -p 5_phages_readmapback
mkdir -p 5_phages_readmapback/map
mkdir -p 5_phages_readmapback/nomap
mkdir -p 6_tax_profile
mkdir -p 6_tax_profile/phage
mkdir -p 6_tax_profile/bact
