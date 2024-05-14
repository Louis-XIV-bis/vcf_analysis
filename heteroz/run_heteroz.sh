#!/bin/bash

#SBATCH --account yeast_neutral_model
#SBATCH --mem 250GB
#SBATCH --partition long

source activate /shared/ifbstor1/projects/yeast_neutral_model/envs

# CAREFUL: give vcf file w/o '.gz' in the path, it will add it afterwards
vcf='/shared/projects/yeast_neutral_model/vcf/2330strains_SNPs_filteredQD10_PASS_repeatMaskerGffJubin.vcf'

# Be careful: you should not give .gz to the python script
gunzip -c $vcf.gz > ./$vcf 
python3 heteroz.py $vcf
rm ./$vcf