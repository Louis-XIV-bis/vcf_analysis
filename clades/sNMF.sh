#!/bin/bash

#SBATCH --account yeast_neutral_model
#SBATCH --mem 250GB
#SBATCH --partition long

conda activate /shared/ifbstor1/projects/yeast_neutral_model/envs/

# vcf_path='/shared/projects/yeast_neutral_model/vcf/2330strains_SNPs_filteredQD10_PASS_repeatMaskerGffJubin.vcf.gz'
vcf='/shared/projects/yeast_neutral_model/vcf/vcf_fanny.vcf'

# gunzip -c $vcf_path > $vcf 




# rm $vcf