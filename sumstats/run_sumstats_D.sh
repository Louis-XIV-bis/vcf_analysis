#!/bin/bash

#SBATCH --account yeast_neutral_model
#SBATCH --mem 250GB
#SBATCH --partition long

conda activate /shared/ifbstor1/projects/yeast_neutral_model/envs/

# vcf_path='/shared/projects/yeast_neutral_model/vcf/2330strains_SNPs_filteredQD10_PASS_repeatMaskerGffJubin.vcf.gz'
vcf='/shared/projects/yeast_neutral_model/vcf/data/vcf_fanny.vcf'
#vcf='/shared/projects/yeast_neutral_model/vcf/data/test.vcf'
clade='/shared/projects/yeast_neutral_model/vcf/data/dict_cluster_strain.json'
chr_size='/shared/projects/yeast_neutral_model/vcf/data/chromsome_sizes.json'

# Be careful: you should not give .gz to the python script
# gunzip -c $vcf_path > $vcf 

# python3 pi.py $vcf $clade $chr_size
python3 D.py $vcf $clade $chr_size
# python3 W.py $vcf
# python3 het_variant.py $vcf
# other sumstats 

# rm $vcf