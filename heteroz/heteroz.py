#!/usr/bin/env python3

import json
import vcf
import sys
from typing import Dict

def count_genotypes(vcf_file: str) -> Dict[str, Dict[str, int]]:
    """
    Count the number of homozygous and heterozygous genotypes for each sample in a VCF file.
    Do not count the set which contains 0/0 or ./. for instance. 

    Args:
        vcf_file (str): Path to the VCF file.

    Returns:
        Dict[str, Dict[str, int]]: A dictionary containing the counts of homozygous ('hom') 
        and heterozygous ('het') genotypes for each sample. Returns None if the file is not found.
    """
    sample_results = {}

    try:
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    except FileNotFoundError:
        print("File not found.")
        return None

    for record in vcf_reader:
        for sample in record.samples:
            sample_name = sample.sample

            if sample_name not in sample_results:
                sample_results[sample_name] = {'hom': 0, 'het': 0}

            genotype = sample['GT']

            # Skip sites with no SNP / indel (haplo to polyploids): can be 0 or ., a mix of both 
            if genotype in ['0', '.']:
                continue
            else:
                genotype_parts = genotype.split('|') if '|' in genotype else genotype.split('/')
                if all(part in ['0', '.'] for part in genotype_parts): #if contains ., 0 or both only => skip
                    continue

            alleles = genotype.split('|' if '|' in genotype else '/')
            if alleles[0] == alleles[1]:
                sample_results[sample_name]['hom'] += 1
            else:
                sample_results[sample_name]['het'] += 1

    return sample_results

def save_results_to_file(sample_results: Dict[str, Dict[str, int]], output_file: str) -> None:
    """
    Save the genotype counts for each sample to a JSON file.

    Args:
        sample_results (Dict[str, Dict[str, int]]): A dictionary containing the counts of homozygous ('hom') 
        and heterozygous ('het') genotypes for each sample.
        output_file (str): Path to the output JSON file.
    """
    try:
        with open(output_file, 'w') as file:
            json.dump(sample_results, file, indent=4)
        print("Results saved successfully.")
    except Exception as e:
        print(f"Error occurred while saving results: {e}")

def main():
    vcf_file = f'{sys.argv[1]}'
    output_file = 'sample_homhet.json'

    results = count_genotypes(vcf_file)
    if results:
        save_results_to_file(results, output_file)

if __name__ == "__main__":
    main()
