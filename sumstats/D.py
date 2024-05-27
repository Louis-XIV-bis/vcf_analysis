import allel
import numpy as np
import json
import sys
from typing import Dict, Any, List

from functions import load_vcf, extract_genotype_data, save_to_json

# https://scikit-allel.readthedocs.io/en/stable/stats/diversity.html

def compute_D(genotypes: np.ndarray) -> float:
    """
    Compute and return population-wide Tajima's D for the whole population.

    Args:
        genotypes (np.ndarray): Genotype data.

    Returns:
        float: The computed value of Tajima's D statistic.
    Raises:
        RuntimeError: If the computation of Tajima's D statistic results in NaN or if any 
        other error occurs during the computation process.
    """

    try:
        allele_counts = allel.GenotypeArray(genotypes).count_alleles()
        D = allel.tajima_d(allele_counts)
        return D
    except Exception as e:
        raise RuntimeError(f"Error computing population D: {e}")

def main():
    
    if len(sys.argv) < 3:
        print("Usage: python script.py <vcf_file> <clade_file_dict>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    clade_dict = sys.argv[2]

    json_output_file = "tajimasD.json"

    callset = load_vcf(vcf_file)
    genotypes = extract_genotype_data(callset)
    
    results = {}

    # Compute population-wide Tajima's D
    D_pop = compute_D(genotypes)
    results['population'] = D_pop

    # Compute Tajima's D for each clade
    clade_dict = {1: ['AAAA','AAAD'], 2:['AAAB','AAAC']}
    for clade, sample_names in clade_dict.items():
        try:
            clade_genotypes = extract_genotype_data(callset, sample_names)
            D_clade = compute_D(clade_genotypes)
            results[clade] = D_clade
        except RuntimeError as e:
            print(f"Error computing Tajima's D for clade {clade}: {e}")

    # Save results to JSON
    save_to_json(results, json_output_file)

if __name__ == "__main__":
    main()
