import allel
import numpy as np
import json
import sys

from functions import load_vcf, extract_genotype_data, save_to_json

# https://scikit-allel.readthedocs.io/en/stable/stats/diversity.html

def compute_population_D(callset: dict, genotypes: np.ndarray) -> float:
    try:
        allele_counts = allel.GenotypeArray(genotypes).count_alleles()
        pi_pop = allel.tajima_d(allele_counts)
        if np.isnan(pi_pop):
            raise RuntimeError("Tajima's D resulted in NaN.")
        return pi_pop
    except Exception as e:
        raise RuntimeError(f"Error computing population D: {e}")

def main():
    
    if len(sys.argv) < 2:
        print("Usage: python script.py <vcf_file>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    json_output_file = "tajimasD.json"

    callset = load_vcf(vcf_file)
    genotypes = extract_genotype_data(callset)
    
    results = {}

    # Compute population-wide Tajima's D
    try:
        pi_pop = compute_population_D(callset, genotypes)
        results['population'] = pi_pop
    except RuntimeError as e:
        print(e)

    # Save results to JSON
    save_to_json(results, json_output_file)

if __name__ == "__main__":
    main()
