import allel
import numpy as np
import json
import sys
from typing import Dict, List, Any

from functions import load_vcf, extract_genotype_data, save_to_json

def compute_obs_het_variant(genotypes: np.ndarray) -> np.ndarray:
    """
    Compute the observed heterozygosity for each variant.

    Parameters:
    genotypes (np.ndarray): Array of genotypes.

    Returns:
    np.ndarray: Array of observed heterozygosity values for each variant.
    """
    g = allel.GenotypeArray(genotypes)
    het_obs = allel.heterozygosity_observed(g)
    return het_obs

def compute_HW_het_variant(genotypes: np.ndarray, ploidy=2) -> np.ndarray:
    """
    Compute the expected heterozygosity under Hardy-Weinberg equilibrium for each variant.
    WARNING : assume ploidy = 2 ==> run for only haplo / polyploids & new ploidy if wanted

    Parameters:
    genotypes (np.ndarray): Array of genotypes.
    ploidy (int): ploidy used to compute the expected frequencies (Default to 2).

    Returns:
    np.ndarray: Array of expected heterozygosity values for each variant.
    """
    g = allel.GenotypeArray(genotypes)
    af = g.count_alleles().to_frequencies()
    
    # WARNING : assume ploidy = 2 ==> run for only haplo / polyploids & new ploidy if wanted
    het_hw = allel.heterozygosity_expected(af, ploidy=ploidy)
    return het_hw

def compute_inbreed_coef_variant(genotypes: np.ndarray) -> np.ndarray:
    """
    Compute the inbreeding coefficient for each variant.

    Parameters:
    genotypes (np.ndarray): Array of genotypes.

    Returns:
    np.ndarray: Array of inbreeding coefficients for each variant.
    """
    g = allel.GenotypeArray(genotypes)
    inb_coef = allel.inbreeding_coefficient(g)
    return inb_coef

def aggregate_results(obs_het: np.ndarray, HW_het: np.ndarray, inb_coef: np.ndarray, sample_ids: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Aggregate the observed heterozygosity, expected heterozygosity, and inbreeding coefficient for each sample.

    Parameters:
    obs_het (np.ndarray): Array of observed heterozygosity values.
    HW_het (np.ndarray): Array of expected heterozygosity values.
    inb_coef (np.ndarray): Array of inbreeding coefficients.
    sample_ids (List[str]): List of sample IDs.

    Returns:
    Dict[str, Dict[str, Any]]: Dictionary mapping each sample ID to its observed heterozygosity, expected heterozygosity, and inbreeding coefficient.
    """
    results = {}
    for i, sample_id in enumerate(sample_ids):
        results[sample_id] = {
            'observed_het': obs_het[i],
            'HW_het': HW_het[i],
            'inbreeding_coef': inb_coef[i]
        }
    return results

def main():

    if len(sys.argv) < 2:
        print("Usage: python script.py <vcf_file>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    json_output_file = "het_HW.json"

    callset = load_vcf(vcf_file)
    genotypes = extract_genotype_data(callset)
    sample_ids = callset['samples']

    # Compute statistics for each variants
    obs_het = compute_obs_het_variant(genotypes)
    HW_het = compute_HW_het_variant(genotypes)
    inb_coef = compute_inbreed_coef_variant(genotypes)

    # Aggregate results
    results = aggregate_results(obs_het, HW_het, inb_coef, sample_ids)

    # Save results to JSON
    save_to_json(results, json_output_file)

if __name__ == "__main__":
    main()
