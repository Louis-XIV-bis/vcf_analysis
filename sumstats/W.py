import allel
import numpy as np
import json
import sys 

from functions import load_vcf, extract_genotype_data, save_to_json

# https://scikit-allel.readthedocs.io/en/stable/stats/diversity.html

def compute_population_W(callset: dict, genotypes: np.ndarray) -> float:
    """
    Compute and return population-wide Watterson’s estimator (W).

    Args:
        callset (dict): Callset containing VCF data.
        genotypes (np.ndarray): Genotype data.

    Returns:
        float: Population-wide Watterson’s estimator (W).

    Raises:
        RuntimeError: If there is an error computing population Watterson’s estimator (W).
    """
    try:
        variants_pos = callset['variants/POS']
        allele_counts = allel.GenotypeArray(genotypes).count_alleles()
        W_pop = allel.watterson_theta(variants_pos, allele_counts)
        return W_pop
    
    except Exception as e:
        raise RuntimeError(f"Error computing population diversity: {e}")

def compute_sample_W(callset: dict, genotypes: np.ndarray) -> dict:
    """
    Compute Watterson’s estimator (W) for each sample and return as a dictionary.

    Args:
        callset (dict): Callset containing VCF data.
        genotypes (np.ndarray): Genotype data.

    Returns:
        dict: Dictionary with sample IDs as keys and their Watterson’s estimator (W) as values.
    """

    W_dict = {}

    # Get sample IDs & SNP pos
    sample_ids = callset['samples']
    variants_pos = callset['variants/POS']

    for i, sample_id in enumerate(sample_ids):
        try:
            # Extract genotypes for the current sample & add a singleton dimension
            # because it's needed for the tool (the removed dimension is the sample)
            sample_genotypes = genotypes[:, i, :]
            sample_genotypes = sample_genotypes[:, :, np.newaxis]
            
            # Compute allele counts for the current sample
            allele_counts = allel.GenotypeArray(sample_genotypes).count_alleles()
            W_sample = allel.watterson_theta(variants_pos, allele_counts)
            W_dict[sample_id] = W_sample

        except Exception as e:
            print(f"Error computing Watterson’s estimator for sample {sample_id}: {e}")

    return W_dict

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <vcf_file>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    json_output_file = "W.json"

    callset = load_vcf(vcf_file)
    genotypes = extract_genotype_data(callset)
    
    results = {}

    # Compute population-wide Watterson’s estimator (W)
    pi_pop = compute_population_W(callset, genotypes)
    results['population'] = pi_pop

    # Compute sample-specific Watterson’s estimator (W)
    sample_diversity = compute_sample_W(callset, genotypes)
    results.update(sample_diversity)

    # Save results to JSON
    save_to_json(results, json_output_file)

if __name__ == "__main__":
    main()
