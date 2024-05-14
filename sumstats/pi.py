import allel
import numpy as np
import json
import sys 

from functions import load_vcf, extract_genotype_data, save_to_json

# https://scikit-allel.readthedocs.io/en/stable/stats/diversity.html

def compute_population_diversity(callset: dict, genotypes: np.ndarray) -> float:
    """
    Compute and return population-wide genetic diversity (π).

    Args:
        callset (dict): Callset containing VCF data.
        genotypes (np.ndarray): Genotype data.

    Returns:
        float: Population-wide genetic diversity (π).

    Raises:
        RuntimeError: If there is an error computing population diversity.
    """
    try:
        variants_pos = callset['variants/POS']
        allele_counts = allel.GenotypeArray(genotypes).count_alleles()
        pi_pop = allel.sequence_diversity(variants_pos, allele_counts)
        return pi_pop
    
    except Exception as e:
        raise RuntimeError(f"Error computing population diversity: {e}")

def compute_sample_diversity(callset: dict, genotypes: np.ndarray) -> dict:
    """
    Compute genetic diversity (π) for each sample and return as a dictionary.

    Args:
        callset (dict): Callset containing VCF data.
        genotypes (np.ndarray): Genotype data.

    Returns:
        dict: Dictionary with sample IDs as keys and their genetic diversity (π) as values.
    """

    diversity_dict = {}

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
            pi_sample = allel.sequence_diversity(variants_pos, allele_counts)
            diversity_dict[sample_id] = pi_sample

        except Exception as e:
            print(f"Error computing diversity for sample {sample_id}: {e}")

    return diversity_dict

def main():
    
    if len(sys.argv) < 2:
        print("Usage: python script.py <vcf_file>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    json_output_file = "diversity.json"

    callset = load_vcf(vcf_file)
    genotypes = extract_genotype_data(callset)
    
    results = {}

    # Compute population-wide diversity
    pi_pop = compute_population_diversity(callset, genotypes)
    results['population'] = pi_pop

    # Compute sample-specific diversity
    sample_diversity = compute_sample_diversity(callset, genotypes)
    results.update(sample_diversity)

    # Save results to JSON
    save_to_json(results, json_output_file)

if __name__ == "__main__":
    main()
