import allel
import numpy as np
import json
import sys 
import re

from functions import load_vcf, extract_genotype_data, save_to_json

# https://scikit-allel.readthedocs.io/en/stable/stats/diversity.html

def compute_population_diversity(callset: dict, genotypes: np.ndarray) -> dict:
    """
    Compute and return population-wide genetic diversity (π) for each contig.

    Args:
        callset (dict): Callset containing VCF data.
        genotypes (np.ndarray): Genotype data.

    Returns:
        dict: A dictionary with contig names as keys and their respective genetic diversity (π) as values.

    Raises:
        RuntimeError: If there is an error computing population diversity.
    """
    try:
        # Extract contig names and variant positions from the callset
        contig_names = callset['variants/CHROM']
        variants_pos = callset['variants/POS']
        unique_contigs = np.unique(contig_names)

        pi_results = {}
        for contig in unique_contigs:
            # Create a mask to filter positions and genotypes specific to the current contig
            contig_mask = (contig_names == contig)
            contig_positions = variants_pos[contig_mask]
            contig_genotypes = genotypes[contig_mask]

            # Compute allele counts for the current contig
            allele_counts = allel.GenotypeArray(contig_genotypes).count_alleles()
            # Compute genetic diversity (π) for the current contig
            pi_pop = allel.sequence_diversity(contig_positions, allele_counts)
            pi_results[contig] = pi_pop

        return pi_results
    
    except Exception as e:
        raise RuntimeError(f"Error computing population diversity: {e}")

def compute_sample_diversity(callset: dict, genotypes: np.ndarray) -> dict:
    """
    Compute genetic diversity (π) for each sample per chromosome and return as a nested dictionary.

    Args:
        callset (dict): Callset containing VCF data.
        genotypes (np.ndarray): Genotype data.

    Returns:
        dict: Nested dictionary with sample IDs as keys and dictionaries of chromosome-specific π values as values.
    """
    diversity_dict = {}

    # Get sample IDs, contig names, and variant positions
    sample_ids = callset['samples']
    contig_names = callset['variants/CHROM']
    variants_pos = callset['variants/POS']
    unique_contigs = np.unique(contig_names)

    for i, sample_id in enumerate(sample_ids):
        sample_diversity = {}
        
        for contig in unique_contigs:
            try:
                # Create a mask for the current contig
                contig_mask = (contig_names == contig)
                contig_positions = variants_pos[contig_mask]
                # Extract genotypes for the current sample and contig
                contig_genotypes = genotypes[contig_mask, i, :]
                contig_genotypes = contig_genotypes[:, :, np.newaxis]

                # Compute allele counts for the current sample and contig
                allele_counts = allel.GenotypeArray(contig_genotypes).count_alleles()
                # Compute genetic diversity (π) for the current sample and contig
                pi_sample = allel.sequence_diversity(contig_positions, allele_counts)
                sample_diversity[contig] = pi_sample

            except Exception as e:
                print(f"Error computing diversity for sample {sample_id} on contig {contig}: {e}")

        diversity_dict[sample_id] = sample_diversity

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
