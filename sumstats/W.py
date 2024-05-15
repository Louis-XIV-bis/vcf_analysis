import allel
import numpy as np
import json
import sys 

from functions import load_vcf, extract_genotype_data, save_to_json

# https://scikit-allel.readthedocs.io/en/stable/stats/diversity.html

import allel
import numpy as np

def compute_population_W(callset: dict, genotypes: np.ndarray) -> dict:
    """
    Compute and return population-wide Watterson’s estimator (W) for each contig.

    Args:
        callset (dict): Callset containing VCF data.
        genotypes (np.ndarray): Genotype data.

    Returns:
        dict: A dictionary with contig names as keys and their respective Watterson’s estimator (W) as values.

    Raises:
        RuntimeError: If there is an error computing population Watterson’s estimator (W).
    """
    try:
        contig_names = callset['variants/CHROM']
        variants_pos = callset['variants/POS']
        unique_contigs = np.unique(contig_names)

        W_results = {}
        for contig in unique_contigs:
            # Create a mask to filter positions and genotypes specific to the current contig
            contig_mask = (contig_names == contig)
            contig_positions = variants_pos[contig_mask]
            contig_genotypes = genotypes[contig_mask]

            # Compute allele counts for the current contig
            allele_counts = allel.GenotypeArray(contig_genotypes).count_alleles()
            # Compute Watterson’s estimator (W) for the current contig
            W_pop = allel.watterson_theta(contig_positions, allele_counts)
            W_results[contig] = W_pop

        return W_results
    
    except Exception as e:
        raise RuntimeError(f"Error computing population Watterson’s estimator: {e}")

def compute_sample_W(callset: dict, genotypes: np.ndarray) -> dict:
    """
    Compute Watterson’s estimator (W) for each sample per contig and return as a nested dictionary.

    Args:
        callset (dict): Callset containing VCF data.
        genotypes (np.ndarray): Genotype data.

    Returns:
        dict: Nested dictionary with sample IDs as keys and dictionaries of contig-specific Watterson’s estimator (W) values.
    """
    W_dict = {}

    # Get sample IDs, contig names, and variant positions
    sample_ids = callset['samples']
    contig_names = callset['variants/CHROM']
    variants_pos = callset['variants/POS']
    unique_contigs = np.unique(contig_names)

    for i, sample_id in enumerate(sample_ids):
        sample_W = {}
        
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
                # Compute Watterson’s estimator (W) for the current sample and contig
                W_sample = allel.watterson_theta(contig_positions, allele_counts)
                sample_W[contig] = W_sample

            except Exception as e:
                print(f"Error computing Watterson’s estimator for sample {sample_id} on contig {contig}: {e}")

        W_dict[sample_id] = sample_W

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
