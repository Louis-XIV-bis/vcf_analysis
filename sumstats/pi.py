import allel
import numpy as np
import sys 
from typing import Dict, List

from functions import load_vcf, extract_genotype_data, save_to_json, load_json_to_dict

# https://scikit-allel.readthedocs.io/en/stable/stats/diversity.html

def compute_population_diversity(callset: Dict[str, np.ndarray], genotypes: np.ndarray) -> dict:
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

def compute_clade_diversity(callset: Dict[str, np.ndarray], genotypes: np.ndarray, clusters: Dict[int, List[str]]) -> Dict[int, Dict[str, float]]:
    """
    Compute and return population-wide genetic diversity (π) for each contig for each cluster of samples.

    Args:
        callset (dict): Callset containing VCF data.
        genotypes (np.ndarray): Genotype data.
        clusters (dict): A dictionary where keys are cluster numbers and values are lists of sample names.

    Returns:
        dict: A dictionary with cluster numbers as keys and another dictionary as value.
              The inner dictionary has contig names as keys and their respective genetic diversity (π) as values.

    Raises:
        RuntimeError: If there is an error computing population diversity.
        ValueError: If a sample name in the clusters is not found in the callset samples.
    """
    try:
        # Extract contig names and variant positions from the callset
        contig_names = callset['variants/CHROM']
        variants_pos = callset['variants/POS']
        unique_contigs = np.unique(contig_names)
        sample_names = callset['samples']

        # Check if all sample names in clusters exist in the callset sample names
        for cluster, cluster_samples in clusters.items():
            for sample in cluster_samples:
                if sample not in sample_names:
                    raise ValueError(f"Sample name '{sample}' in cluster '{cluster}' is not found in the callset samples.")

        cluster_results = {}

        for cluster, cluster_samples in clusters.items():
            # Find indices of the samples that belong to the current cluster
            cluster_sample_indices = [sample_names.tolist().index(sample) for sample in cluster_samples]

            pi_results = {}
            for contig in unique_contigs:
                # Create a mask to filter positions and genotypes specific to the current contig
                contig_mask = (contig_names == contig)
                contig_positions = variants_pos[contig_mask]
                contig_genotypes = genotypes[contig_mask][:, cluster_sample_indices]

                # Compute allele counts for the current contig
                allele_counts = allel.GenotypeArray(contig_genotypes).count_alleles()

                # Compute genetic diversity (π) for the current contig
                pi_pop = allel.sequence_diversity(contig_positions, allele_counts)
                pi_results[contig] = pi_pop

            cluster_results[cluster] = pi_results

        return cluster_results

    except ValueError as ve:
        raise ve
    except Exception as e:
        raise RuntimeError(f"Error computing population diversity: {e}")
    
def add_genome_wide_pi(data: Dict[str, Dict[str, float]], weights: Dict[str, float]) -> Dict[str, Dict[str, float]]:
    """
    Adds a 'genome-wide' key to each entry in the provided dictionary, with the value being the
    weighted mean of the other values in the entry. The weights are related to the chromosome size.
    
    Parameters:
    data (dict): A dictionary where each key maps to another dictionary of values.
    weights (dict): A dictionary where keys are the same as those in the inner dictionaries of `data`,
                    and values are the weights associated with each key.
    
    Returns:
    dict: The input dictionary with the 'genome-wide' key added to each entry.
    """

    # Loop through each entry in the original data dictionary
    for key, values in data.items():
        weighted_sum = 0
        total_weight = 0
        for k, v in values.items():
            if k in weights:
                weight = weights[k]
                weighted_sum += v * weight
                total_weight += weight
        if total_weight > 0:
            genome_wide_mean = weighted_sum / total_weight
        else:
            genome_wide_mean = 0  # or handle case where there are no valid weights
        data[key]['genome-wide'] = genome_wide_mean

    return data

def main():
    
    if len(sys.argv) < 4:
        print("Usage: python script.py <vcf_file> <clade_file_dcit> <chromosome_size_dict>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    clade_dict = sys.argv[2]
    chr_dict = sys.argv[3]

    json_output_file = "diversity.json"

    # Load the file and extract the genotypes
    callset = load_vcf(vcf_file)
    genotypes = extract_genotype_data(callset)
    
    results = {}

    # Compute population-wide diversity
    pi_pop = compute_population_diversity(callset, genotypes)
    results['population'] = pi_pop

    # Compute clade-specific diversity
    clades = load_json_to_dict(clade_dict)
    # clades = {1: ['AAAA','AAAD'], 2:['AAAB','AAAC']}

    sample_diversity = compute_clade_diversity(callset, genotypes, clades)
    results.update(sample_diversity)

    # Compute genome-wide pi for each clade (weigthed mean of chromosome pi)
    chr_size = load_json_to_dict(chr_dict)
    results = add_genome_wide_pi(results, chr_size)
    print(results)

    # Save results to JSON
    save_to_json(results, json_output_file)

if __name__ == "__main__":
    main()
