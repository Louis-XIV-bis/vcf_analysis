import allel
import numpy as np
import json
import sys 

def load_vcf(file_path: str) -> dict:
    """
    Load VCF file and return the callset.

    Args:
        file_path (str): Path to the VCF file.

    Returns:
        dict: Callset containing VCF data.

    Raises:
        IOError: If there is an error loading the VCF file.
    """
    try:
        callset = allel.read_vcf(file_path)
        return callset
    
    except Exception as e:
        raise IOError(f"Error loading VCF file: {e}")

def extract_genotype_data(callset: dict) -> np.ndarray:
    """
    Extract and return genotype data from the callset.

    Args:
        callset (dict): Callset containing VCF data.

    Returns:
        np.ndarray: Genotype data.

    Raises:
        KeyError: If genotype data is not found in the callset.
    """
    try:
        genotypes = callset['calldata/GT']
        if genotypes.dtype != 'i1':
            genotypes = genotypes.astype('i1')
        return genotypes
    
    except KeyError as e:
        raise KeyError(f"Genotype data not found in callset: {e}")

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
        allele_counts = allel.GenotypeArray(genotypes).count_alleles()
        pi_pop = allel.sequence_diversity(callset['variants/POS'], allele_counts)
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

def save_to_json(data: dict, file_path: str) -> None:
    """
    Save dictionary data to a JSON file.

    Args:
        data (dict): Data to save.
        file_path (str): Path to the JSON file.
    """
    with open(file_path, 'w') as json_file:
        json.dump(data, json_file, indent=4)

def main():
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
