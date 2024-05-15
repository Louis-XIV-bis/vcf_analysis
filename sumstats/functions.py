import allel
import numpy as np
import json

def load_vcf(file_path: str) -> dict:
    try:
        callset = allel.read_vcf(file_path)
        return callset
    except Exception as e:
        raise IOError(f"Error loading VCF file: {e}")

def extract_genotype_data(callset: dict) -> np.ndarray:
    try:
        genotypes = callset['calldata/GT']
        if genotypes.dtype != 'i1':
            genotypes = genotypes.astype('i1')
        return genotypes
    except KeyError as e:
        raise KeyError(f"Genotype data not found in callset: {e}")
    
def save_to_json(data: dict, file_path: str) -> None:
    """
    Save dictionary data to a JSON file.

    Args:
        data (dict): Data to save.
        file_path (str): Path to the JSON file.
    """
    with open(file_path, 'w') as json_file:
        json.dump(data, json_file, indent=4)