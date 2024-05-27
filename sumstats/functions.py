import allel
import numpy as np
import json
from typing import Any, Dict, List, Optional


def load_vcf(file_path: str) -> Dict[str, Any]:
    """
    Load a VCF (Variant Call Format) file and return its contents as a dictionary.

    Parameters:
    file_path (str): Path to the VCF file.

    Returns:
    Dict[str, Any]: Dictionary containing the VCF data.

    Raises:
    IOError: If there is an error loading the VCF file.
    """
    try:
        callset = allel.read_vcf(file_path)
        return callset
    except Exception as e:
        raise IOError(f"Error loading VCF file: {e}")
    
def extract_genotype_data(callset: Dict[str, Any], sample_names: Optional[List[str]] = None) -> np.ndarray:
    """
    Extract genotype data from a callset dictionary optionally for specified sample names.

    Parameters:
    callset (Dict[str, Any]): Dictionary containing VCF data.
    sample_names (Optional[List[str]]): Optional list of sample names to filter genotype data.

    Returns:
    np.ndarray: Numpy array containing the genotype data. If sample_names is provided, genotypes for specified samples are returned, otherwise, all genotypes are returned.

    Raises:
    KeyError: If genotype data or sample names are not found in the callset.
    """
    try:
        genotypes = callset['calldata/GT']
        if genotypes.dtype != 'i1':
            genotypes = genotypes.astype('i1')
        
        if sample_names is None:
            return genotypes
        
        # Extract sample indices
        all_sample_names = callset['samples']
        sample_set = set(sample_names)
        sample_indices = [i for i, name in enumerate(all_sample_names) if name in sample_set]
        
        # Filter genotype data for the specified samples
        filtered_genotypes = genotypes[:, sample_indices, :]
        
        return filtered_genotypes
    except KeyError as e:
        raise KeyError(f"Required data not found in callset: {e}")
    except ValueError as e:
        raise ValueError(f"One or more sample names not found in callset: {e}")
    
    
def save_to_json(data: Dict[str, Any], file_path: str) -> None:
    """
    Save dictionary data to a JSON file.

    Args:
        data (dict): Data to save.
        file_path (str): Path to the JSON file.
    """
    with open(file_path, 'w') as json_file:
        json.dump(data, json_file, indent=4)

def load_json_to_dict(file_path: str) -> Optional[Dict[str, Any]]:
    """
    Loads a JSON file into a dictionary.
    
    Args:
    - file_path (str): The path to the JSON file.
    
    Returns:
    - Optional[Dict[str, Any]]: The contents of the JSON file as a dictionary,
                                or None if the file is not found or invalid.
    """
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
            return data
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        return None
    except json.JSONDecodeError:
        print(f"Error: The file {file_path} is not a valid JSON.")
        return None