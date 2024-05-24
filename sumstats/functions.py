import allel
import numpy as np
import json
from typing import Any, Dict, Optional


def load_vcf(file_path: str) -> Dict[str, Any]:
    try:
        callset = allel.read_vcf(file_path)
        return callset
    except Exception as e:
        raise IOError(f"Error loading VCF file: {e}")

def extract_genotype_data(callset: Dict[str, Any]) -> np.ndarray:
    try:
        genotypes = callset['calldata/GT']
        if genotypes.dtype != 'i1':
            genotypes = genotypes.astype('i1')
        return genotypes
    except KeyError as e:
        raise KeyError(f"Genotype data not found in callset: {e}")
    
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