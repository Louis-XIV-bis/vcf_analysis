import json
import pandas as pd
from typing import Dict, List, Any

def load_json(file_path: str) -> Dict[str, Any]:
    """
    Load JSON data from a file and return it as a dictionary.

    Parameters:
    file_path (str): The path to the JSON file.

    Returns:
    dict: The JSON data as a dictionary.
    """
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def json_to_dataframe(data: Dict[str, Any], column_names: List[str] = ['ID', 'hom', 'het']) -> pd.DataFrame:
    """
    Convert a dictionary to a pandas DataFrame with specific columns.

    Parameters:
    data (dict): The input dictionary with nested dictionaries.
    column_names (List[str]): The list of column names for the DataFrame. Default to ['ID', 'hom', 'het'].

    Returns:
    pd.DataFrame: The resulting DataFrame with specified columns.
    """
    df = pd.DataFrame.from_dict(data, orient='index').reset_index()
    df.columns = column_names
    return df

def main():
    file_path = 'sample_homhet.json'

    data = load_json(file_path)

    df = json_to_dataframe(data)
    del data

if __name__ == "__main__":
    main()
