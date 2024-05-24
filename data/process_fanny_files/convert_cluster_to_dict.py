import json

input_file="sample_cluster_name_90.txt"
output_file="../dict_cluster_strain.json"

# Initialize a dictionary to store cluster IDs and corresponding sample IDs
cluster_dict = {}

# Read the input file
with open(input_file, 'r') as file:
    # Skip the header
    next(file)
    for line in file:
        sample_id, cluster_id, _ = line.strip().split('\t')
        # If the cluster ID is not in the dictionary, add it with an empty list
        if cluster_id not in cluster_dict:
            cluster_dict[cluster_id] = [sample_id]
        else:
            # If the cluster ID is already in the dictionary, append the sample ID to its list
            cluster_dict[cluster_id].append(sample_id)

# Write the dictionary to a JSON file
with open(output_file, 'w') as json_file:
    json.dump(cluster_dict, json_file, indent=4)

print("Data has been written to", output_file)
