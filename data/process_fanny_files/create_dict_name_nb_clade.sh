#!/bin/bash

input_file="cluster_names_90.txt"
output_file="../dict_clade_nb.json"

echo "{" > $output_file

index=1
while IFS= read -r line; do
    echo -e "\t\"$index\": \"$line\"," >> $output_file
    ((index++))
done < $input_file

# Remove the last comma and add the closing brace
sed -i '$ s/,$//' $output_file
echo "}" >> $output_file

echo "Data has been written to $output_file"
