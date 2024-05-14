import gzip
import vcf 

# # Specify the path to your gz file
gz_file_path = 'head_vcf.vcf.gz'

# Open the gz file in binary mode ('rb')
with gzip.open(gz_file_path, 'rb') as gz_file:
    # Read the contents of the gz file
    gz_content = gz_file.read()

# Decode the content if necessary (assuming it contains text)
decoded_content = gz_content.decode('utf-8')

# Now you can work with the content as needed
# print(decoded_content)

# Write the decoded content to a new file
# with open('test.cf', 'w') as output_file:
#     output_file.write(decoded_content)
# decoded_content = 'test.vcf'

vcf_reader = vcf.Reader(open(decoded_content))
for record in vcf_reader:
    for sample in record.samples:
        sample_name = sample.sample
        print(sample_name)



