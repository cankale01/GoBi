from Bio import SeqIO

# Path to your genome file
genome_file_path = "/Users/cankale/Desktop/gen_dataset/apis_florea/Apis_florea_genomic.fna"

# The ID of the sequence you're interested in
sequence_id = "NW_003790927.1"

start_position = 302732
end_position = 302954

# Flag to indicate if the sequence was found
sequence_found = False

# Iterate over each record in the fasta file
for record in SeqIO.parse(genome_file_path, "fasta"):
    if record.id == sequence_id:
        # Extract the region of interest
        region_of_interest = record.seq[start_position - 1:end_position]
        print(f"Extracted region from {sequence_id}: {region_of_interest}")
        sequence_found = True

        # Optionally, save this region to a new file
        with open("extracted_region.fasta", "w") as output_handle:
            output_handle.write(f">{sequence_id}_{start_position}_{end_position}\n")
            output_handle.write(str(region_of_interest) + "\n")

        break  # Exit the loop once the sequence is found and processed

if not sequence_found:
    print(f"Sequence ID {sequence_id} not found in the file.")

