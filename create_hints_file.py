fasta_input_path = '/Users/cankale/Desktop/gen_dataset/Apis_mellifera/apis_mellifera_exons.txt'  # Modify with the path to your input file
hints_output_path = '/Users/cankale/Desktop/gen_dataset/AUGUSTUS/Hints/mellifera_hints'  # Modify with the path to your output file

def process_fasta_to_hints(fasta_path, hints_path):
    try:
        with open(fasta_path, 'r') as fasta_file, open(hints_path, 'w') as hints_file:
            for line in fasta_file:
                if line.startswith('>'):
                    # Extract the necessary parts from the header
                    parts = line.strip().split('|')
                    seq_id = parts[1]  # This will get something like XM_396824.6
                    position_part = parts[2].split(':')[1]
                    position_part_last = position_part.split(' ')[0]
                    start, end = position_part_last.split('-')

                    # Define other fixed parts of the hints line
                    source = 'B'
                    feature = 'exonpart'
                    score = '0'  # Placeholder score, as actual scores are not provided
                    strand = '+'  # Placeholder strand, adjust as needed
                    phase = '.'
                    attributes = "src=M;pri=4"

                    # Write the formatted line to the hints file
                    hints_line = f"{seq_id}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n"
                    hints_file.write(hints_line)

        print(f"Hints file successfully created at: {hints_path}")
    except Exception as e:
        print(f"Error creating hints file: {e}")

# Run the function with the specified paths
process_fasta_to_hints(fasta_input_path, hints_output_path)
