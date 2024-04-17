import os
import pandas as pd
from Bio import SeqIO

def update_fasta(bed_dir, fasta_dir, output_dir):
    # Iterate over each BED file in the directory
    for bed_file in os.listdir(bed_dir):
        if bed_file.endswith(".bed"):
            bed_path = os.path.join(bed_dir, bed_file)
            bed_df = pd.read_csv(bed_path, header=None, sep="\t")
            bed_df.columns = ["chrom", "start", "end", "name"]

            # Parse the corresponding FASTA file
            fasta_file = os.path.join(fasta_dir, bed_file.replace(".bed", ".fa"))
            fasta_sequences = SeqIO.parse(fasta_file, "fasta")

            # Create a dictionary to store the sequences
            sequences = {}

            # Populate the dictionary with the sequences
            for fasta in fasta_sequences:
                sequences[fasta.id] = str(fasta.seq)

            # Create the output directory if it doesn't exist
            os.makedirs(output_dir, exist_ok=True)

            # Create a new file to write the updated sequences
            output_file = os.path.join(output_dir, bed_file.replace(".bed", "_updated.fasta"))
            with open(output_file, "w") as f:
                # Match the IDs and write the sequences to the new file
                for i, row in bed_df.iterrows():
                    if row["name"] in sequences:
                        f.write(f">{row['name']}\n{sequences[row['name']]}\n")
                    else:
                        print(f"Error: no sequence found for ID {row['name']}")

            print("Sequences updated and written to:", output_file)



update_fasta("/Users/cankale/Desktop/Genespace_data/bed", "/Users/cankale/Desktop/Genespace_data/peptide", "/Users/cankale/Desktop/Genespace_data/updated_fasta")
