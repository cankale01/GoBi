from Bio import SeqIO
import pandas as pd

# Parse the FASTA file
proteins = SeqIO.parse('/Users/cankale/Desktop/Genespace_data/peptide/proteins_tribolium_castaneum.fa', 'fasta')

# Create an empty DataFrame to store the BED data
bed_data = pd.DataFrame(columns=['chrom', 'start', 'end', 'name'])

# Iterate over each protein
for protein in proteins:
    # Get the GeneID from the protein description
    gene_id = protein.description.split('[GeneID=')[-1].split(']')[0]

    # Search for the GeneID in the GFF file
    with open('/Users/cankale/Desktop/gen_dataset/Tribolium_castaneum/Tribolium_castaneum.gff') as gff_file:
        for line in gff_file:
            if gene_id in line:
                # Extract the chrom, start, and end from the GFF line
                chrom, _, _, start, end, _, _, _, _ = line.split('\t')
                start = int(start) - 1  # BED files are 0-based
                end = int(end)

                # Add the data to the DataFrame
                new_row = pd.DataFrame([{'chrom': chrom, 'start': start, 'end': end, 'name': protein.id}])
                bed_data = pd.concat([bed_data, new_row], ignore_index=True)
                break

# Write the DataFrame to a BED file
bed_data.to_csv('/Users/cankale/Desktop/Genespace_data/bed/tribolium_castaneum.bed', sep='\t', index=False)