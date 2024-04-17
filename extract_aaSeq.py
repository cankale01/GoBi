import sys
import pandas as pd
from Bio import Entrez


Entrez.email = "pasichye@gmail.com"


csv_file = sys.argv[1]
gff_file = sys.argv[2]
output_file = sys.argv[3]

blast_results = pd.read_csv(csv_file, header=None, usecols=[1, 8, 9], names=['id', 'start', 'end'], delimiter='\t')


cds_info = []
with open(gff_file, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split('\t')
        if fields[2] == 'CDS':
            seq_id = fields[0]
            if seq_id in blast_results['id'].values:
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]
                gene_id = None
                for attr in attributes.split(';'):
                    if attr.startswith('ID=cds-'):
                        gene_id = attr.split('-')[1]
                        break
                if gene_id:
                    cds_info.append((seq_id, start, end, gene_id))

tolerance = 10
ordered_gene_ids = []
for _, blast_row in blast_results.iterrows():
    for seq_id, start, end, gene_id in cds_info:
        start_match = start - tolerance <= blast_row['start'] <= end + tolerance
        end_match = start - tolerance <= blast_row['end'] <= end + tolerance

        if seq_id == blast_row['id'] and (start_match or end_match):
            if gene_id not in ordered_gene_ids:
                ordered_gene_ids.append(gene_id)
                break



sequences = {}
for rna_id in ordered_gene_ids:
    try:
        handle = Entrez.efetch(db="protein", id=rna_id, rettype="fasta", retmode="text")
        seq_record = handle.read()
        sequences[rna_id] = seq_record
        handle.close()
    except Exception as e:
        print(f"Error fetching sequence for ID {rna_id}: {e}")





if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <csv file> <gff file> <output>")
        sys.exit(1)
    else:
        with open(output_file, 'w') as out:
            for s in sequences.values():
                out.write(s)
        print(f'{output_file} successfully created')
