import ssl

from Bio import Entrez, SeqIO

# Always provide your email when using NCBI E-utilities
Entrez.email = "your.email@example.com"
ssl._create_default_https_context = ssl._create_unverified_context

# List of accession numbers for which to fetch CDS sequences
accession_numbers = [
    "XM_396824.6", "XM_006558943", "XM_003249378.4", "XM_006558928.3",
    "NM_001011601.1", "NM_001011579.1", "XM_026443532.1", "NM_001011622.1",
    "XM_026443531.1", "XM_026443530.1", "NM_001014429.1", "NM_001011564.2",
    "NM_001024697.1", "XM_006558929.3"
]

# Initialize a list to hold all CDS records
all_cds_records = []

for accession in accession_numbers:
    # Fetching the sequence
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta_cds_na", retmode="text")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    # Add the fetched records to the all_cds_records list
    all_cds_records.extend(records)

# Saving all CDS sequences to a single file
filename = "all_cds_sequences.fasta"
with open(filename, "w") as output_handle:
    SeqIO.write(all_cds_records, output_handle, "fasta")
print(f"Saved all CDS sequences in {filename}")

# Note: This code fetches the nucleotide sequences for the CDS regions and saves them in a single FASTA file.
