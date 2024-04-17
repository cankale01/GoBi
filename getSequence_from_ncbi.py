import ssl
import sys
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"
ssl._create_default_https_context = ssl._create_unverified_context


def read_ids(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]


def resolve_locus_tag_to_gene_id(locus_tag):
    try:
        handle = Entrez.esearch(db="gene", term=locus_tag, retmode="xml")
        results = Entrez.read(handle)
        handle.close()
        gene_ids = results['IdList']
        return gene_ids[0] if gene_ids else None
    except Exception as e:
        print(f"Failed to resolve locus tag {locus_tag} to Gene ID: {e}")
        return None



def link_gene_id_to_nuccore_accessions(gene_id):
    try:
        links = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id, linkname="gene_nuccore")
        link_results = Entrez.read(links)
        nuccore_ids = [link['Id'] for link in link_results[0]['LinkSetDb'][0]['Link']]
        return nuccore_ids
    except Exception as e:
        print(f"Failed to link Gene ID {gene_id} to nuccore accession numbers: {e}")
        return []



def fetch_cds_sequence(seq_id):
    if seq_id.startswith("LOC"):
        gene_id = resolve_locus_tag_to_gene_id(seq_id)
        if gene_id:
            nuccore_ids = link_gene_id_to_nuccore_accessions(gene_id)
            if nuccore_ids:
                seq_id = nuccore_ids[0]
            else:
                print(f"No linked nucleotide sequences found for {seq_id}")
                return None, None
        else:
            print(f"No Gene ID found for {seq_id}")
            return None, None

    # Existing sequence fetching logic for seq_id (now possibly resolved to a nuccore ID)
    try:
        with Entrez.efetch(db="nuccore", id=seq_id, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
            for feature in record.features:
                if feature.type == "CDS":
                    description = feature.qualifiers.get('product', ['No description'])[0]
                    header = f">{seq_id} | {description}"
                    sequence_data = str(feature.location.extract(record.seq))
                    return header, sequence_data
    except Exception as e:
        print(f"An error occurred while fetching data for {seq_id}: {e}")
    return None, None



def write_sequences(all_sequences, output_path):
    with open(output_path, 'w') as file:
        for header, seq in all_sequences:
            if header and seq:  # Ensure both header and sequence are not None
                file.write(f"{header}\n{seq}\n\n")


def main(input_file, output_file):
    ids = read_ids(input_file)
    all_sequences = []
    for seq_id in ids:
        header, sequence = fetch_cds_sequence(seq_id)
        all_sequences.append((header, sequence))
        if header and sequence:
            print(f"Sequence extraction for {seq_id} successful")
        else:
            print(f"No sequence extracted for {seq_id}")

    write_sequences(all_sequences, output_file)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file_path> <output_file_path>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
