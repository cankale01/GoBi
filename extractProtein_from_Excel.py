import pandas as pd
from Bio import Entrez, SeqIO
import sys
from io import StringIO

Entrez.email = "your_email@example.com"

def fetch_protein_info(loc_id):
    if loc_id:
        try:
            handle = Entrez.esearch(db="gene", term=loc_id, retmax=10)
            gene_record = Entrez.read(handle)
            handle.close()

            if not gene_record["IdList"]:
                return None, "Gene ID not found"

            gene_uid = gene_record["IdList"][0]
            protein_idss = []
            protein_sequences = []
            species = ''

            handle = Entrez.elink(dbfrom="gene", db="protein", id=gene_uid, linkname="gene_protein_refseq")
            link_record = Entrez.read(handle)
            handle.close()

            if link_record and link_record[0]["LinkSetDb"]:
                protein_ids = [link["Id"] for link in link_record[0]["LinkSetDb"][0]["Link"]]

                for prot_id in protein_ids:
                    handle = Entrez.efetch(db="protein", id=prot_id, rettype="gb", retmode="text")
                    protein_data = handle.read()
                    handle.close()

                    protein_data_io = StringIO(protein_data)
                    for protein_record in SeqIO.parse(protein_data_io, "gb"):
                        species = protein_record.annotations.get('organism', 'Unknown')
                        protein_id = protein_record.id
                        protein_idss.append(protein_id)
                        protein_sequence = str(protein_record.seq)
                        protein_sequences.append(protein_sequence)

            return species,loc_id, protein_idss, protein_sequences
        except Exception as e:
            return None, str(e)

    else:
        return f'{loc_id} can\'t be found \n'


def process_excel(file_path, output_path):
    printed_locs = []
    df = pd.read_excel(file_path, header=None, skiprows=1)

    with open(output_path, "w") as out_file:
        for col in df.columns[3:]:
            for cell in df[col]:
                if pd.notna(cell):
                    if str(cell).startswith(("LOC", "Dmel")):
                        cell_id = str(cell).split()[0]
                    elif 'NCBI:' in str(cell):
                        cell_id = str(cell).split('NCBI:')[-1].strip(')\'')
                    else:
                        continue  # Use continue to skip to the next iteration


                    try:
                        organism, loc_header, protein_id, aa_Seq = fetch_protein_info(cell_id)
                        if aa_Seq and loc_header not in printed_locs:
                            printed_locs.append(loc_header)
                            out_file.write(f"{organism}:\n")
                            out_file.write(f">{loc_header}:\n\n")
                            print(f"ID written: {loc_header}")  # Log when loc_header is written
                            for header, seq in zip(protein_id, aa_Seq):
                                out_file.write(f'{header}:\n')
                                for i in range(0, len(seq), 70):
                                    out_file.write(f"{seq[i:i + 70]}\n")
                                out_file.write('\n')
                    except ValueError as e:
                        print(f"Error processing cell_id {cell_id}: {e}")
                        # Handle the error as needed


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_excel_file_path> <output_file_path>")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_excel(input_file, output_file)
    print('done')
