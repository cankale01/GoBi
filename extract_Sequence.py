def extract_sequences_by_ids(file, identify, out):
    with open(file, 'r') as fasta_file, open(out, 'w') as output_file:
        for genomic_id, start, end in identify:
            fasta_file.seek(0)
            sequence = ''
            match_found = False
            for line in fasta_file:
                if line.startswith('>'):
                    line_header = line.strip().split()[0][1:]
                    if line_header == genomic_id:
                        header = line.strip()
                        match_found = True
                    else:
                        match_found = False
                elif match_found:
                    sequence += line.strip()
                    if len(sequence) >= end:
                        break
            if match_found:
                extracted_sequence = sequence[start-1:end]
                output_file.write(f"{header}\n{extracted_sequence}\n\n")
    print(f"Regions successfully extracted to: {output_file_path}")


file_path = '/Users/cankale/Desktop/sequence.fasta'
identifiers = [('NC_037648.1',2230673,2342854)]
output_file_path = '/Users/cankale/Desktop/gen_dataset/Regions/region_mellifera.fna'

extract_sequences_by_ids(file_path, identifiers, output_file_path)