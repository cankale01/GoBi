

blast_result = '/Users/cankale/Desktop/gen_dataset/EXON_blasts/blast_exon_tribolium_castaneum'
output_file = '/Users/cankale/Desktop/gen_dataset/EXON_best_hits/exon_best_castaneum'
Eval = None
threshold = float(0.05)
queries = []

with open(blast_result,'r',encoding='utf-8') as file, open(output_file,'w',encoding='utf-8') as out:
    lines = file.readlines()
    for line in lines:
        Eval = float(line.split()[10])
        if Eval < threshold:
            queries.append(line)
    out.write(''.join(queries))
    print(f'Successfully written to {output_file}')