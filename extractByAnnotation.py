from Bio import SeqIO
import os


result_file = "/DATA/annotation_statistics.txt"#statistic information
output_dir = "/DATA/int8T/users/symsx/coe_research/ONLY_coe/out"#output folder
os.makedirs(output_dir, exist_ok=True)


sequences = {}
with open(result_file, 'r') as f:
    next(f)
    for line in f:
        category, annotation, count, files_count = line.strip().split('\t')
        count = int(count)
        files_count = int(files_count)
        if files_count == 229:
            sequences[(category, annotation)] = []


input_dir = '/DATA/ONLY_coe'#input folder which has gb files
for file in os.listdir(input_dir):
    if file.endswith('.gb'):
        gb_file_name = os.path.splitext(file)[0]
        with open(os.path.join(input_dir, file), 'r') as handle:
            for record in SeqIO.parse(handle, 'genbank'):
                dic = {}
                for feature in record.features:
                    for key, value in feature.qualifiers.items():
                        for category_annotation in sequences:
                            category, annotation = category_annotation
                            if key == category and value[0] == annotation:
                                start = feature.location.start
                                end = feature.location.end
                                seq_fragment = record.seq[start:end]
                                location_str = str(feature.location)
                                if annotation in dic:
                                    pass
                                else:
                                    if feature.strand == -1:  
                                        seq_fragment = seq_fragment.reverse_complement()
                                    seq_info = (gb_file_name, str(seq_fragment))  
                                    sequences[category_annotation].append(seq_info)
                                    dic[annotation] = 0
                                    

for category_annotation, seqs in sequences.items():
    category, annotation = category_annotation
    fasta_file = os.path.join(output_dir, f"{category}_{annotation}.fasta")
    with open(fasta_file, 'w') as f:
        for gb_file_name, seq in seqs:
            f.write(f">{gb_file_name}\n")
            f.write(f"{seq}\n")

print("All results saved.")
