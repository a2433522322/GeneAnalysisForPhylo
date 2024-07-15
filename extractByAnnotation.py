from Bio import SeqIO
import os

# need change
result_file = os.getcwd()+"/annotation_statistics.txt"
output_dir = os.getcwd()+"/out"
os.makedirs(output_dir, exist_ok=True)

sequences = {}
with open(result_file, 'r') as f:
    next(f)  
    for line in f:
        category, annotation, count, files_count = line.strip().split('\t')
        count = int(count)
        files_count = int(files_count)
        if count == 253 and files_count == 253:#condition
            sequences[(category, annotation)] = []

input_dir = os.getcwd()
for file in os.listdir(input_dir):
    if file.endswith('.gb'):
        gb_file_name = os.path.splitext(file)[0]
        added_sequences = set()  
        with open(os.path.join(input_dir, file), 'r') as handle:
            for record in SeqIO.parse(handle, 'genbank'):
                for feature in record.features:
                    for key, value in feature.qualifiers.items():
                        for category_annotation in sequences:
                            category, annotation = category_annotation
                            if key == category and value[0] == annotation:
                                start = feature.location.start
                                end = feature.location.end
                                seq_fragment = record.seq[start:end]
                                location_str = str(feature.location)
                                if 'complement' in location_str:  
                                    seq_fragment = seq_fragment.reverse_complement()
                                    seq_info = (gb_file_name, str(seq_fragment)) 
                                else:
                                    seq_fragment = seq_fragment
                                    seq_info = (gb_file_name, str(seq_fragment))  
                                if seq_info not in added_sequences: 
                                    sequences[category_annotation].append(seq_info)
                                    added_sequences.add(seq_info) 


for category_annotation, seqs in sequences.items():
    category, annotation = category_annotation
    fasta_file = os.path.join(output_dir, f"{category}_{annotation}.fasta")
    with open(fasta_file, 'w') as f:
        for gb_file_name, seq in seqs:
            f.write(f">{gb_file_name}\n")
            f.write(f"{seq}\n")

print("All files have saved")
