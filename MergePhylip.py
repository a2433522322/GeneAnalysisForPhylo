import os
from Bio import AlignIO
from collections import defaultdict

def read_phylip(file):
    alignment = AlignIO.read(file, "phylip")
    return alignment

def merge_alignments(directory):
    sequences = defaultdict(list)
    file_positions = []
    current_position = 1
    
    phylip_files = sorted([f for f in os.listdir(directory) if f.endswith('.fasta.phylip')])

    all_sample_ids = set()
    alignment_lengths = []

    for file in phylip_files:
        file_path = os.path.join(directory, file)
        alignment = read_phylip(file)
        alignment_lengths.append(alignment.get_alignment_length())
        
        for record in alignment:
            all_sample_ids.add(record.id)
    
    for sample_id in all_sample_ids:
        sequences[sample_id] = []

    for i, file in enumerate(phylip_files):
        file_path = os.path.join(directory, file)
        alignment = read_phylip(file)
        seq_length = alignment_lengths[i]
        
        present_sample_ids = set(record.id for record in alignment)
        missing_sample_ids = all_sample_ids - present_sample_ids

        for sample_id in missing_sample_ids:
            sequences[sample_id].append('-' * seq_length)
        
        for record in alignment:
            sequences[record.id].append(str(record.seq))
        
        start_position = current_position
        end_position = current_position + seq_length - 1
        identifier = file.replace('.fasta.phylip', '')
        file_positions.append(f"DNA, {identifier} = {start_position} - {end_position}")
        current_position = end_position + 1
    
    merged_records = []
    for sample_id, seqs in sequences.items():
        merged_seq = ''.join(seqs)
        merged_records.append((sample_id, merged_seq))
    
    return merged_records, file_positions

def write_phylip(output_file, merged_records):
    seq_length = len(merged_records[0][1])
    num_samples = len(merged_records)

    with open(output_file, 'w') as out_file:
        out_file.write(f"{num_samples} {seq_length}\n")
        for sample_id, seq in merged_records:
            out_file.write(f"{sample_id} {seq}\n")

def write_position_file(position_file, file_positions):
    with open(position_file, 'w') as pos_file:
        for line in file_positions:
            pos_file.write(f"{line}\n")


def main():
    directory = "/DATA/int8T/users/symsx/PLA_research/anno/out/mafft"  # the phylip files after mafft
    output_phylip_file = "merge.phylip"
    output_position_file = "gene_partition.txt"

    merged_records, file_positions = merge_alignments(directory)


    write_phylip(output_phylip_file, merged_records)

    write_position_file(output_position_file, file_positions)

    print(f"successful! Result has been saved in {output_phylip_file}")
    print(f"Partion in {output_position_file}")

if __name__ == "__main__":
    main()
