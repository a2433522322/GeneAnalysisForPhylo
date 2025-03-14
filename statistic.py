from collections import defaultdict
from Bio import SeqIO
import glob
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s', handlers=[
    logging.FileHandler("process_log.txt"),
    logging.StreamHandler()
])

def process_gb_file(file_path):
    annotations_in_file = defaultdict(lambda: defaultdict(int))
    invalid_annotations = [] 
    

    try:
        with open(file_path, 'r') as handle:
            for record in SeqIO.parse(handle, 'genbank'):
                for feature in record.features:
                    if feature.type == 'gene':
                        try:
                            
                            if feature.location.end < feature.location.start:
                                invalid_annotations.append((file_path, feature))
                                continue
                            
                            annotation_type = feature.type
                            for key, value in feature.qualifiers.items():
                                for item in value:
                                    annotations_in_file[annotation_type][item] += 1
                        except ValueError as e:
                            invalid_annotations.append((file_path, feature, str(e)))
    except Exception as e:
        logging.info(f"Error processing file {file_path}: {e}")
    
    return annotations_in_file, invalid_annotations

annotations_count = defaultdict(lambda: defaultdict(int))  
annotations_files = defaultdict(lambda: defaultdict(int))  
invalid_annotations_list = []  

input_dir = '/DATA/'# the folder which has gb files

gb_files = glob.glob(os.path.join(input_dir, '*.gb'))  


for file in gb_files:
    file_annotations, invalid_annotations = process_gb_file(file)
    invalid_annotations_list.extend(invalid_annotations)

    for annotation_type in file_annotations:
        for annotation, count in file_annotations[annotation_type].items():
            annotations_count[annotation_type][annotation] += count
            annotations_files[annotation_type][annotation] += 1

output_file = 'annotation_statistics.txt'
invalid_output_file = 'invalid_annotations.txt'

with open(output_file, 'w') as f:
    f.write("type\tannotation\tTotal counts\tFile counts\n")
    for annotation_type in annotations_count:
        for annotation in annotations_count[annotation_type]:
            f.write(f"{annotation_type}\t{annotation}\t{annotations_count[annotation_type][annotation]}\t{annotations_files[annotation_type][annotation]}\n")



logging.info(f"All results saved in {output_file}")
