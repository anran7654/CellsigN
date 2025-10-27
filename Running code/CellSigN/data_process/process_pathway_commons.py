import os
import numpy as np
from scipy.sparse import csr_matrix


lr_file = 'D:/Personal issues/CSNet/interaction data/human/ligand-receptor.txt'
tf_file = 'D:/Personal issues/CSNet/interaction data/human/TF-target.txt'
PathwayCommons = 'D:/Personal issues/CSNet/interaction data/human/PathwayCommons/pc-hgnc.txt'

re2lg = {}
if os.path.isfile(lr_file):  # Check if it's a file (not a directory)
    with open(lr_file, 'r', encoding='utf-8') as file:
        next(file)
        for line in file:
            elements = line.strip().split('\t')  # Split each line into elements based on delimiter \t
            if len(elements) >= 2:
                lg = elements[0]
                re = elements[1]
                re2lg[re] = lg
print(len(re2lg))

tf2tg = {}
if os.path.isfile(tf_file):  # Check if it's a file (not a directory)
    with open(tf_file, 'r', encoding='utf-8') as file:
        next(file)
        for line in file:
            elements = line.strip().split('\t')  # Split each line into elements based on delimiter \t
            if len(elements) >= 2:
                tf = elements[0]
                tg = elements[1]
                tf2tg[tf] = tg
print(len(tf2tg))

result = []
dict = {}
if os.path.isfile(PathwayCommons):  # Check if it's a file (not a directory)
    with (open(PathwayCommons, 'r', encoding='utf-8') as file):
        for line in file:
            elements = line.strip().split('\t')
            if len(elements) >= 4:
                re = elements[0]
                type = elements[1]
                tf = elements[2]
                source = elements[3]
                if type == 'controls-phosphorylation-of' or type == 'controls-state-change-of' or type == 'interacts-with' or type == 'in-complex-with' or type == 'controls-transport-of'  or type == 'catalysis-precedes'  or type == 'controls-expression-of':
                    if type not in dict:
                        dict[type] = set()
                    dict[type].add(source)
                    result.append(f'{re}\t{tf}\t{type}\t{source}')

with open("D:/dictionary.txt", "w") as file:
    # Iterate over the dictionary items
    for key, value in dict.items():
        # Convert the set to a sorted list to ensure consistent order
        sorted_values = sorted(value)
        # Write the key and the set elements to the file
        file.write(f"{key}: {', '.join(sorted_values)}\n")

with open("D:/Pathway_Commons.txt", 'w') as file:
    file.write('Participant_1\tParticipant_2\tType\n')
    for element in result:
        file.write(f"{element}\n")