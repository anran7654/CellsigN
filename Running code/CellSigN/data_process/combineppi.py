import os
from collections import defaultdict


def list_files(folder_path):
    # List all files in the directory
    result = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            file_path = os.path.join(root, file)
            result.append(file_path)
    return result


def count_overlaps(dicts):
    # Use a defaultdict to store the count of each key-value pair
    overlap_count = defaultdict(int)
    # Iterate over each dictionary
    for d in dicts:
        for key, values in d.items():
            # Iterate over each value in the set
            for value in values:
                # Create a tuple (key, value) to use as a unique identifier in the defaultdict
                overlap_count[(key, value)] += 1
    return overlap_count


folder_path = 'D:/Personal issues/CSNet/tool/interaction data/human/PPI/each/'
files = list_files(folder_path)

dicts = []
for file_path in files:
    print(file_path)
    dictionary = {}
    with open(file_path, 'r', encoding='utf-8') as file:
        first_line = next(file)
        for line in file:
            elements = line.strip().split('\t')  # Split each line into elements based on delimiter \t
            gene1, gene2 = elements[0], elements[1]
            if gene1 <= gene2:
                if gene1 not in dictionary:
                    dictionary[gene1] = set()
                dictionary[gene1].add(gene2)
            else:
                if gene2 not in dictionary:
                    dictionary[gene2] = set()
                dictionary[gene2].add(gene1)
    dicts.append(dictionary)


print(len(dicts))
overlap_counts = count_overlaps(dicts)

with open(f"D:/Personal issues/CSNet/tool/interaction data/human/PPI/combined_ppi.txt", "w") as f:
    f.write(f"Protein1\tProtein2\tScore\n")
    for key_value, count in overlap_counts.items():
        key, value = key_value  # Unpack the tuple
        f.write(f"{key}\t{value}\t{count}\n")

