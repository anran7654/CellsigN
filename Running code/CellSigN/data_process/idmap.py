from Bio import Entrez


def check_species(gene_symbol):  # whether a gene symbol belongs to homo species
    Entrez.email = "free1234hm@163.com"  # Set your email for NCBI API
    handle = Entrez.esearch(db="gene", term=gene_symbol + "[Gene Name] AND Homo sapiens[Organism]")
    record = Entrez.read(handle)
    return bool(record["IdList"])  # If the list is not empty, gene symbol is associated with Homo sapiens


def uppercase(string):  # whether a string is uppercase
    for char in string:
        if char.isalpha() and not char.isupper():
            return False
    return True


filepath = ('D:/Personal issues/CSNet/tool/interaction '
            'data/mouse/Regulator-target/miRNA-target/miRDB/idmap_mus.txt')
dictionary = {}
with open(filepath, 'r', encoding='utf-8') as file:
    # next(file)  # Read and discard the first line
    for line in file:
        elements = line.strip().split('\t')  # Split each line into elements based on delimiter \t
        if len(elements) > 1:
            gene1, gene2 = elements[0], elements[1]
            if gene1 not in dictionary:
                dictionary[gene1] = set()  # Create a new set if it doesn't exist
            dictionary[gene1].add(gene2)


filepath = ('D:/Personal issues/CSNet/tool/interaction '
            'data/mouse/Regulator-target/miRNA-target/miRDB/miRDB_mus.txt')
result = []
with open(filepath, 'r', encoding='utf-8') as file:
    line1 = next(file)  # Read and discard the first line
    for line in file:
        elements = line.strip().split('\t')  # Split each line into elements based on delimiter \t
        gene1, gene2 = elements[0], elements[1]
        if gene2 in dictionary:
            set2 = dictionary[gene2]
            for g2 in set2:
                result.append(gene1 + '\t' + g2)

with open(f'D:/Personal issues/CSNet/tool/interaction '
          'data/mouse/Regulator-target/miRNA-target/miRDB/miRDB_mus2.txt', "w") as f:
    f.write(f"miRNA\tTarget\n")
    for row in result:
        f.write(f"{row}\n")
