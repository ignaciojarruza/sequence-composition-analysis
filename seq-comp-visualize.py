import matplotlib.pyplot as plt
from collections import Counter
import sys
from Bio import SeqIO

def nuc_bar(file_path: str, file_format: str):
    try:
        sequences = list(SeqIO.parse(file_path, file_format))
        nucleotide_counts = Counter()
        for record in sequences:
            nucleotide_counts.update(record.seq)

        nucleotides = ['A', 'T', 'C', 'G']
        counts = [nucleotide_counts.get('A'), nucleotide_counts.get('T'), nucleotide_counts.get('G'), nucleotide_counts.get('C')]
        plt.bar(nucleotides, counts, color=['blue', 'green', 'orange', 'red'])
        plt.title("Nucleotide Composition")
        plt.xlabel("Nucleotide")
        plt.ylabel("Count")
        for i, count in enumerate(counts):
            plt.text(i, count, f'{count:,}', ha='center', va='bottom')
        plt.show()
    except FileNotFoundError:
        print("File not found.")
        sys.exit(1)
    except ValueError:
        print("Format not supported. Please execute -h for supported file formats.")

nuc_bar(sys.argv[1], sys.argv[2])
