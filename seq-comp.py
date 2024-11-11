from Bio import SeqIO
from collections import Counter
import numpy as np
import sys

def sequence_statistics_overall(file_path: str, file_format: str):
    try:
        sequences = list(SeqIO.parse(file_path, file_format))

        # Basic statistics
        total_sequences = len(sequences)
        sequence_lengths = [len(record.seq) for record in sequences]
        gc_contents = [(record.seq.count("G") + record.seq.count("C")) / len(record.seq) * 100 for record in sequences]
        # calculate nucleotide composition across all sequences
        nucleotide_counts = Counter()
        for record in sequences:
            nucleotide_counts.update(record.seq)

        # General Stats Print-out 
        print(f"Total Sequences: {total_sequences}")
        print(f"Average Sequence Length: {np.mean(sequence_lengths):.2f}")
        print(f"Minimum Sequence Length: {np.min(sequence_lengths)}")
        print(f"Maximum Sequence Length: {np.max(sequence_lengths)}")
        print(f"Average GC Content: {np.mean(gc_contents):.2f}%")
        print("Nucleotide Count:")
        for nucleotide, count in nucleotide_counts.items():
            print(f"{nucleotide}: {count}")
    except FileNotFoundError:
        print("File not found.")
        sys.exit(1)
    except ValueError:
        print("Format not supported. Please execute -h for supported file formats.")

def sequence_statistics_iteratively(file_path: str, format: str):
    return ""


sequence_statistics_overall(sys.argv[1], sys.argv[2])
sequence_statistics_iteratively()