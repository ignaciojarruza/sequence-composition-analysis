import numpy as np
from collections import Counter
from Bio import SeqIO

sequences = list(SeqIO.parse("samples/human.protein.fasta", "fasta"))

# Basic statistics
total_sequences = len(sequences)
sequence_lengths = [len(record.seq) for record in sequences]
gc_contents = [(record.seq.count("G") + record.seq.count("C")) / len(record.seq) * 100 for record in sequences]

# general stats
print(f"Total Sequences: {total_sequences}")
print(f"Average Sequence Length: {np.mean(sequence_lengths):.2f}")
print(f"Minimum Sequence Length: {np.min(sequence_lengths)}")
print(f"Maximum Sequence Length: {np.max(sequence_lengths)}")
print(f"Average GC Content: {np.mean(gc_contents):.2f}%")

# calculate nucleotide composition across all sequences
nucleotide_counts = Counter()
for record in sequences:
    nucleotide_counts.update(record.seq)

for nucleotide, count in nucleotide_counts.items():
    print(f"{nucleotide}: {count}")