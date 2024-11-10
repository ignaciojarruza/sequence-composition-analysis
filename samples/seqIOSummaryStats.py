import numpy as np
from collections import Counter
from Bio import SeqIO
sequences = list(SeqIO.parse("samples/human.protein.fasta", "fasta"))
nucleotide_counts = Counter()
for record in sequences:
    nucleotide_counts.update(record.seq)