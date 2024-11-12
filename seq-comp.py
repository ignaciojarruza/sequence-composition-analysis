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

def sequence_statistics_iteratively(file_path: str, file_format: str):
    try:
        with open(file_path) as handle:
            for record in SeqIO.parse(handle, file_format):
                sequence = str(record.seq)
                # Sequence Count 
                g_count = sequence.count("G")
                t_count = sequence.count("T")
                a_count = sequence.count("A")
                c_count = sequence.count("C")
                # G/C and A/T content
                gc_content = ((g_count + c_count) / len(sequence)) * 100
                at_content = ((a_count + t_count) / len(sequence)) * 100
                # Interest Ratios
                at_ratio = a_count / t_count if t_count != 0 else 0
                gc_ratio = g_count / c_count if c_count != 0 else 0
                purine_pyrimidine_ratio = (a_count + g_count) / (c_count + t_count) if (c_count + t_count) != 0 else 0
                gc_ac_ratio = (g_count + c_count) / (a_count + t_count) if (a_count + t_count) != 0 else 0

                # Print-Out
                print(f"Sequence ID: {record.id}")
                print(f"A: {a_count}, T {t_count}, C: {c_count}, G:{g_count}")
                print(f"GC Content: {gc_content:.2f}%")
                print(f"AT Content: {at_content:.2f}%")
                print(f"A/T Ratio: {at_ratio}")
                print(f"G/C Ratio: {gc_ratio}")
                print(f"Purine/Pyrimidine Ratio: {purine_pyrimidine_ratio}")
                print(f"GC/AT Ratio: {gc_ac_ratio}")
                print("--------------------------------------------------------")
    except FileNotFoundError:
        print("File not found.")
        sys.exit(1)
    except ValueError:
        print("Format not supported. Please execute -h for supported file formats.")

if len(sys.argv) < 4:
    print("Invalid command length, usage: python3 seq-comp.py command file_path file_format")
    sys.exit(1)
if sys.argv[1] == "iter":
    sequence_statistics_iteratively(sys.argv[2], sys.argv[3])
elif sys.argv[1] == "whole":
    sequence_statistics_overall(sys.argv[2], sys.argv[3])
else:
    print("Invalid command, usage: python3 seq-comp.py command file_path file_format")
    print("Supported seq-comp commands are: [iter, whole]")
    sys.exit(1)
