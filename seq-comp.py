import typer
import sys
from Bio import SeqIO
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

app = typer.Typer()

@app.command()
def per_seq_stats(file_path: str, file_format: str):
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

@app.command()
def overall_stats(file_path: str, file_format: str):
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

@app.command()
def visualize_nuc_bar(file_path: str, file_format: str):
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

@app.command()
def visualize_gc_histo(file_path: str, file_format: str):
    try:
        sequences = list(SeqIO.parse(file_path, file_format))
        gc_contents = [(record.seq.count("G") + record.seq.count("C")) / len(record.seq) * 100 for record in sequences]
        sns.histplot(gc_contents, bins=20, kde=True, color='skyblue')
        plt.title("GC Content Distribution")
        plt.xlabel("GC Content (%)")
        plt.ylabel("Frequency")
        plt.show()
    except FileNotFoundError:
        print("File not found.")
        sys.exit(1)
    except ValueError:
        print("Format not supported. Please execute -h for supported file formats.")

@app.command()
def visualize_seq_len_dist(file_path: str, file_format: str):
    try:
        sequences = list(SeqIO.parse(file_path, file_format))
        sequence_lengths = [len(record.seq) for record in sequences]
        sns.histplot(sequence_lengths, bins=30, kde=True, color='salmon')
        plt.title("Sequence Length Distribution")
        plt.xlabel("Sequence Length (bp)")
        plt.ylabel("Frequency")
        plt.show()
    except FileNotFoundError:
        print("File not found.")
        sys.exit(1)
    except ValueError:
        print("Format not supported. Please execute -h for supported file formats.")

@app.command()
def visualize_dinuc_heatmap(file_path: str, file_format: str):
    try:
        sequences = list(SeqIO.parse(file_path, file_format))
        dinucleotides = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
        dinuc_counts = Counter({dinuc: 0 for dinuc in dinucleotides})
        for record in sequences:
            seq = str(record.seq).upper()  # Convert to uppercase to avoid case issues
            for i in range(len(seq) - 1):
                dinuc = seq[i:i+2]
                if dinuc in dinuc_counts:
                    dinuc_counts[dinuc] += 1

        # Calculate total dinucleotides counted
        total_dinucleotides = sum(dinuc_counts.values())

        # Convert counts to frequencies
        dinucleotide_frequencies = {dinuc: count / total_dinucleotides for dinuc, count in dinuc_counts.items()}

        # Reshape frequencies into a 4x4 matrix for plotting
        freq_matrix = np.array([dinucleotide_frequencies[dinuc] for dinuc in dinucleotides]).reshape(4, 4)

        # Plot heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(freq_matrix, annot=True, xticklabels=['A', 'T', 'G', 'C'], yticklabels=['A', 'T', 'G', 'C'])
        plt.title("Dinucleotide Frequency Heatmap")
        plt.xlabel("Second Nucleotide")
        plt.ylabel("First Nucleotide")
        plt.show()
    except FileNotFoundError:
        print("File not found.")
        sys.exit(1)
    except ValueError:
        print("Format not supported. Please execute -h for supported file formats.")

if __name__ == "__main__":
    app()