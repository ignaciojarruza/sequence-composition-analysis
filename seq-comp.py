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
    '''
    Calculate per-sequence nucleotide statistics for a given file.

    This command reads a DNA sequence file in the specified format and calculates nucleotide counts (A, T, C, G), GC and AT content 
    percentages, and various nucleotide ratios.
    Each sequence's statistics are printed separately, including:
    - GC and AT content percentages
    - A/T and G/C ratios
    - Purine/Pyrimidine ratio
    - GC/AT ratio

    Parameters:
    file_path : str
        The path to the DNA sequence file (e.g., 'human.protein.fasta').
        
    file_format : str
        The format of the sequence file (e.g., 'fasta', 'fastq'). Refer to Biopython's SeqIO module for supported formats.

    Raises:
    FileNotFoundError:
        If the specified file path does not exist.
    
    ValueError:
        If the specified file format is unsupported.
    '''
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
                typer.echo(f"Sequence ID: {record.id}")
                typer.echo(f"A: {a_count}, T {t_count}, C: {c_count}, G:{g_count}")
                typer.echo(f"GC Content: {gc_content:.2f}%")
                typer.echo(f"AT Content: {at_content:.2f}%")
                typer.echo(f"A/T Ratio: {at_ratio}")
                typer.echo(f"G/C Ratio: {gc_ratio}")
                typer.echo(f"Purine/Pyrimidine Ratio: {purine_pyrimidine_ratio}")
                typer.echo(f"GC/AT Ratio: {gc_ac_ratio}")
                typer.echo("--------------------------------------------------------")
    except FileNotFoundError:
        typer.echo("File not found.", err=True)
        raise typer.Exit(code=1)
    except ValueError:
        typer.echo("Format not supported. Please execute -h for supported file formats.", err=True)

@app.command()
def overall_stats(file_path: str, file_format: str):
    '''
    Calculate and display overall nucleotide statistics for a given file.

    This command provides an overview of the DNA sequences in the specified file, including:
    - Total number of sequences
    - Average, minimum, and maximum sequence lengths
    - Average GC content percentage across all sequences
    - Nucleotide counts across all sequences (A, T, C, G, etc.)

    Parameters:
    file_path : str
        The path to the DNA sequence file (e.g., 'human.protein.fasta').
        
    file_format : str
        The format of the sequence file (e.g., 'fasta', 'fastq'). Refer to Biopython's SeqIO module for supported formats.

    Raises:
    FileNotFoundError:
        If the specified file path does not exist.
    
    ValueError:
        If the specified file format is unsupported.
    '''
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
        typer.echo(f"Total Sequences: {total_sequences}")
        typer.echo(f"Average Sequence Length: {np.mean(sequence_lengths):.2f}")
        typer.echo(f"Minimum Sequence Length: {np.min(sequence_lengths)}")
        typer.echo(f"Maximum Sequence Length: {np.max(sequence_lengths)}")
        typer.echo(f"Average GC Content: {np.mean(gc_contents):.2f}%")
        typer.echo("Nucleotide Count:")
        for nucleotide, count in nucleotide_counts.items():
            typer.echo(f"{nucleotide}: {count}")
    except FileNotFoundError:
        typer.echo("File not found.", err=True)
        raise typer.Exit(code=1)
    except ValueError:
        typer.echo("Format not supported. Please execute -h for supported file formats.", err=True)

@app.command()
def visualize_nuc_bar(file_path: str, file_format: str):
    '''
    Parameters:
    file_path : str
        The path to the DNA sequence file (e.g., 'human.protein.fasta').
        
    file_format : str
        The format of the sequence file (e.g., 'fasta', 'fastq'). Refer to Biopython's SeqIO module for supported formats.

    Raises:
    FileNotFoundError:
        If the specified file path does not exist.
    
    ValueError:
        If the specified file format is unsupported.
    '''
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
        typer.echo("File not found.", err=True)
        raise typer.Exit(code=1)
    except ValueError:
        typer.echo("Format not supported. Please execute -h for supported file formats.", err=True)

@app.command()
def visualize_gc_histo(file_path: str, file_format: str):
    '''
    Generate and display a bar chart of nucleotide composition in the DNA sequences.

    This command reads DNA sequences from the specified file, counts occurrences of each nucleotide (A, T, C, G), and 
    visualize the counts in a bar chart. Each nucleotide is displayed with a distinct color.

    Parameters:
    file_path : str
        The path to the DNA sequence file (e.g., 'human.protein.fasta').
        
    file_format : str
        The format of the sequence file (e.g., 'fasta', 'fastq'). Refer to Biopython's SeqIO module for supported formats.

    Raises:
    FileNotFoundError:
        If the specified file path does not exist.
    
    ValueError:
        If the specified file format is unsupported.
    '''
    try:
        sequences = list(SeqIO.parse(file_path, file_format))
        gc_contents = [(record.seq.count("G") + record.seq.count("C")) / len(record.seq) * 100 for record in sequences]
        sns.histplot(gc_contents, bins=20, kde=True, color='skyblue')
        plt.title("GC Content Distribution")
        plt.xlabel("GC Content (%)")
        plt.ylabel("Frequency")
        plt.show()
    except FileNotFoundError:
        typer.echo("File not found.", err=True)
        raise typer.Exit(code=1)
    except ValueError:
        typer.echo("Format not supported. Please execute -h for supported file formats.", err=True)

@app.command()
def visualize_seq_len_dist(file_path: str, file_format: str):
    '''
    Generate and display a histogram of sequence lengths in the DNA sequences.

    This command reads DNA sequences from the specified file, calculates their lengths,
    and visualizes the distribution of sequence lengths in a histogram. A kernel density
    estimate (KDE) curve is also displayed to show the distribution.

    Parameters:
    file_path : str
        The path to the DNA sequence file (e.g., 'human.protein.fasta').
        
    file_format : str
        The format of the sequence file (e.g., 'fasta', 'fastq'). Refer to Biopython's SeqIO module for supported formats.

    Raises:
    FileNotFoundError:
        If the specified file path does not exist.
    
    ValueError:
        If the specified file format is unsupported.'''
    try:
        sequences = list(SeqIO.parse(file_path, file_format))
        sequence_lengths = [len(record.seq) for record in sequences]
        sns.histplot(sequence_lengths, bins=30, kde=True, color='salmon')
        plt.title("Sequence Length Distribution")
        plt.xlabel("Sequence Length (bp)")
        plt.ylabel("Frequency")
        plt.show()
    except FileNotFoundError:
        typer.echo("File not found.", err=True)
        raise typer.Exit(code=1)
    except ValueError:
        typer.echo("Format not supported. Please execute -h for supported file formats.", err=True)

@app.command()
def visualize_dinuc_heatmap(file_path: str, file_format: str):
    '''
    Generate and display a heatmap of dinucleotide frequencies in DNA sequences.

    This command reads DNA sequences from the specified file, calculates the frequency
    of each dinucleotide pair (e.g., AA, AT, GC, etc.), and visualizes the frequencies
    as a heatmap. The heatmap shows the relationships between the first and second nucleotides
    in each dinucleotide pair.
    
    Parameters:
    file_path : str
        The path to the DNA sequence file (e.g., 'human.protein.fasta').
        
    file_format : str
        The format of the sequence file (e.g., 'fasta', 'fastq'). Refer to Biopython's SeqIO module for supported formats.

    Raises:
    FileNotFoundError:
        If the specified file path does not exist.
    
    ValueError:
        If the specified file format is unsupported.
    '''
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

        total_dinucleotides = sum(dinuc_counts.values())
        dinucleotide_frequencies = {dinuc: count / total_dinucleotides for dinuc, count in dinuc_counts.items()}
        freq_matrix = np.array([dinucleotide_frequencies[dinuc] for dinuc in dinucleotides]).reshape(4, 4)
        plt.figure(figsize=(8, 6))
        sns.heatmap(freq_matrix, annot=True, xticklabels=['A', 'T', 'G', 'C'], yticklabels=['A', 'T', 'G', 'C'])
        plt.title("Dinucleotide Frequency Heatmap")
        plt.xlabel("Second Nucleotide")
        plt.ylabel("First Nucleotide")
        plt.show()
    except FileNotFoundError:
        typer.echo("File not found.", err=True)
        raise typer.Exit(code=1)
    except ValueError:
        typer.echo("Format not supported. Please execute -h for supported file formats.", err=True)

if __name__ == "__main__":
    app()