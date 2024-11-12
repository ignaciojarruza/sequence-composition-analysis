# sequence-composition-analysis

This command-line (CLI) application provides tools to analyze DNA sequences, calculate nucleotide content, and generate various visualizations. Built
with Typer, this application is designed for ease of use and effective sequence data analysis.

## Features

- Calculate per-sequence nucleotide composition and content ratios
- Generate overall statistics across multiple sequences
- Visualize nucleotide composition, GC content distribution, sequence length distribution, and dinucleotide frequencies

## Prerequisites

- Python >= 3.8
- Packages:
- - Typer
- - Biopython
- - NumPy
- - Matplotlib
- - Seaborn
    You can quickly install all packages by running the following command:

```
pip3 install typer biopython numpy matplotlib seaborn
```

## Installation

1. Clone the repo:

```
git clone <repo>
```

2. Navigate to the project directory:

```
cd repo
```

## Usage

The CLI application supports multiple commands to perform various types of analyses. Run each command with `--help` to see available options and detailed information.

```
python3 seq-comp.py --help
```

## Commands

Below are descriptions and examples for each command available in the application:

1. `per-seq-stats`: calculates nucleotide composition and content ratios for each sequence in the provided file
2. `overall-stats`: provides general statistics across all sequences in the file, including total sequences, average length, and nucleotide counts
3. `visualize-nuc-bar`: generates a bar chart of nucleotide composition across all sequences
4. `visualize-gc-histo`: plots the GC content distribution across sequences
5. `visualize-seq-len-dist`: displays a histogram of sequence length distribution
6. `visualize-dinuc-heatmap`: creates a heatmap of dinucleotide frequencies in the sequences
