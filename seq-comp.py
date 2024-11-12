import typer
import sys
from Bio import SeqIO

app = typer.Typer()

@app.command()
def overall_stats(file_path: str, file_format: str):
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
def per_seq_stats():
    print("")

if __name__ == "__main__":
    app()