from Bio import SeqIO

with open("samples/human.protein.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        sequence = str(record.seq)

        g_count = sequence.count("G")
        t_count = sequence.count("T")
        a_count = sequence.count("A")
        c_count = sequence.count("C")

        gc_content = ((g_count + c_count) / len(sequence)) * 100
        at_content = ((a_count + t_count) / len(sequence)) * 100

        at_ratio = a_count / t_count if t_count != 0 else 0
        gc_ratio = g_count / c_count if c_count != 0 else 0

        purine_pyrimidine_ratio = (a_count + g_count) / (c_count + t_count) if (c_count + t_count) != 0 else 0

        gc_ac_ratio = (g_count + c_count) / (a_count + t_count) if (a_count + t_count) != 0 else 0

        print(f"Sequence ID: {record.id}")
        print(f"A: {a_count}, T {t_count}, C: {c_count}, G:{g_count}")
        print(f"GC Content: {gc_content:.2f}%")
        print(f"AT Content: {at_content:.2f}%")
        print(f"A/T Ratio: {at_ratio}")
        print(f"G/C Ratio: {gc_ratio}")
        print(f"Purine/Pyrimidine Ratio: {purine_pyrimidine_ratio}")
        print(f"GC/AT Ratio: {gc_ac_ratio}")
