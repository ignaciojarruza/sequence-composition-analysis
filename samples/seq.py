from Bio.Seq import Seq

# Creating sequence object
sequence = Seq("AGTACACTGGT")

# .find() example
sequence.find("ACT") # 5
sequence.find("TAG") # -1

# .count() example
sequence.count("A")
sequence.count("G")

# count_overlap() example
sequenceOverlap = Seq("AAAA")
sequenceOverlap.count_overlap("AA") # returns 3 instead of 2 (nonoverlap)

# complement() and reversecomplement()
sequence.complement()
sequence.reverse_complement()

# transcribe() and back_transcribe()
sequence.transcribe()
sequence.transcribe().back_transcribe()

