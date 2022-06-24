from Bio import SeqIO

records = SeqIO.parse("spikeprot0125.fasta", "fasta")

l = []
for r in records:
    if len(r.seq) >= 1260:
        l.append(r)

SeqIO.write(l, "spike-1260.fasta2", "fasta-2line")
