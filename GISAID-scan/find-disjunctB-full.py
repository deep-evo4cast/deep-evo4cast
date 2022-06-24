import re2
from Bio import SeqIO
import sys

records = SeqIO.parse("spike-1260.fasta2", "fasta")

with open(f'full-wc.pats') as f:
    global dpat
    rpats = f.readlines()
    pats = [rpat.strip() for rpat in rpats]
    disj = '|'.join(pats)
    dpat = re2.compile(disj,max_mem=900000000)

lhits = []
for r in records:
    m = dpat.search(r.seq.__str__()[380:515])
    if m:
        lhits.append(r)

SeqIO.write(lhits, f'full-wc.hits', "fasta-2line")
