#najwieksza zawartosc gc
from Bio import SeqIO

records = list(SeqIO.parse("rosalind_gc.txt", "fasta"))

best_id = ""
max_gc_content = 0

for record in records:
    sequence = str(record.seq)
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100

    if gc_content > max_gc_content:
        max_gc_content = gc_content
        best_id = record.id

print(best_id)
print(max_gc_content)