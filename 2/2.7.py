from Bio import SeqIO

records = list(SeqIO.parse("rosalind_tran.txt", "fasta"))

s1 = str(records[0].seq)
s2 = str(records[1].seq)

transitions = 0
transversions = 0

purines = {'A', 'G'}
pyrimidines = {'C', 'T'}

for a, b in zip(s1, s2):
    if a == b:
        continue

    if (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
        transitions += 1
    else:
        transversions += 1

#transition/transversion ratio
print (transitions / transversions)