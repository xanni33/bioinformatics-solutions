#overlaped graphs dla sekwencji dna
#ostatnie 3 nukleotydy z pierwszej sekwencji pokrywaja sie z pierwszymi 3 nukleotydami drugiej

from Bio import SeqIO

records = list(SeqIO.parse("rosalind_grph.txt", "fasta"))

#string s polaczony ze stringiem t krawedzia skierowana jesli istnieje k-dlugosci suffix z s ktory koresponduje z k-dl. prefixem t
for s in records:
    for t in records:
        if s.id != t.id and str(s.seq)[-3:] == str(t.seq)[:3]:
            print(s.id, t.id)