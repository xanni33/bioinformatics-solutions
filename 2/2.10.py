from Bio import SeqIO

record = next(SeqIO.parse("rosalind_orf.txt", "fasta"))
dna = record.seq

strands = [dna, dna.reverse_complement()]
proteins = set()

for strand in strands:
    for frame in range(3):
        sequence = strand[frame:]
        sequence = sequence[:len(sequence) - len(sequence) % 3]
        amino_acids = str(sequence.translate())

        start = 0
        while True:
            start = amino_acids.find("M", start)

            if start == -1:
                break

            stop = amino_acids.find("*", start)

            if stop != -1:
                proteins.add(amino_acids[start:stop])

            start += 1

for protein in proteins:
    print(protein)