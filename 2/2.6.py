from Bio import SeqIO

sequences = [str(record.seq) for record in SeqIO.parse("rosalind_lcsm.txt", "fasta")]

shortest = min(sequences, key=len)

for length in range(len(shortest), 0, -1):
    for start in range(len(shortest) - length + 1):
        substring = shortest[start:start + length]

        if all(substring in seq for seq in sequences):
            print(substring)
            exit()