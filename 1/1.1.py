dna = input().strip()

#counting the number of occurrences of each nucleotide
A = dna.count('A')
C = dna.count('C')
G = dna.count('G')
T = dna.count('T')

print(A, C, G, T)