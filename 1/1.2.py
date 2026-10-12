#odwrotna sekwencja komplementarna dna
dna = input().strip()

#pierwsze zamiana na komplementarne nukleotydy potem odwrócenie sekwencji
complement = dna.translate(str.maketrans("ACGT", "TGCA"))
reverse_complement = complement[::-1]

print(reverse_complement)
