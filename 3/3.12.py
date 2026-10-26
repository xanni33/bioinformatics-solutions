#changing data formats
from Bio import Entrez, SeqIO

ids = input().split()
Entrez.email = "ania_gruca@interia.pl"

with Entrez.efetch(
    db = "nucleotide",
    id = ids,
    rettype = "fasta",
    retmode = "text"
) as handle:
    records = list(SeqIO.parse(handle, "fasta"))

shortest = min(records, key = lambda record: len(record.seq))

print(shortest.format("fasta"), end = "")
#sample dataset: FJ817486 JX069768 JX469983