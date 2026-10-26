#genbank introduction
from Bio import Entrez

Entrez.email = "ania_gruca@interia.pl"

genus = input().strip()
date_start = input().strip()
date_end = input().strip()

query = (
    f'"{genus}"[Organism] AND '
    f'("{date_start}"[Publication Date] : '
    f'"{date_end}"[Publication Date])'
)
#print(query)

with Entrez.esearch(
    db = "nucleotide",
    term = query
) as handle:
    record = Entrez.read(handle)

print(record["Count"])

""" sample input:
Anthoxanthum
2003/7/25
2005/12/27 
"""