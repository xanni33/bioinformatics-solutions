#pipeline dla wyszukania i pobrania informacji dotyczących mutacji w ludzkim genie BRCA1
from Bio import Entrez, SeqIO
import xml.etree.ElementTree as ET

Entrez.email = "ania_gruca@interia.pl"

def read_entrez(function, **kwargs):
    with function(**kwargs) as handle:
        return Entrez.read(handle)


def get_linked_ids(dbfrom, db, gene_id, term=None):
    args = {
        "dbfrom": dbfrom,
        "db": db,
        "id": gene_id
    }
    if term is not None:
        args["term"] = term

    result = read_entrez(Entrez.elink, **args)
    ids = []

    for link_set in result:
        for link_db in link_set.get("LinkSetDb", []):
            ids.extend(link["Id"] for link in link_db["Link"])

    return list(dict.fromkeys(ids))


def get_text(record, key, default="brak danych"):
    value = record.get(key, default)
    return str(value) if value else default


# krok 1 - wyszukiwanie genu brca1 u czlowieka
print("\n1. GENE: BRCA1")

result = read_entrez(
    Entrez.esearch,
    db="gene",
    term="BRCA1[Gene Name] AND Homo sapiens[Organism]"
)

gene_ids = result["IdList"]

if not gene_ids:
    raise RuntimeError("Nie znaleziono genu BRCA1.")

summaries = read_entrez(
    Entrez.esummary,
    db="gene",
    id=gene_ids
)

gene_id = None

for gene in summaries["DocumentSummarySet"]["DocumentSummary"]:
    print("\nName:", get_text(gene, "Name"))
    print("Chromosome:", get_text(gene, "Chromosome"))
    print("MapLocation:", get_text(gene, "MapLocation"))
    print("Description:", get_text(gene, "Description"))

    if get_text(gene, "Name") == "BRCA1":
        gene_id = str(gene.attributes["uid"])

if gene_id is None:
    raise RuntimeError("Nie udało się ustalić ID genu BRCA1.")

print("\nID genu BRCA1:", gene_id)


# krok 2 - rekordy omim zwiazane z genem
print("\n2. OMIM")

try:
    omim_ids = get_linked_ids("gene", "omim", gene_id)

    if omim_ids:
        omim_records = read_entrez(
            Entrez.esummary,
            db="omim",
            id=omim_ids
        )

        for record in omim_records:
            print(get_text(record, "Title", str(record)))
    else:
        print("Brak dostępnych powiązanych rekordów OMIM.")
except Exception as error:
    print("Nie udało się pobrać OMIM:", error)


# krok 3 - bialka powiazane z genem
print("\n3. PROTEIN")

protein_ids = get_linked_ids("gene", "protein", gene_id)

print("Liczba powiązanych rekordów:", len(protein_ids))

if protein_ids:
    protein_summaries = read_entrez(
        Entrez.esummary,
        db="protein",
        id=protein_ids[:100]
    )

    for protein in protein_summaries:
        accession = get_text(protein, "AccessionVersion")
        gi = protein.get("Gi")

        if gi:
            print("GI:", gi, "| Accession:", accession)
        else:
            print("GI: niedostępne | Accession:", accession)


# krok 4 - pobranie bialka o id 121949022 w formacie genbank
print("\n4. BIAŁKO 121949022")

try:
    with Entrez.efetch(
        db="protein",
        id="121949022",
        rettype="gb",
        retmode="text"
    ) as handle:
        protein_record = SeqIO.read(handle, "genbank")

    print("Nazwa:", protein_record.description)
    print("Sekwencja aminokwasów:")
    print(protein_record.seq)

except Exception as error:
    print("Nie udało się pobrać rekordu:", error)


# krok 5 - warianty snp powiazane z genem
print("\n5. SNP")

try:
    snp_ids = get_linked_ids(
        "gene",
        "snp",
        gene_id,
        term="Homo sapiens"
    )

    print("Liczba znalezionych rekordów SNP:", len(snp_ids))

    if snp_ids:
        with Entrez.efetch(
            db="snp",
            id=snp_ids[:50],
            retmode="xml"
        ) as handle:
            xml_data = handle.read()

        root = ET.fromstring(xml_data)

        def find_xml_value(element, field):
            for child in element.iter():
                if child.tag.split("}")[-1].upper() == field:
                    if child.text:
                        return child.text.strip()
            return "brak danych"

        records = [
            element for element in root.iter()
            if element.tag.split("}")[-1] == "DocumentSummary"
        ]

        if not records:
            records = [
                element for element in root.iter()
                if element.tag.split("}")[-1] in ("Rs", "SNP")
            ]

        for record in records[:50]:
            print(
                "SNP_ID:", find_xml_value(record, "SNP_ID"),
                "| SNP_CLASS:", find_xml_value(record, "SNP_CLASS"),
                "| GENES:", find_xml_value(record, "GENES"),
                "| CHRPOS:", find_xml_value(record, "CHRPOS")
            )

except Exception as error:
    print("Nie udało się pobrać danych SNP:", error)


# krok 6 - wnioski

print("\n6. WNIOSKI")
print(
    "BRCA1 jest genem człowieka związanym z naprawą DNA. "
    "Powiązania z bazami OMIM, Protein i dbSNP pozwalają "
    "analizować jego znaczenie kliniczne, produkty białkowe "
    "oraz zmienność genetyczną."
)
