import pandas as pd
from itertools import combinations


# Input files
FILES = {
    "Phenotype1": r"Path/To/File1.xlsx",
    "Phenotype2": r"Path/To/File2.xlsx",
    "Phenotype3": r"Path/To/File3.xlsx"
}

COL_GENE = "GENE_SYMBOL"
COL_ENTREZ = "ENTREZ_ID"


def load_data(files: dict) -> dict:
    data = {}
    for phenotype, path in files.items():
        df = pd.read_excel(path)
        df = df[[COL_GENE, COL_ENTREZ]].drop_duplicates(subset=[COL_GENE])
        data[phenotype] = df
    return data


def build_gene_map(data: dict) -> dict:
    gene_map = {}

    for phenotype, df in data.items():
        for _, row in df.iterrows():
            gene = row[COL_GENE]
            entrez = row[COL_ENTREZ]

            if gene not in gene_map:
                gene_map[gene] = {
                    "ENTREZ_ID": entrez,
                    "phenotypes": set()
                }

            gene_map[gene]["phenotypes"].add(phenotype)

    return gene_map


def generate_overlap_table(gene_map: dict) -> pd.DataFrame:
    records = []

    for gene, info in gene_map.items():
        if len(info["phenotypes"]) >= 2:
            records.append({
                "GENE_SYMBOL": gene,
                "ENTREZ_ID": info["ENTREZ_ID"],
                "PHENOTYPES": ", ".join(sorted(info["phenotypes"]))
            })

    return pd.DataFrame(records)


# Compute pairwise and higher-order overlaps
def compute_overlap_counts(gene_map: dict) -> dict:
    counts = {}

    for info in gene_map.values():
        phenotypes = sorted(info["phenotypes"])

        if len(phenotypes) >= 2:
            for r in range(2, len(phenotypes) + 1):
                for combo in combinations(phenotypes, r):
                    counts[combo] = counts.get(combo, 0) + 1

    return counts


def save_outputs(overlap_df: pd.DataFrame, counts: dict):
    overlap_df.to_excel("Gene_Overlap_Phenotypes.xlsx", index=False)

    with open("Gene_Overlap_Report.txt", "w", encoding="utf-8") as f:
        total = sum(counts.values())
        f.write(f"Total genes with overlap (>=2 phenotypes): {total}\n\n")

        f.write("Triple overlap (3 phenotypes):\n")
        for k, v in counts.items():
            if len(k) == 3:
                f.write(f"{' + '.join(k)}: {v} genes\n")

        f.write("\nPairwise overlaps:\n")
        for k, v in counts.items():
            if len(k) == 2:
                f.write(f"{' + '.join(k)}: {v} genes\n")


def main():
    data = load_data(FILES)
    gene_map = build_gene_map(data)

    overlap_df = generate_overlap_table(gene_map)
    counts = compute_overlap_counts(gene_map)

    save_outputs(overlap_df, counts)

    print("Completed.")


if __name__ == "__main__":
    main()