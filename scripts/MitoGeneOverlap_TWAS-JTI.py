import pandas as pd


# Input files

FILES = {
    "Phenotype1": r"Path/To/File1.txt",
    "Phenotype2": r"Path/To/File2.txt",
    "Phenotype3": r"Path/To/File3.txt"
}


# Column names

COL_GENE = "genename"
COL_ID = "gene"


# Output directory

OUTPUT_DIR = r"Path/To/Output/"


def load_data(files):

    data = {}

    for phenotype, path in files.items():

        df = pd.read_csv(
            path,
            sep="\t"
        )

        df = df.drop_duplicates(
            subset=[COL_GENE]
        )

        data[phenotype] = df

    return data


def build_gene_map(data):

    gene_map = {}

    for phenotype, df in data.items():

        for _, row in df.iterrows():

            gene = row[COL_GENE]

            if gene not in gene_map:

                gene_map[gene] = {
                    "gene_id": row[COL_ID],
                    "phenotypes": set()
                }

            gene_map[gene]["phenotypes"].add(phenotype)

    return gene_map


def generate_overlap_groups(gene_map):

    groups = {
        "Shared_All": set(),

        "Shared_Phenotype1_Phenotype2": set(),
        "Shared_Phenotype1_Phenotype3": set(),
        "Shared_Phenotype2_Phenotype3": set(),

        "Exclusive_Phenotype1": set(),
        "Exclusive_Phenotype2": set(),
        "Exclusive_Phenotype3": set()
    }


    for gene, info in gene_map.items():

        phenotypes = info["phenotypes"]


        if len(phenotypes) == 3:

            groups["Shared_All"].add(gene)


        elif len(phenotypes) == 2:

            phenotype_list = sorted(phenotypes)

            key = "Shared_" + "_".join(phenotype_list)

            groups[key].add(gene)


        elif len(phenotypes) == 1:

            phenotype = list(phenotypes)[0]

            groups[
                f"Exclusive_{phenotype}"
            ].add(gene)


    return groups


def save_gene_file(genes, filename, data):

    rows = []

    for gene in sorted(genes):

        gene_id = None

        for df in data.values():

            result = df[
                df[COL_GENE] == gene
            ]

            if not result.empty:

                gene_id = result.iloc[0][COL_ID]
                break


        rows.append({
            "gene": gene_id,
            "genename": gene
        })


    pd.DataFrame(rows).to_csv(
        filename,
        sep="\t",
        index=False
    )


def save_full_file(genes, filename, data):

    output = []

    for phenotype, df in data.items():

        subset = df[
            df[COL_GENE].isin(genes)
        ].copy()

        subset["phenotype"] = phenotype

        output.append(subset)


    if output:

        result = pd.concat(output)

        result = result.sort_values(
            [COL_GENE, "phenotype"]
        )

        result.to_csv(
            filename,
            sep="\t",
            index=False,
            float_format="%.10e"
        )


def save_report(groups, data):

    with open(
        OUTPUT_DIR + "Overlap_Report.txt",
        "w",
        encoding="utf-8"
    ) as file:


        file.write("Gene overlap summary\n\n")

        file.write("Total genes per phenotype:\n")

        for phenotype, df in data.items():

            total = df[COL_GENE].nunique()

            file.write(
                f"{phenotype}: {total}\n"
            )


        file.write("\nShared genes:\n")

        for name, genes in groups.items():

            if name.startswith("Shared"):

                file.write(
                    f"{name}: {len(genes)}\n"
                )


        file.write("\nExclusive genes:\n")

        for name, genes in groups.items():

            if name.startswith("Exclusive"):

                file.write(
                    f"{name}: {len(genes)}\n"
                )


def save_outputs(groups, data):

    for name, genes in groups.items():

        save_gene_file(
            genes,
            OUTPUT_DIR + "Genes_" + name + ".txt",
            data
        )


        save_full_file(
            genes,
            OUTPUT_DIR + "Full_" + name + ".txt",
            data
        )


    save_report(
        groups,
        data
    )


def main():

    data = load_data(FILES)

    gene_map = build_gene_map(data)

    groups = generate_overlap_groups(
        gene_map
    )

    save_outputs(
        groups,
        data
    )

    print("Overlap analysis completed.")


if __name__ == "__main__":
    main()