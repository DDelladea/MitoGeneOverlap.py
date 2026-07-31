import pandas as pd
import os


# Input files

FILES = {

    "Phenotype1": r"Path/To/File1.txt",
    "Phenotype2": r"Path/To/File2.txt",
    "Phenotype3": r"Path/To/File3.txt"

}


COL_GENE = "genename"
COL_ID = "gene"


# Output directory

OUTPUT_DIR = r"Path/To/Output/"


os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)



def load_data(files):

    data = {}

    for phenotype, path in files.items():

        df = pd.read_csv(
            path,
            sep="\t"
        )


        df[COL_GENE] = (
            df[COL_GENE]
            .astype(str)
            .str.strip()
        )


        df = df.drop_duplicates(
            subset=[COL_GENE]
        )


        data[phenotype] = df


    return data



def generate_overlap_groups(data):


    gene_sets = {

        phenotype: set(
            df[COL_GENE]
        )

        for phenotype, df in data.items()

    }


    phenotypes = list(gene_sets.keys())


    groups = {}


    shared_all = set.intersection(
        *gene_sets.values()
    )


    groups["Shared_All"] = shared_all



    for i in range(len(phenotypes)):

        for j in range(i + 1, len(phenotypes)):

            p1 = phenotypes[i]
            p2 = phenotypes[j]


            overlap = (
                gene_sets[p1]
                &
                gene_sets[p2]
            )


            overlap = overlap - shared_all


            groups[
                f"Shared_{p1}_{p2}"
            ] = overlap



    for phenotype in phenotypes:

        exclusive = (
            gene_sets[phenotype]
            -
            set.union(
                *[
                    gene_sets[p]
                    for p in phenotypes
                    if p != phenotype
                ]
            )
        )


        groups[
            f"Exclusive_{phenotype}"
        ] = exclusive



    return groups, gene_sets



def save_gene_file(genes, filename, data):


    rows = []


    for gene in sorted(genes):

        gene_id = None


        for df in data.values():

            hit = df[
                df[COL_GENE] == gene
            ]


            if not hit.empty:

                gene_id = hit.iloc[0][COL_ID]

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



def save_binary_matrix(gene_sets):


    all_genes = sorted(
        set.union(
            *gene_sets.values()
        )
    )


    matrix = pd.DataFrame({

        "genename": all_genes

    })


    for phenotype, genes in gene_sets.items():

        matrix[phenotype] = (

            matrix["genename"]
            .isin(genes)
            .astype(int)

        )


    matrix.to_csv(

        OUTPUT_DIR +
        "Consensus_BinaryMatrix.txt",

        sep="\t",

        index=False

    )



def save_report(groups, gene_sets):


    with open(

        OUTPUT_DIR +
        "Consensus_Overlap_Report.txt",

        "w",

        encoding="utf-8"

    ) as file:


        file.write(
            "Consensus gene overlap summary\n\n"
        )


        file.write(
            "Total genes per phenotype:\n"
        )


        for phenotype, genes in gene_sets.items():

            file.write(
                f"{phenotype}: {len(genes)}\n"
            )


        file.write(
            "\nShared genes:\n"
        )


        for name, genes in groups.items():

            if name.startswith("Shared"):

                file.write(
                    f"{name}: {len(genes)}\n"
                )


        file.write(
            "\nExclusive genes:\n"
        )


        for name, genes in groups.items():

            if name.startswith("Exclusive"):

                file.write(
                    f"{name}: {len(genes)}\n"
                )



def save_outputs(groups, data, gene_sets):


    for name, genes in groups.items():


        save_gene_file(

            genes,

            OUTPUT_DIR +
            "Genes_" +
            name +
            ".txt",

            data

        )


    save_binary_matrix(
        gene_sets
    )


    save_report(
        groups,
        gene_sets
    )



def main():


    data = load_data(
        FILES
    )


    groups, gene_sets = generate_overlap_groups(
        data
    )


    save_outputs(
        groups,
        data,
        gene_sets
    )


    print(
        "Consensus overlap analysis completed."
    )



if __name__ == "__main__":

    main()