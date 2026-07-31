import pandas as pd
import os


# Input files

FILES = {

    "TWAS": r"Path/To/TWAS_File.txt",

    "MAGMA": r"Path/To/MAGMA_File.txt"

}


# Column names

TWAS_GENE = "genename"

TWAS_ID = "gene"


MAGMA_GENE = "genename"



# Output directory

OUTPUT_DIR = r"Path/To/Output/"


os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)



def load_twas(path):

    df = pd.read_csv(
        path,
        sep="\t"
    )

    df[TWAS_GENE] = (
        df[TWAS_GENE]
        .astype(str)
        .str.strip()
    )

    df = df.drop_duplicates(
        subset=[TWAS_GENE]
    )

    return df



def load_magma(path):

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[MAGMA_GENE]
    )


    df[MAGMA_GENE] = (
        df[MAGMA_GENE]
        .astype(str)
        .str.strip()
    )


    return df



def identify_consensus(twas, magma):

    magma_genes = set(
        magma[MAGMA_GENE]
    )


    consensus = twas[
        twas[TWAS_GENE].isin(magma_genes)
    ].copy()


    consensus = consensus.sort_values(
        "pvalue"
    )


    return consensus



def save_outputs(consensus):


    consensus.to_csv(

        OUTPUT_DIR +
        "Consensus_Magma_TWAS.txt",

        sep="\t",

        index=False,

        float_format="%.10e"

    )


    consensus[[
        TWAS_ID,
        TWAS_GENE

    ]].to_csv(

        OUTPUT_DIR +
        "Consensus_Genes.txt",

        sep="\t",

        index=False

    )



def save_report(consensus, twas, magma):


    with open(

        OUTPUT_DIR +
        "Consensus_Report.txt",

        "w",

        encoding="utf-8"

    ) as file:


        file.write(
            "Consensus MAGMA and TWAS-JTI analysis\n\n"
        )


        file.write(
            f"TWAS-JTI significant genes: {len(twas)}\n"
        )


        file.write(
            f"MAGMA significant genes: {len(magma)}\n"
        )


        file.write(
            f"Consensus genes: {len(consensus)}\n\n"
        )


        if len(consensus) > 0:


            file.write(
                "Top consensus genes:\n\n"
            )


            for _, row in consensus.head(20).iterrows():


                file.write(

                    f"{row[TWAS_GENE]}"

                )


                if "pvalue" in row:

                    file.write(

                        f" | pvalue={row['pvalue']:.10e}"

                    )


                if "Bonferroni" in row:

                    file.write(

                        f" | Bonferroni={row['Bonferroni']:.10e}"

                    )


                file.write("\n")



def main():


    twas = load_twas(
        FILES["TWAS"]
    )


    magma = load_magma(
        FILES["MAGMA"]
    )


    consensus = identify_consensus(
        twas,
        magma
    )


    save_outputs(
        consensus
    )


    save_report(
        consensus,
        twas,
        magma
    )


    print(
        "Consensus analysis completed."
    )



if __name__ == "__main__":

    main()