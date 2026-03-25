# MitoGeneOverlap.py
Python script to identify and quantify gene overlap between phenotypes, including triple and pairwise intersections.

## Usage

### 1. Requirements
- Python 3.x
- pandas

Install dependencies:
pip install pandas openpyxl

### 2. Input files

- Excel files containing gene lists for each phenotype
- Required columns:
  - GENE_SYMBOL
  - ENTREZ_ID

### 3. Configure input files

Edit the FILES dictionary in the script:

FILES = {
    "Phenotype1": r"Path/To/File1.xlsx",
    "Phenotype2": r"Path/To/File2.xlsx",
    "Phenotype3": r"Path/To/File3.xlsx"
}

### 4. Run the script

python GeneOverlapMap.py

### 5. Output files

The script generates:

- Gene_Overlap_Phenotypes.xlsx  
  -> Genes shared across two or more phenotypes

- Gene_Overlap_Report.txt  
  -> Summary of overlap counts:
     - Pairwise overlaps
     - Triple overlaps (if applicable)

### 6. Notes

- Each gene is considered once per phenotype (duplicates are removed)
- Gene identifiers must be consistent across input files
