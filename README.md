# MitoGeneOverlap

A Python-based tool for identifying overlapping genes across multiple biological datasets.

## Overview

MitoGeneOverlap is a Python tool developed to identify shared and unique genes across different biological datasets. The tool was designed to support integrative analyses involving mitochondrial genes associated with complex phenotypes.

This repository contains different versions of the script adapted for distinct analytical steps, including:

* overlap analysis among phenotype-specific gene lists;
* integration of results obtained from different statistical approaches;
* identification of consensus genes;
* identification of shared consensus genes across multiple phenotypes.

The tool was developed and applied in the investigation of mitochondrial genes potentially involved in sarcopenia by integrating results from Gene-Based Analysis (GBA) and Transcriptome-Wide Association Study using Joint-Tissue Imputation (TWAS-JTI).

---

# Analysis workflow

The scripts available in this repository correspond to different stages of the analytical pipeline:

```
Significant mitochondrial genes GBA / TWAS-JTI
            |
            ↓
Phenotype-specific overlap analysis GBA / TWAS-JTI
            |
            ↓
Integration between MAGMA GBA and TWAS-JTI
            |
            ↓
Consensus mitochondrial genes
            |
            ↓
Overlap analysis among consensus genes
```

---

# Script versions

## 1. MitoGeneOverlap_GBA.py

### Purpose

Identifies overlapping genes among multiple phenotype-specific gene lists obtained from Gene-Based Analysis.

This version was designed for gene sets obtained after mitochondrial gene mapping and statistical filtering.

### Input files

Three tabular files containing gene identifiers.

Required columns:

```
GENE_SYMBOL
ENTREZ_ID
```

### Output files

The script generates:

* genes shared among all phenotypes;
* genes shared between phenotype pairs;
* phenotype-specific genes;
* overlap summary report.

---

## 2. MitoGeneOverlap_TWAS.py

### Purpose

Identifies overlapping genes among phenotype-specific gene lists obtained from TWAS-JTI analysis.

### Input files

Three tabular files containing TWAS-JTI significant genes.

Required columns:

```
gene
genename
```

Additional columns from TWAS-JTI results are preserved in complete output files.

### Output files

The script generates:

* shared genes among phenotypes;
* pairwise shared genes;
* exclusive genes;
* complete overlap tables;
* overlap summary report.

---

## 3. MitoGeneOverlap_Consensus.py

### Purpose

Identifies consensus genes detected by two independent analytical approaches:

* MAGMA Gene-Based Analysis;
* TWAS-JTI analysis.

Consensus genes are defined as genes present in both significant gene sets.

### Input files

Two files are required:

1. TWAS-JTI significant genes.
2. MAGMA significant genes.

### Output files

The script generates:

* complete consensus gene table;
* simplified consensus gene list;
* consensus summary report.

---

## 4. MitoGeneOverlap_ConsensusOverlap.py

### Purpose

Identifies shared consensus genes across multiple phenotypes after integration of MAGMA GBA and TWAS-JTI results.

### Input files

Three consensus gene lists.

Required columns:

```
gene
genename
```

### Output files

The script generates:

* genes shared among all phenotypes;
* genes shared between phenotype pairs;
* phenotype-specific consensus genes;
* binary presence/absence matrix;
* overlap summary report.

---

# Requirements

Python version:

```
Python ≥ 3.8
```

Required packages:

```
pandas
openpyxl
```

Install dependencies:

```bash
pip install pandas openpyxl
```

---

# Usage

1. Clone this repository.

2. Install dependencies.

3. Edit the input file paths inside the selected script:

```python
FILES = {
    "Phenotype1": "Path/To/File1",
    "Phenotype2": "Path/To/File2",
    "Phenotype3": "Path/To/File3"
}
```

4. Execute:

```bash
python MitoGeneOverlap_GBA.py
```

or the desired script version.

---

# Citation

If you use MitoGeneOverlap, please cite:

**------------**



