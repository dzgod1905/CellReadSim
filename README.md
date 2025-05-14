# CellReadSim

**CellReadSim** is a tool to simulate raw reads for **single-cell RNA-seq (scRNA-seq)** data.  
It was developed based on research around the question of how to simulate raw data for scRNA-seq experiments.

---

## Installation

**CellReadSim** is currently available for **Linux only**.  
The software is already compiled for Linux, and installation takes less than 2 minutes.  
Please follow the instructions below to install and use CellReadSim.

### Dependencies

The following packages are required to run CellReadSim:

- `gcc`
- `R` (with the following libraries):
  - `BiocManager`
  - `R.utils`
  - `GenomicFeatures`
  - `Biostrings`
  - `ShortRead`

> ⚠️ You only need to install **R** and **Rscript** manually.  
> The required R libraries will be automatically installed the first time the script runs.

---

## Quick Tutorial on CellReadSim

### Command to run:

```bash
Rscript gen.R "reference.Fasta" "transcript_expression_matrix.csv"
```

### Arguments

- `<reference.Fasta>` — Path to the reference cDNA FASTA file.
- `<transcript_expression.csv>` — CSV file containing the transcript-by-cell expression matrix.

The result will be saved in the "output" folder.  
Read1 will contain **barcodes and UMIs**, while Read2 will contain **fragmented transcript sequences**.

---

## Test Dataset

The test dataset includes:
- 10 cells
- 8 genes (26 total transcripts) with 69,177 simulated reads
- Transcript-level reads range from 0 to 5000
- Each cell has about 7000 reads on average

### To run the test dataset:

1. Unzip the FASTA file in the folder:

```bash
unzip input/test_dataset/Homo_sapiens.GRCh38.106.cdna.all.clean.zip -d input/test_dataset/
```

2. Run the command:
```bash
Rscript gen.R input/test_dataset/Homo_sapiens.GRCh38.106.cdna.all.clean.fa input/test_dataset/transcript_expression.csv input/test_dataset/ output/
```

