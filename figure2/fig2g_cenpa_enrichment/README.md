# Modified-base density analysis

This directory contains a refactored workflow for calculating per-read modified-
base densities in centromere dip regions and matched non-CDR intervals. It replaces
the original notebook's hard-coded paths, global BAM handle, repeated export code,
and cell-order dependencies with reusable Python modules and one command-line entry
point.

## What the workflow does

1. Reads CDR and active-array intervals from BED files
2. Subtracts CDR intervals from active arrays to define non-CDR intervals
3. Splits both groups into fixed-size windows
4. Maps high-confidence modified bases from each read to reference coordinates
5. Calculates per-read mA or CpG mC density in each window
6. Writes tidy density tables
7. Calculates chromosome-level variation
8. Fits overall and per-chromosome Gini thresholds

All genomic intervals use the BED coordinate convention. Starts are zero-based and
ends are exclusive.

## Files

- `config.py` defines validated analysis and sample settings
- `regions.py` reads, validates, subtracts, flanks, and segments genomic intervals
- `modifications.py` extracts mA and mC measurements from BAM modified-base tags
- `analysis.py` calculates variation and Gini thresholds
- `run_analysis.py` runs the complete workflow

## Requirements

- Python 3.9 or newer
- An indexed reference FASTA
- Coordinate-sorted and indexed BAM files containing MM and ML tags
- BED files describing CDRs and active arrays

Create a virtual environment and install the Python dependencies.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Create missing FASTA and BAM indexes with samtools.

```bash
samtools faidx data/reference/hg002v1.1.fasta
samtools index data/bam/cenpa.bam
samtools index data/bam/cenpc.bam
samtools index data/bam/h3k9me3.bam
```

## Recommended input layout

```text
project/
├── analysis/
│   ├── analysis.py
│   ├── config.py
│   ├── modifications.py
│   ├── regions.py
│   ├── requirements.txt
│   └── run_analysis.py
└── data/
    ├── bam/
    │   ├── cenpa.bam
    │   ├── cenpa.bam.bai
    │   ├── cenpc.bam
    │   ├── cenpc.bam.bai
    │   ├── h3k9me3.bam
    │   └── h3k9me3.bam.bai
    ├── reference/
    │   ├── hg002v1.1.fasta
    │   └── hg002v1.1.fasta.fai
    └── regions/
        ├── active_arrays.bed
        └── cdr_regions.bed
```

Large reference and BAM files should normally be excluded from Git and documented
with stable download or controlled-access instructions.

## Run the complete analysis

Run this command from the directory containing `run_analysis.py`. Replace the
example paths with the real local or cluster paths.

```bash
python run_analysis.py \
  --reference ../data/reference/hg002v1.1.fasta \
  --cdr-bed ../data/regions/cdr_regions.bed \
  --active-bed ../data/regions/active_arrays.bed \
  --sample CENPA=../data/bam/cenpa.bam \
  --sample CENPC=../data/bam/cenpc.bam \
  --sample H3K9me3=../data/bam/h3k9me3.bam \
  --output-dir ../results \
  --segment-size 5000 \
  --non-cdr-limit 15 \
  --min-mapq 0
```

The default run processes both mA and mC. To process one modification only, add
either `--modification mA` or `--modification mC`.

The original score cutoffs are retained as defaults. They can be changed with
`--ma-threshold` and `--mc-threshold`.

Use the built-in help to see every option.

```bash
python run_analysis.py --help
```

## Outputs

For each sample and modification, the workflow creates

- a CDR density CSV
- a non-CDR density CSV
- a per-chromosome threshold CSV

It also creates `analysis_summary.json`, which records modification settings,
record counts, variation summaries, and overall thresholds.

## Scientific assumptions to verify

- mA opportunities are aligned adenines
- mC opportunities are aligned CpG motifs
- modification scores use the 0 to 255 scale exposed by pysam
- duplicate, secondary, supplementary, and unmapped alignments are excluded
- reads without the requested modified-base tag are treated as missing, not zero
- only complete fixed-size windows are analyzed
- the first configured number of non-CDR windows is retained per chromosome

These assumptions are explicit so they can be reviewed against the experimental
protocol before publication or reuse.

## Reproducibility note

The full sequencing datasets may be too large or restricted for inclusion in Git.
The repository should still provide exact data accessions, expected filenames,
software versions, and a small synthetic example that exercises the complete
workflow.
