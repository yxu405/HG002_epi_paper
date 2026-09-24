# Haplotype-resolved DiMeLo-seq maps centromeric chromatin in a complete diploid human genome

Code and analysis workflows accompanying the published study by Xu et al. in
*Cell Genomics*.

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.xgen.2026.101324-blue)](https://doi.org/10.1016/j.xgen.2026.101324)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## Overview

This study combines the complete diploid T2T-HG002 genome with ultra-long,
adaptive-sampling DiMeLo-seq to map CENP-A, H3K9me3, and CpG methylation across
haplotype-resolved human centromeres. The analyses show that CENP-A occupies
multiple discrete domains within centromere dip regions and that changes in DNA
methylation remodel centromeric chromatin boundaries while broadly preserving
CENP-A dosage between homologous chromosomes.

This repository contains the scripts, notebooks, browser-track workflows, and
supplemental tools used to process the sequencing data and reproduce analyses
for Figures 1 through 3. Large sequencing files are not stored in Git and must
be obtained through the data-access information associated with the publication.

## My contribution

As first author, I developed and applied computational workflows for processing
DiMeLo-seq reads, quantifying modified-base signals, characterizing centromere
dip regions, and generating publication figures. My contributions represented
in this repository include the figure-specific analyses, modified-base density
workflows, CDR analysis, and organization of the code into documented and
reusable components.

## Paper

Xu Y, Loucks H, Menendez J, et al. Haplotype-resolved DiMeLo-seq maps
centromeric chromatin in a complete diploid human genome. *Cell Genomics*.
2026;6(8):101324.

- [Published article](https://doi.org/10.1016/j.xgen.2026.101324)
- [PubMed record](https://pubmed.ncbi.nlm.nih.gov/42561950/)

## Analysis index

| Analysis | Repository location | Purpose |
| --- | --- | --- |
| Figure 1 | [`figure1/`](figure1/) | Haplotype-resolved CenSat and CDR overview |
| Figure 2C | [`figure2/fig2c_profiles/`](figure2/fig2c_profiles/) | Aggregate CENP-A, mCpG, and H3K9me3 profiles |
| Figure 2D | [`figure2/fig2d_cdr_variability/`](figure2/fig2d_cdr_variability/) | Chromosome-level CDR signal variability |
| Figure 2E | [`figure2/fig2e_chromosome_comparison/`](figure2/fig2e_chromosome_comparison/) | Maternal and paternal chromosome comparisons |
| Figure 2G | [`figure2/fig2g_cenpa_enrichment/`](figure2/fig2g_cenpa_enrichment/) | Per-read modified-base densities and Gini thresholds |
| Figure 3 | [`figure3/`](figure3/) | Passage- and cell-state-associated changes in centromeric chromatin |
| Modified-base density pipeline | [`intermediate_scripts/`](intermediate_scripts/) | Windowed 6mA and CpG 5mC density calculation |
| Supplemental analyses | [`Supplemental/`](Supplemental/) | CDR calling and external analysis workflows |
| DiMeLo-seq raw-data processing | [`Supplemental/HG002_dimelo_raw_data_processing/`](Supplemental/HG002_dimelo_raw_data_processing/) | POD5 basecalling, alignment, filtering, and BAM preparation |
| Genome browser hub | [`hub/`](hub/) | UCSC browser configuration and track-building workflows |

## Repository structure

```text
HG002_epi_paper/
├── figure1/
│   ├── Figure1A.py
│   ├── HG002v1.1_ONT_hmmCDR_subCDR.hc.merge100kb.bed
│   └── hg002v1.1.cenSatv2.0.filteredCT_CHRY.bed
├── figure2/
│   ├── fig2c_profiles/
│   ├── fig2d_cdr_variability/
│   ├── fig2e_chromosome_comparison/
│   └── fig2g_cenpa_enrichment/
│       ├── notebooks/
│       ├── src/
│       ├── README.md
│       ├── requirements.txt
│       └── run_analysis.py
├── figure3/
├── intermediate_scripts/
│   ├── README.md
│   └── dimelo_density.py
├── Supplemental/
│   ├── cdr_caller_v2.py
│   └── external Git submodules
├── hub/
├── LICENSE
└── README.md
```

## Getting started

Clone the repository and initialize its external workflow submodules.

```bash
git clone --recurse-submodules \
  https://github.com/yxu405/HG002_epi_paper.git

cd HG002_epi_paper
```

If the repository has already been cloned, initialize the submodules separately.

```bash
git submodule update --init --recursive
```

Create a Python environment for the common analysis dependencies.

```bash
python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install \
  biopython \
  jupyter \
  matplotlib \
  numpy \
  pandas \
  pysam \
  scikit-learn \
  scipy \
  seaborn \
  tabulate
```

The Figure 2C notebook additionally requires the DiMeLo analysis package used
for the study. Command-line tools for raw-read processing and browser-track
construction are described in their respective workflow directories.

## Reproducing Figure 1

The two small BED inputs required by `Figure1A.py` are included in the
repository. Run the following commands from the repository root.

```bash
mkdir -p Fig1A

python figure1/Figure1A.py \
  --bed figure1/hg002v1.1.cenSatv2.0.filteredCT_CHRY.bed \
  --cdr figure1/HG002v1.1_ONT_hmmCDR_subCDR.hc.merge100kb.bed
```

The script writes the chromosome overview and legend into `Fig1A/`.

## Reproducing Figure 2

The Figure 2C, 2D, and 2E analyses are provided as notebooks in their
panel-specific directories. Start Jupyter from the repository root.

```bash
jupyter lab
```

These notebooks require the sequencing-derived inputs described in their
configuration cells. Paths from the original cluster execution should be
replaced with paths appropriate for the local or computing-cluster environment.

### Figure 2G reusable workflow

Figure 2G includes a modular command-line implementation under
[`figure2/fig2g_cenpa_enrichment/`](figure2/fig2g_cenpa_enrichment/). Install its
requirements and display the complete interface.

```bash
cd figure2/fig2g_cenpa_enrichment
python -m pip install -r requirements.txt
python run_analysis.py --help
```

An example invocation is shown below.

```bash
python run_analysis.py \
  --reference /path/to/hg002v1.1.fasta \
  --cdr-bed /path/to/cdr_regions.bed \
  --active-bed /path/to/active_arrays.bed \
  --sample CENPA=/path/to/cenpa.bam \
  --sample CENPC=/path/to/cenpc.bam \
  --sample H3K9me3=/path/to/h3k9me3.bam \
  --output-dir results/fig2g \
  --segment-size 5000 \
  --non-cdr-limit 15
```

See the [Figure 2G README](figure2/fig2g_cenpa_enrichment/README.md) for input
requirements, scientific assumptions, and output descriptions.

## Reproducing Figure 3

Figure 3 notebooks analyze changes associated with extended cell passage and
pluripotent-cell reprogramming.

| Notebook | Paper panels |
| --- | --- |
| [`island_foldchange_fig3B.ipynb`](figure3/island_foldchange_fig3B.ipynb) | Figure 3B |
| [`flanking_foldchange_box_fig3C.ipynb`](figure3/flanking_foldchange_box_fig3C.ipynb) | Figure 3C |
| [`Density_dot_lineplot_fig3D_E.ipynb`](figure3/Density_dot_lineplot_fig3D_E.ipynb) | Figures 3D and 3E |

The notebooks require sequencing-derived density tables and BED files that are
not stored in Git. Configure the input and output paths before running the
notebooks in Jupyter.

## Modified-base density pipeline

`intermediate_scripts/dimelo_density.py` computes windowed modified-base
densities from an indexed BAM file and BED regions.

```bash
python intermediate_scripts/dimelo_density.py \
  --bam /path/to/input.bam \
  --bed /path/to/regions.bed \
  --ref /path/to/hg002v1.1.fasta \
  --mod-tag A \
  --threshold 0.5 \
  --output results/density.tsv \
  --threads 4 \
  --window-size 1000
```

Use `--mod-tag A` for 6mA or `--mod-tag CG` for CpG 5mC. Additional details
are available in the
[intermediate workflow README](intermediate_scripts/README.md).

## Raw DiMeLo-seq processing

The
[`HG002_dimelo_raw_data_processing`](Supplemental/HG002_dimelo_raw_data_processing/)
submodule contains the SLURM workflow used to convert raw POD5 data into
analysis-ready modified-base BAM files. The workflow covers Dorado basecalling,
Winnowmap alignment, modification-score filtering, primary-alignment filtering,
sorting, indexing, and preparation of modified-base tags for downstream
DiMeLo-seq analysis.

```bash
sbatch dimelo_raw_data_processing.slurm \
  SAMPLE_NAME \
  POD5_DIRECTORY \
  OUTPUT_DIRECTORY \
  REFERENCE_FASTA \
  REPETITIVE_KMER_FILE \
  DORADO_MODEL
```

The SLURM resource directives and software paths may require adjustment for a
different computing cluster.

## Supplemental tools and submodules

The supplemental directory includes the CDR caller and pinned versions of
external tools used in the study.

| Tool | Purpose | Upstream repository |
| --- | --- | --- |
| `cdr_caller_v2.py` | Identification and scoring of centromere dip regions | This repository |
| `alphaAnnotation` | Centromeric satellite annotation | [kmiga/alphaAnnotation](https://github.com/kmiga/alphaAnnotation) |
| `horhap_tool` | HOR haplotype annotation | [fedorrik/horhap_tool](https://github.com/fedorrik/horhap_tool) |
| `superHOR_HG002` | HG002 superHOR analysis | [fedorrik/superHOR_HG002](https://github.com/fedorrik/superHOR_HG002) |
| `HumAS-HMMER` | Alpha-satellite monomer annotation | [enigene/HumAS-HMMER](https://github.com/enigene/HumAS-HMMER) |
| `HumAS-HMMER_for_AnVIL` | HumAS-HMMER workflow for AnVIL | [fedorrik/HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) |
| `chm13_hsat` | Classical human-satellite k-mer database | [altemose/chm13_hsat](https://github.com/altemose/chm13_hsat) |
| `HG002_dimelo_raw_data_processing` | Raw DiMeLo-seq processing | [yxu405/HG002_dimelo_raw_data_processing](https://github.com/yxu405/HG002_dimelo_raw_data_processing) |

## Genome browser hub

The [`hub/`](hub/) directory contains UCSC browser-hub configuration files,
AutoSQL schemas, and track-building scripts for HG002v1.1. These workflows rely
on UCSC command-line utilities and large sequencing-derived inputs that are not
stored in Git.

## Data and reproducibility

- Small annotation files needed for Figure 1 are included in this repository.
- Raw POD5 files, BAM files, reference FASTA files, and large intermediate
  results are excluded because of their size.
- Data access and accession information should be taken from the published
  article and its data-availability statement.
- BAM inputs must be coordinate sorted and indexed where required.
- Reference FASTA inputs must be indexed where required.
- Legacy notebooks preserve the original analysis record and may contain
  cluster-specific paths that must be configured for another environment.

The command-line Figure 2G workflow is the most portable entry point for
reviewing the modified-base analysis architecture without opening the original
notebooks.

## Citation

If this repository or its workflows are useful in your work, please cite

> Xu Y, Loucks H, Menendez J, et al. Haplotype-resolved DiMeLo-seq maps
> centromeric chromatin in a complete diploid human genome. *Cell Genomics*.
> 2026;6(8):101324. https://doi.org/10.1016/j.xgen.2026.101324

## License

This repository is distributed under the terms of the [MIT License](LICENSE).
External submodules retain their own licenses and attribution requirements.

## Questions and feedback

For questions about the code or reproducibility, open an issue in this
repository with the relevant figure, script, command, and error message.
