
## Repository Structure

```
HG002_epi_paper/
├── figure1/                         # Scripts for Figure 1
│   └── Figure1A.py
│
├── figure2/                         # Scripts for Figure 2
│   ├── CENPA_clustering_vis_figure2.ipynb
│   ├── CENPA_mCpG_H3K9me3_profile_plots_fig2C.ipynb
│   ├── CENPA_CDR_variability_heatmap_fig2D.ipynb
│   ├── chr_dot_plot_fig2E.ipynb
│   ├── decision_tree_gini_threshold_fig2G_1.ipynb
│   └── single_mol_CDR_CENPA_enrichment_fig2G_2.ipynb
│
├── figure3/                         # Scripts for Figure 3
│   ├── island_foldchange_fig3B.ipynb
│   ├── flanking_foldchange_box_fig3C.ipynb
│   └── Density_dot_lineplot_fig3D_E.ipynb
│
└── intermediate_scripts/            # Data processing pipeline scripts
    └── dimelo_density.py
```

---

## Scripts

### Figure 1
| Script | Description |
|--------|-------------|
| `Figure1A.py` | Generates a genome-wide visualization of centromeric satellite annotations (CenSat) and centromere dip regions (CDRs) across all chromosomes for both maternal and paternal haplotypes. Takes a CenSat BED file and CDR predictions as input and outputs a haplotype-resolved ideogram plot with a color-coded legend. |

### Figure 2
| Script | Description |
|--------|-------------|
| `CENPA_clustering_vis_figure2.ipynb` | CENP-A clustering visualization |
| `CENPA_mCpG_H3K9me3_profile_plots_fig2C.ipynb` | Profile plots of CENP-A, mCpG, and H3K9me3 signals |
| `CENPA_CDR_variability_heatmap_fig2D.ipynb` | Heatmap of CENP-A CDR variability across chromosomes |
| `chr_dot_plot_fig2E.ipynb` | Chromosome-level dot plot |
| `decision_tree_gini_threshold_fig2G_1.ipynb` | Decision tree with Gini threshold for CENP-A classification |
| `single_mol_CDR_CENPA_enrichment_fig2G_2.ipynb` | Single-molecule CDR CENP-A enrichment analysis |

### Figure 3
| Script | Description |
|--------|-------------|
| `island_foldchange_fig3B.ipynb` | Island-level fold change analysis |
| `flanking_foldchange_box_fig3C.ipynb` | Flanking region fold change box plots |
| `Density_dot_lineplot_fig3D_E.ipynb` | Density dot and line plots |

### Intermediate / Pipeline Scripts
| Script | Description |
|--------|-------------|
| `dimelo_density.py` | Computes region-level modified base density from DiMeLo-seq BAM files. Tiles genomic regions into fixed-size windows and calculates 6mA or CpG 5mC methylation density per window. Outputs a TSV of chrom, coordinates, density, and coverage. Supports parallel processing. |

---

##Tools & Workflows

The following external tools and workflows were used in this study. Please refer to their respective repositories for installation and usage instructions:

| Tool | Description | Repository |
|------|-------------|------------|
| DiMeLo-seq raw data processing | Batch SLURM script for BAM conversion, alignment, and post-processing | [HG002_dimelo_raw_data_processing](https://github.com/yxu405/HG002_dimelo_raw_data_processing) |
| alphaAnnotation | CenSat annotation workflow | [alphaAnnotation](https://github.com/kmiga/alphaAnnotation) |
| HORhap | HOR haplotype annotation tool | [horhap_tool](https://github.com/fedorrik/horhap_tool) |
| superHOR annotation | SuperHOR annotation notebook for HG002 | [superHOR_HG002](https://github.com/fedorrik/superHOR_HG002) |
| HumAS-HMMER | Monomeric annotations of alpha satellite HORs | [HumAS-HMMER](https://github.com/enigene/HumAS-HMMER) |
| HumAS-HMMER for AnVIL | HumAS-HMMER adapted for the AnVIL platform | [HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) |
| chm13_hsat | Classical satellite kmer database | [chm13_hsat](https://github.com/altemose/chm13_hsat) |

---

## Dependencies

- Python 3.8+
- `pysam`
- `numpy`
- `pandas`
- `biopython`
- `matplotlib`
- `jupyter`

Install dependencies with:
```bash
pip install pysam numpy pandas biopython matplotlib jupyter
```

---

## Usage — `dimelo_density.py`

```bash
python intermediate_scripts/dimelo_density.py \
  --bam input.bam \
  --bed regions.bed \
  --ref reference.fasta \
  --mod-tag CG \
  --threshold 0.5 \
  --output output.tsv \
  --threads 4 \
  --window-size 1000
```

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--bam` | Input BAM file with modified base tags |
| `--bed` | BED file with regions of interest |
| `--ref` | Reference genome FASTA |
| `--mod-tag` | Modification type: `A` (6mA) or `CG` (CpG 5mC) |
| `--threshold` | Score threshold; values below are set to 0 |
| `--output` | Output TSV file path |
| `--threads` | Number of parallel worker processes (default: 1) |
| `--window-size` | Window size in bp for tiling (default: 1000) |

---

> Xu Y, Loucks H, Menendez J, et al. *Haplotype-resolved centromeric chromatin organization from a complete diploid human genome.* Cell Genomics (submitted).
