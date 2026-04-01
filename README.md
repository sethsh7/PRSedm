# PRSedm
[![PyPI Version](https://img.shields.io/pypi/v/prsedm.svg)](https://pypi.org/project/prsedm/)
[![Conda Version](https://img.shields.io/conda/v/sethsh7/prsedm.svg)](https://anaconda.org/sethsh7/prsedm)
[![Paper DOI](https://img.shields.io/badge/DOI-10.2337%2Fdc25--0142-orange)](https://doi.org/10.2337/dc25-0142)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17903985.svg)](https://doi.org/10.5281/zenodo.17903985)

![Graphical Abstract](https://github.com/sethsh7/PRSedm/raw/main/misc/graphical_abstract_v2.png)

## Overview

**PRSedm (Polygenic Risk Score Extension for Diabetes Mellitus)** is a flexible and extendable open-source package for efficient local and remote (All of Us, UK Biobank, etc.) generation of published Polygenic Risk Scores (PRS) for Diabetes Mellitus (DM) and related cardiometabolic phenotypes. 

PRSedm introduces a new parallelized "one-liner" method to generate standardized PRS and pPS for DM robust to variables such as genotyping method, quality control, and imputation panel.  

## Updates (v1.3.0)
- Improved multiallelic SNP handling and fixed related bugs
- Significant performance improvements from optimized SNP batching and parallelism
- Simplified command-line interface and argument structure
- Fixed bugs and syntax issues in the SNP database backend
- PRS metadata is now fetched automatically alongside the SNP database
- Per-PRS variant log files are now generated, capturing metrics such as INFO/R², missing variants, and allele frequency
- Renamed "imputation" feature to "estimate" to distinguish from genotype imputation
- Added support for custom proxy variant substitution via `--proxy`

## New proxy feature

PRSedm supports optional variant substitution via a user-supplied proxy file (`--proxy`).
Required format (whitespace-delimited):

<pre style="overflow-x: auto; white-space: pre;">
target_rsid target_contig_id target_position target_effect_allele sub_rsid sub_contig_id sub_position sub_effect_allele
rs12345     chr1             1234567         A                     rs54321  chr1           1234999      A
rs23456     chr2             7654321         G                     rs65432  chr2           7654000      G
</pre>

## Installation

### Dependencies

PRSedm requires the following packages:

- Python (>=3.9), Joblib (>=1.3.2), Pandas (>=2.2.3), Pysam (>=0.22.0), Numpy\* (2.x/1.x)

\*Build with 1.x when deploying to RAP platforms with 1.x dependencies.

### User Installation
PIP: `pip install prsedm`  
Anaconda: `conda install sethsh7::prsedm`  
Build from source: `python -m build`

## Usage
PRSEDM can be called from the command line:

```prsedm --vcf <path_to_vcf_file> [options]```

PRSEDM can be also be called from Python:
```python
import prsedm
df = prsedm.gen_dm(
    vcf=vcf,
    ... 
)
```
### Options
- `--vcf` *(required)*: Path to an indexed VCF or BCF file, or a text file mapping one VCF/BCF per contig.
- `--col`: Genotype column to score (default: `GT`, options: `GT` for WGS, `GP` for imputed data).
- `--build`: Genome build to use (default: `hg38`, options: `hg19`, `hg38`).
- `--scores` *(required)*: Comma-separated list of PRS to generate, e.g., `t1dgrs2-luckett25,t2dp-udler18`.
- `--estimate` *(optional)*: Path to indexed reference VCF/BCF, or text file mapping one VCF/BCF per contig. Used to estimate missing variants and enable normalization when variants are absent.
- `--ntasks` *(optional)*: Number of tasks to use (default: `1`).
- `--batch-size` *(optional)*: Number of variants per batch (default: `5000`).
- `--output`: Path to save the output file (default: `results.csv`).
- `--full`: Include individual variant scores with PRS name prepended.
- `--getsql`: Download or locate the PRS SQL database (`variants.db`) and metadata JSON (`prs_meta.json`) and exit.


### Single file per-chromosome loading
For `--vcf` and `--estimate` you can point to a single text file mapping per contiguous region formatted as such (whitespace delimited):

```text
chr1   file1.chr1.vcf.gz
chr2   file2.chr2.vcf.gz
...
```

## PRS Database and metadata
The database containing PRS designs and metadata is hosted at
https://zenodo.org/records/17903390 and downloaded automatically.
PRSedm will first check the environment variables `PRSEDM_SQL_PATH` and `PRSEDM_META_PATH` for custom databases.

## Research Analysis Platforms (RAP's)

Remote deployment to remote Research Analysis Platforms (RAP's) is possible via notebook wrappers:
- All of Us (WGS) - [Notebook Here](https://github.com/sethsh7/PRSedm/blob/main/notebooks/PRSedm-AllofUS-notebook.ipynb)
- UK Biobank (imputed WGS) - [Notebook Here](https://github.com/sethsh7/PRSedm/blob/main/notebooks/PRSedm-DNAnexus-WGS-notebook.ipynb)
- UK Biobank (imputed array) - [Notebook Here](https://github.com/sethsh7/PRSedm/blob/main/notebooks/PRSedm-DNAnexus-imputed-notebook.ipynb)

## List of available PRS
### Type 1 Diabetes

| Flag | Method | Variants | Description | PMID |
| --- | --- | ---: | --- | --- |
| `t1dgrs2-luckett25` | HLA Interaction + Partitioned | 67 | "GRS2x" updated PRS with widest compatibility and HLA-based risk pPS. | [40267362](https://pubmed.ncbi.nlm.nih.gov/40267362/) |
| `t1dgrs2-qu22` | HLA Interaction + Partitioned | 71 | Original "GRS2" PRS with the addition of 4 African ancestry SNPs from Onengut, proposed in Qu et al and utilized in eMERGE. | [34997821](https://pubmed.ncbi.nlm.nih.gov/34997821/) |
| `t1dgrs2-sharp21` | HLA Interaction + Partitioned | 67 | Version of "GRS2" PRS designed for "TOPMED-R2" from 2021 GitHub. | [35312757](https://pubmed.ncbi.nlm.nih.gov/35312757/) |
| `t1d-onengut19-afr` | Additive | 6 | African-ancestry PRS proposed by Onengut in 2019, updated for modern compatibility. | [30659077](https://pubmed.ncbi.nlm.nih.gov/30659077/) |
| `t1dgrs2-sharp19` | HLA Interaction + Partitioned | 67 | Original 1000 Genomes version of "GRS2" PRS as published, with limited modern compatibility. | [30655379](https://pubmed.ncbi.nlm.nih.gov/30655379/) |

### Type 2 Diabetes

| Flag | Method | Variants | Description | PMID |
| --- | --- | ---: | --- | --- |
| `t2d-suzuki24-prscsx-ma` | Additive + Partitioned | ~1 M | Full genome-wide multi-ancestry PRS for Suzuki (PRS-CSx meta). | [38374256](https://pubmed.ncbi.nlm.nih.gov/38374256/) |
| `t2d-suzuki24-prscsx-<ancestry>` | Additive + Partitioned | >500k | Full genome-wide ancestry-specific PRS for Suzuki (PRS-CSx), where `<ancestry>` is one of `eur`, `afr`, `eas`, `sas`, `safr`, or `his`. | [38374256](https://pubmed.ncbi.nlm.nih.gov/38374256/) |
| `t2dp-suzuki24-ma` | Additive + Partitioned | 1289 | Multiancestry weighted Suzuki T2D index variant PRS, and pPS from hard-clustering analyses. | [38374256](https://pubmed.ncbi.nlm.nih.gov/38374256/) |
| `t2dp-suzuki24-<ancestry>` | Additive + Partitioned | 1128 - 1285 | As above but weighted for specific ancestries `<eur/afr/safr/eas/sas/his>`. | [38374256](https://pubmed.ncbi.nlm.nih.gov/38374256/) |
| `t2dp-smith24-ma` | Additive + Partitioned | 353 | Multiancestry cluster-weighted Smith T2D index variant PRS, and pPS from soft-clustering analyses. | [38443691](https://pubmed.ncbi.nlm.nih.gov/38443691/) |
| `t2dp-smith24-<ancestry>` | Additive + Partitioned | 25 - 490 | As above but from ancestry-specific soft clustering `<eur/afr/eas/amr>`. | [38443691](https://pubmed.ncbi.nlm.nih.gov/38443691/) |
| `t2d-mahajan22-ma` | Additive | 338 | Older PRS from Mahajan et al composed of multiancestry index variants. | [35551307](https://pubmed.ncbi.nlm.nih.gov/35551307/) |
| `t2d-mahajan22-prscsx-eur` | Additive + Partitioned | >500k | Genome-wide European ancestry PRS for Mahajan (PRS-CSx). | [38374256](https://pubmed.ncbi.nlm.nih.gov/38374256/) |
| `t2dp-udler18` | Additive + Partitioned | 67 | T2D pPS from first soft-clustering analysis. | [30240442](https://pubmed.ncbi.nlm.nih.gov/30240442/) |

### Other

| Flag | Phenotype | Method | Variants | Description | PMID |
| --- | --- | --- | ---: | --- | --- |
| `cdgrs-sharp25` | Celiac Disease | HLA Interaction + Partitioned | 42 | Modernized Celiac disease PRS and pPS with similar model to "GRS2x", utilized for combined screening. | [32790217](https://pubmed.ncbi.nlm.nih.gov/32790217/) |

## Features

### HLA Interaction PRS (+GRS2x Update)

PRSedm features a complete algorithm for GRS which incorporate HLA interaction terms as previously published by us such as T1D-GRS2 (or just GRS2). A number of advancements have been added to improve the generation of HLA interaction, described as GRS2x.

#### HLA Type Estimation and LD Tiebreak

HLA alleles can be estimated by proxy (or tag) single nucleotide polymorphisms alone and predictions are output e.g. `(DR3-DQ2.5/DR3-DQ2.5)`. Due to imperfect proxy SNPs, >2 HLA calls can be made in interaction scores such as GRS2, and a probabilistic tiebreaker algorithm using HLA reference frequencies (Klitz et al) now resolves impossible numbers of calls without excluding any samples.

### Missing variant mean effect estimation (optional)

PRSedm optionally uses Hardy-Weinberg Equilibrium with a reference VCF/BCF legend (ensure you have variant frequency coded as `AF`, genotypes not required) to estimate the mean effect size for missing SNPs, handle missing variants, and enable static normalization.

- dbSNP hg38 - [TOPMED Bravo Freeze 8](https://legacy.bravo.sph.umich.edu/freeze8/hg38/downloads), or [NCBI](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/) (`AF` field added) are recommended.
- dbSNP hg19 - [1000 Genomes](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html), [Haplotype Reference Consortium](https://www.sanger.ac.uk/collaboration/haplotype-reference-consortium/), [NCBI GRCh37](https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh37p13/VCF/) recommended.

### Minimum and Maximum Normalization

PRSedm hardcodes static normalization of minimum and maximum potential risk contribution (no risk alleles vs all risk alleles) creating a scale of 0-1. Static normalization with estimation ensures that PRS values translate to a common relative risk scale across datasets. If variants are missing and no estimation reference is supplied, normalization is skipped automatically and a warning is emitted.

## Development

Developed and maintained by Seth A. Sharp ([ssharp@stanford.edu](mailto:ssharp@stanford.edu)) at the [Translational Genomics of Diabetes](https://med.stanford.edu/genomics-of-diabetes.html), Stanford University, with collaboration from colleagues at the [University of Exeter](https://www.diabetesgenes.org) and [MGH/Broad Institute](https://www.udlerlab.org). Lu Zhang and Han Sun contributed to v1.1.0 onwards.

## Citation
If you use **PRSedm** in your research, please cite both the software release ([10.5281/zenodo.17903985](https://doi.org/10.5281/zenodo.17903985)) and the accompanying article ([10.2337/dc25-0142](https://doi.org/10.2337/dc25-0142)).

## License

This project is licensed under the **MIT License (Non-Commercial)**.  
- Academic, research, and personal use are allowed.  
- Commercial use is prohibited without prior permission.  

See the [LICENSE](https://raw.githubusercontent.com/sethsh7/PRSedm/main/LICENSE) file for full details.
