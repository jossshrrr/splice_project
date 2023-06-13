# Splice_Project

This project contains the code developed for my Bioinformatics Master's project. The primary aim of this project is to identify intronic variants that are predicted to impact splicing. 

**Please note that the code is currently in development. At present, the code relies on locally installed software, but we are in the process of containerizing it.**

This repository is primarily set up for my supervisors to inspect the overall structure of the pipeline.

## Table of Contents

- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)

## Description

This pipeline utilises several existing splice prediction tools. Variants are annotated with prediction scores from the tools below: 

- [SpliceAI](https://github.com/Illumina/SpliceAI)
- [CADD-Splice](https://github.com/kircherlab/CADD-scripts)
- [Pangolin](https://github.com/tkzeng/Pangolin)
- [SQUIRLS](https://github.com/TheJacksonLaboratory/Squirls)

A decision-tree model is then applied to the annotated variants. This ensemble approach was shown to outperform any individual tool when applied to test data.

## Installation

Currently, the code is not yet set up for installation and execution on systems other than the development environment. This section will be updated when the code is ready for general use.

## Usage

The pipeline is implemented in [*Nextflow*](https://github.com/nextflow-io/nextflow). The pipeline is run via the command below:

```sh
nextflow run variant_workflow.nf -c workflow_info.config
```

## Parameters

This workflow uses a configuration file called `workflow_info.config` to set various parameters. Below is a brief explanation of each parameter in the file:

| Parameter          | Description                                                                                                     | Example Value                                       |
|--------------------|-----------------------------------------------------------------------------------------------------------------|-----------------------------------------------------|
| `id`               | Identifier for the cohort.                                                                                      | `'Cohort_ID'`                                       |
| `input_vcf`        | Path to the input VCF file (compressed via bgzip and indexed via tabix).                                        | `'/path/to/variants.vcf.gz'`                        |
| `filter`           | A boolean indicating whether to apply gene-list filtering.                                                      | `true`                                              |
| `gene_bed`         | Path to the BED file containing genes to filter for.                                                            | `'/pipeline_data/EpilepsyGenes_v2022-09_230210.bed'`|
| `reference_fasta`  | Path to the reference FASTA file.                                                                               | `'/path/to/Homo_sapiens_assembly38.fasta'`          |
| `gnomad_count`     | Variants that have a gnomAD count > value are filtered out of the VCF.                                          | `100`                                               |
| `spliceai_distance`| The distance within which SpliceAI should analyze variants. Optional; if not provided, uses default distance.   | `500`                                               |
| `pangolin_distance`| The distance within which Pangolin should analyze variants. Optional; if not provided, uses default distance.   | `500`                                               |

## Additional Files Required

This workflow requires several additional files that are essential for its operation. Below is a list of these files along with a brief description, the file path, and the source where they can be obtained:

| Description                            | Path                                             | Source                                                      |
|----------------------------------------|--------------------------------------------------|-------------------------------------------------------------|
| Pangolin annotation file.              | `/pipeline_data/pangolin/gencode.v38.annotation.db` | [Pangolin data](https://www.dropbox.com/sh/6zo0aegoalvgd9f/AADWN_cGIWpvVN9BYJ37vGmZa?dl=0) |
| SpliceAI annotation file.              | `/pipeline_data/spliceai/gencode.v38.annotation.txt` | [SpliceAI-lookup data](https://spliceailookup-api.broadinstitute.org/annotations) |
| SQUIRLS command-line interface.        | `/pipeline_data/squirls/squirls-cli-2.0.0.jar`     | [Squirls GitHub](https://github.com/TheJacksonLaboratory/squirls) |
| SQUIRLS data.                          | `/pipeline_data/squirls/data`                      | [SQUIRLS data](https://squirls.readthedocs.io/en/master/setup.html#squirls-downloadable-resources) |
| CADD v1.6 Script.                      | `/pipeline_data/cadd/CADD.sh`                      | [CADD Script](https://github.com/kircherlab/CADD-scripts)   |
| VEP cache directory.                   | `/path/to/vep-cache`                              | [Ensembl website](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html) |
| gnomAD data.                           | `/pipeline_data/gnomad/gnomad.genomes.r3.0.sites.vcf.bgz` | [gnomAD FTP](https://gnomad.broadinstitute.org/downloads/#v3) |
| ClinVar data.                          | `/pipeline_data/clinvar/clinvar.vcf.gz`            | [ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) |

Please ensure that these files are in the correct locations specified above as they are essential for the proper functioning of the workflow.

## vcfanno Annotation

The pipeline utilizes the tool [*vcfanno*](https://github.com/brentp/vcfanno) to annotate variants with information from additional sources. *vcfanno* requires configuration files. All files below are in the directory `'/pipeline_data/vcfanno/'`:

| Parameter      | Description                                         | Example Value                  |
|----------------|-----------------------------------------------------|--------------------------------|
| `gnomAD`  | Path to the *vcfanno* configuration file for gnomAD.| `gnomad.toml` |
| `ClinV` | Path to the *vcfanno* configuration file for ClinVar.| `clinvar.toml`|


The pipeline also uses *vcfanno* to process splice prediction scores. This requires the additional configuration files written in both TOML and Lua scripting language below:

| Tool      | Lua File             | TOML File                   |
|-----------|----------------------|-----------------------------|
| Pangolin  | `pangolin.lua`       | `pangolin_postprocess.toml` |
| SpliceAI  | `spliceai.lua`       | `spliceai_postprocess.toml` |
| SQUIRLS   | `squirls.lua`        | `squirls_postprocess.toml`  |

Please note *CADD* is already provided as a single PHRED score during annotation, so no additional postprocessing is required for this score. 

## Acknowledgments

I would like to express my gratitude to my supervisors for their support throughout the development of this project.
