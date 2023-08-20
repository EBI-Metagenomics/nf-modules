# Microbiome Informatics Nextflow modules & subworkflows

Welcome to the Microbiome Informatics Nextflow modules and subworkflows repository. This repository uses the same tools and conventions as [nf-core modules](https://nf-co.re/).

## Getting Started

To use the modules and subworkflows from this repository make sure you have [nf-core/tools](https://github.com/nf-core/tools) installed:

```bash
pip install nf-core
```

## Modules Usage

List the available modules:

```bash
nf-core modules --git-remote git@github.com:EBI-Metagenomics/nf-modules.git list remote
```

Install a module in your pipeline:

```bash
nf-core modules --git-remote git@github.com:EBI-Metagenomics/nf-modules.git install <tool>
```

Then import the desired module in your pipeline script:

```groovy
include { module_name } from '../modules/ebi-metagenomics/tool/main.nf'
```

## Subworkflows Usage

Subworkflows in this repository are self-contained and can be included in your main pipeline:

Install a module in your pipeline:

```bash
nf-core subworkflows --git-remote git@github.com:EBI-Metagenomics/nf-modules.git install <subworkflow>
```

Then import the desired module in your pipeline script:

```groovy
include { <subworkflow_name> } from '../subworkflows/ebi-metagenomics/<subworkflow_name>.nf'
```

## References

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
