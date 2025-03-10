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

## nf-core modules

The [nf-core](https://nf-co.re/) team supports a large number of high-quality modules, and our team contributes whenever we can. At the moment, the [nf-core tools](https://github.com/nf-core/tools/) don't support subworkflows that install modules from different repos ([#3083](https://github.com/nf-core/tools/pull/3083)). That is why we decided to copy some modules from nf-core into this repo (a nasty hack, but it works). The nf-core team has been making impressive progress on supporting this use case (subworkflows with modules from different repos), and we will remove the duplicated modules once they reach that point. In the meantime, you will find duplicated modules from nf-core here.

## References

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
