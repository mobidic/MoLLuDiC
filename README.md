# MoLLuDIC
- [MoLLuDiC - a CNV Exome pipeline adapted from clamms ](#MoLLuDiC)
	- [Overview](#overview)
		- [Citing MoLLuDiC](#citing-molludic)
		- [Input](#input)
		- [Output](#output)
	- [Installation](#installation)
		- [Requirements](#requirements)
	- [Quick start](#quick-start)

--------------------------------------------------------------------------------
![MoLLuDiC logo](https://raw.githubusercontent.com/mobidic/MoLLuDiC/master/molludic.png =150x150)

The **MoBiDiC capture CNV** calling, annotation and interpretation tools based on clamms workflow [https://github.com/rgcgithub/clamms](https://github.com/rgcgithub/clamms). 

This workflow is easily accessible with a Singularity image, low memory usage, handle batch effects and annotate with useful databases for human diagnostics.

# Overview

NV calling on NGS Capture library technologies are difficult to implement in bioinformatics pipeline and their performance are doubtful due to capture bias. 
Moreover, since there is no standard for CNV caller output, CNV are difficult to annotate. 

We propose an **exportable WDL workflow** based on open source tools and open source MoBiDiC script for CNV calling, annotation and interpretation interface.

- MoLLuDiC can be easily installed with a **singularity** container.
- MoLLuDiC is powered by **WDL** and **Cromwell** from Broad Institute. This pipeline are adaptable and easy-to-use via JSON input file .
- CNV calling and removing batch effects are based on **clamms** [https://github.com/rgcgithub/clamms](https://github.com/rgcgithub/clamms).
- CNV annotation is based on **bedtools and MoBiDiC-made master annotation file**.
- CNV **familial segregation** is made with bedtools. 
- CNV call could be interpreted in TSV viewer and SNV/CNV common interpretation could be realized via **Captain ACHAB** https://github.com/mobidic/Captain-ACHAB.


## Citing MoLLuDiC

> MoLLuDiC : Exportable CNV calling, annotating and interpretating workflow for NGS Capture sequencing (2019).  https://github.com/mobidic/MoLLuDiC

## Input

- Metrics from Picard Tools (CollectHsMetrics and InsertSizeMetrics)
- Coverage from samtools bedcov or GATK DepthOfCoverage

## Output

- An annotated CNV calling bed file (with or without familial segregation)

# Installation

To download MoLLuDiC, please use git to download the most recent development tree.
Currently, the tree is hosted on github, and can be obtained via:

```bash
$ git clone https://github.com/mobidic/MoLLuDiC.git
```
## Create Singularity Image
**  First, build**
```bash
singularity build <filename.simg> Singulairity 
```
**Then run**
```bash
singularity run <filename.simg> -i workflow_inputs.json
```
**Singularity help**
```bash
singularity help <filename.simg>
```
## Requirements 

Need to install the FNN package install.packages("FNN")

### Install Singularity


# Quick start

## From Scratch workflow


## NGS routine workflow




