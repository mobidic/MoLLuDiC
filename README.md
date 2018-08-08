# MoLLuDIC
- [MoLLuDiC - a CNV Exome pipeline adapted from clamms ](#MoLLuDiC)
	- [Overview](#overview)
		- [Main Workflow](#main-workflow)
		- [Citing MoLLuDiC](#citing-molludic)
		- [Input](#input)
		- [Output](#output)
	- [Installation](#installation)
		- [Requirements](#requirements)
	- [MoLLuDiC command](#molludic-command)
		- [MoLLuDiC Container](#molludic-container)
	- [Troubleshooting](#troubleshooting)

--------------------------------------------------------------------------------
<img src="img/molludic.png" width="150">

The **MoBiDiC capture CNV** calling, annotation and interpretation tools based on clamms workflow [https://github.com/rgcgithub/clamms](https://github.com/rgcgithub/clamms). 

This workflow is easily accessible with a Singularity image, low memory usage, handle batch effects and annotate with useful databases for human diagnostics.

# Overview

NV calling on NGS Capture library technologies are difficult to implement in bioinformatics pipeline and their performance are doubtful due to capture bias. 
Moreover, since there is no standard for CNV caller output, CNV are difficult to annotate. 

We propose an **exportable WDL workflow** based on open source tools and open source MoBiDiC script for CNV calling, annotation and interpretation interface. Performance of MoLLuDIC are increasing with data collection enlargement. 

- MoLLuDiC can be easily installed with a **singularity** container.
- MoLLuDiC is powered by **WDL** and **Cromwell** from Broad Institute. This pipeline are adaptable and easy-to-use via JSON input file .
- CNV calling and removing batch effects are based on **clamms** [https://github.com/rgcgithub/clamms](https://github.com/rgcgithub/clamms).
- CNV **familial segregation** is made with bedtools. 
- CNV annotation is based on **bedtools and MoBiDiC-made master annotation file**.
- CNV call could be interpreted in TSV viewer and SNV/CNV common interpretation could be realized via **Captain ACHAB** https://github.com/mobidic/Captain-ACHAB.

## Main workflow

![molludic workflow description](img/molludic_workflow.svg)

molludic.sh
- Select mode : "from scratch" to process all your data or "routine" to process new samples
- Select a library : capture from panel to whole exome sequencing
- Remove batch effect : Selection of a cluster of X most identical sample within 7 technical parameters (SEX, AT and GC Dropout, Mean Insert Size, Percentage of targeted base covered with 10X and 50X...) via a statistical method named k-d Tree.
The X number of sample is scalable depending on your data (KNN).  
- Remove relatives from CNV calling : add a family list file (tabulated file with sample identifier) remove relative from calling
- CNV calling : Read Depth Coverage statistical analysis
- CNV familial segregation : bedtools intersect and merge with all CNV call from relatives
- CNV annotation : bedtools intersect with a master annotation file containing cytoband, OMIM, ExAc CNV population frequency and metrics, in silico predictions tools. Instructions for master annotation file creation are described below. 


## Citing MoLLuDiC

> MoLLuDiC : Exportable CNV calling, annotating and interpretating workflow for NGS Capture sequencing (2019).  https://github.com/mobidic/MoLLuDiC

## Input

- Library bed file
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

## Requirements 

- Linux OS
- Cromwell
- C
- git
- python 3
- bedtools (v2.27.1)
- R software and the FNN package install.packages("FNN")

# MoLLuDiC command

```bash
  MoLLuDiC (version ${VERSION}) is a CNV workflow for calling and annotation !
  Usage : /.molludic.sh
 General arguments : 
      help : show this help message
    -v : decrease of increase verbosity level (ERROR : 1 | WARNING : 2 | INFO [default] : 3 | DEBUG : 4)

MoLLuDiC is composed of several functions. You print help for each module by typing help after function name.
    Example : ./molludic.sh install help
List of MoLLuDiC's functions : 
  dirpreparation <OPTION> : Create folders to use correctly clamms
  install <CLAMM_DIRECTORY> : install Clamms in specific directory
  mapinstall <CLAMM_DIRECTORY> <BigWigToWig_PATH> : create Mapability bed
  windowsBed <CLAMMS_DIRECTORY> <INSERT_SIZE> <INTERVALBEDFILE> <REFFASTA> <CLAMMS_SPECIAL_REGIONS> <LIBRARY_DIRECTORY> : run clamms annotate windows
  normalizeFS <COVERAGE_PATH> <CLAMMS_DIRECTORY> <WINDOWS_BED> <LIBRARY_DIRECTORY> : normalize bed files from scratch
  normalize <CLAMMS_DIRECTORY> <SAMPLEID> <CLAMMSCOVERAGEFILE> <WINDOWS_BED> <LIBRARY_DIRECTORY> : normalize one bed file
  metricsMatrixFS <LIBRARY_DIRECTORY> <HS_FOLDER> <PYTHON_PATH> <MATCH_METRICS> : create kd tree metrics from scratch
  metricsMatrix : <LIBRARY_DIRECTORY> <SAMPLEID> <HSMETRICSTXT> <INSERT_SIZE_METRICS_TXT> <PYTHON_PATH> <MATCH_METRICS> : create kd tree metric for 1 sample
  removeRelatives <ALLKDTREE> <FAMILYLIST> <LIBRARY_DIRECTORY> : remove relatives from all kd tree file
  makekdtree <RSCRIPT_PATH> <RSCRIPT_FILE> <KNN> <ALL_TREE> <LIBRARY_DIRECTORY> <FROM_SCRATCH> : use Rscript to do kd tree
  cnvCallingFS <CLAMMS_DIRECTORY> <LIBRARY_DIRECTORY> <LIST_KDTREE> <WINDOWS_BED> <KNN> : do calling from scratch
  cnvCalling <CLAMMS_DIRECTORY> <LIBRARY_DIRECTORY> <NORMCOVBED> <LIST_KDTREE> <WINDOWS_BED> <KNN> : do calling for 1 sample
  annotation <LIBRARY_DIRECTORY> <SAMPLEID> <BEDTOOLS_PATH> <HGBED> <HEADER_FILE> <CNV_BED> <DAD> (optional) <MUM> (optional) : annotate cnv bed file
```
## MoLLuDiC Container

Soon, you will be able to launch MoLLuDiC via a singularity container.

# Troubleshooting

## Using GATK

GATK needs that bed file contains "chr" before chromosome number. Clamms needs no chr before chromosome number. Be careful with your bed data ! 

```bash
sed 's/^chr//g' yourbedwithchr > bedforclamms.bed"
sed 's/^/chr/g' yourbedwithoutchr.bed > bedforgatk.bed"
```

## For panel capture

Creation of windowsBed need that your capture library.bed got the same chromosome that in the hg19.fa genome.
To create a specific hg19.fa genome without selected chromosome, please find a shell script that should work (example here with chromosome 13, 21 and 22 removed).

```bash
sed '/chr13/,/chr14/{//!d}' /usr/local/share/refData/genome/hg19/hg19.fa | grep -v "chr13" |  sed '/chr21/,/chr22/{//!d}' |  sed '/chr22/,/chrX/{//!d}' | grep -v "chr21" | grep -v "chr22" | sed 's/chr//g' > hg19_moins132122_nochr.fa

```
--------------------------------------------------------------------------------

**Montpellier Bioinformatique pour le Diagnostique Clinique (MoBiDiC)**

*CHU de Montpellier*

France

![MoBiDiC](logos/logo-mobidic.png)

[Visit our website](https://neuro-2.iurc.montp.inserm.fr/mobidic/)

--------------------------------------------------------------------------------
