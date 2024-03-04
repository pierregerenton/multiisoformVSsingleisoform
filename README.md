# multiisoformVSsingleisoform

Understand the information gain by the multiple-isoform gene annotation compared to single-isoform annotation.

## Goal

In the manuscript of the [geneidx](https://github.com/guigolab/geneidx), an analysis of the annotation of multiple vertebrates species has been done to understand the gain of functional information given by multiple-isoform annotation. That was done by extracting the longest isoform for coding genes and comparing the semantic similarity of GO terms of the longest isoform and the gene with all isoform. 

The purpose of this project is to extended this analysis to understand the mechanisms which genes can benefit from multiple-isoform annotation and why.

**Warning** : This repository is just a way to store and share script for people who want to learn about the project. It is not intended to be reproductible as the data used will not be provided in this git. However, the file [`data/get_data.md`](./data/get_data.md) will explain how generate these data locally.

## Requirement

If you want to get input data by yourselves :

- [**samtools**](https://github.com/samtools/samtools)

