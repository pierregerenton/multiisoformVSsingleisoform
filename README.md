# multiisoformVSsingleisoform

Understand the information gain by the multiple-isoform gene annotation compared to single-isoform annotation.

## Goal

In the manuscript of the [geneidx](https://github.com/guigolab/geneidx), an analysis of the annotation of multiple vertebrates species has been done to understand the gain of functional information given by multiple-isoform annotation. That was done by extracting the longest isoform for coding genes and comparing the semantic similarity of GO terms of the longest isoform and the gene with all isoform. 

The purpose of this project is to extended this analysis to understand the mechanisms which genes can benefit from multiple-isoform annotation and why.

**Warning** : This repository is just a way to store and share script for people who want to learn about the project. It is not intended to be reproductible as the data used will not be provided in this git. However, the file [`data/get_data.md`](./data/get_data.md) will explain how generate these data locally.

## Requirement

If you want to get input data by yourselves :

- [samtools](https://github.com/samtools/samtools)
- [AGAT](https://github.com/NBISweden/AGAT)

To compute GO semantic similarity :

- [GOGO](https://github.com/zwang-bioinformatics/GOGO) : to run this software, your current working directory needs to be the installation GOGO directory

You'll also need [python3](https://www.python.org/downloads/) and some packages written in [`requierement.txt`](requierement.txt).

You can install all with :

```sh
pip install -r requierement.txt 
```

## Organization of the repo

- [`data`](data) : data location with [`data/get_data.md`](./data/get_data.md) with information to retreive them
- [`src`](src) : scripts.
    - [`src/choice_go_set`](src/choice_go_set) : subdirectory with a notebook to explore different dataset we can utilize for our analysis
    - [`src/analysis_notebook`](src/analysis_notebook) : some notebook used as a base for scripts


## Run the code

### Computation of the number of GO term with different annotation between files

To get an UpSetPlot of the number of genes which differ in term of GO annotation between input, you can run this command :

```sh
python3 src/number_genes_with_different_go_term_between_files.py -i data/pannzer_output/human.all.nr_off.out data/pannzer_output/human.long.nr_off.out data/pannzer_output/human.mane.nr_off.out -b data/pannzer_output/human.all.nr_off.out -o res/human_gene_count.pdf
```

- `-i` : pannzer output as input file (at least 2, or 1 with -`b` option)
- `-b` : path of a pannzer output where the best isoform will be kept (and this new element add to the plot)
- `-a` : infer GO term ancestry (longer)
- `-o` : name of the output file


python3 ./similarity_genes_between_files.py -i git_data/pannzer_output/human.all.nr_off.out git_data/pannzer_output/human.long.nr_off.out git_data/pannzer_output/human.mane.nr_off.out -g ~/software/GOGO/ -fb -o human_gene_sim