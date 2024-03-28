# multiisoformVSsingleisoform

Understand the information gain by the multiple-isoform gene annotation compared to single-isoform annotation.

## Goal

In the manuscript of the [geneidx](https://github.com/guigolab/geneidx), an analysis of the annotation of multiple vertebrates species has been done to understand the gain of functional information given by multiple-isoform annotation. That was done by extracting the longest isoform for coding genes and comparing the semantic similarity of GO terms of the longest isoform and the gene with all isoform. 

The purpose of this project is to extended this analysis to understand the mechanisms which genes can benefit from multiple-isoform annotation and why.

**Warning** : This repository is just a way to store and share script for people who want to learn about the project. It is not intended to be reproductible as the data used will not be provided in this git. However, the file [`data/get_data.md`](./data/get_data.md) will explain how generate these data locally.

## Summary

- [Goal](#goal)
- [Summary](#summary)
- [Requierement](#requirement)
- [Organization of the repo](#organization-of-the-repo)
- [Run the code](#run-the-code)
    - [Computation of the number of GO terms with different annotation between files](#computation-of-the-number-of-go-terms-with-different-annotation-between-files)
    - [Computation of GOGO Similarity of genes between files](#computation-of-gogo-similarity-of-genes-between-files)
    - [Print different plot to compare a reference annotation with an alternative one](#print-different-plot-to-compare-a-reference-annotation-with-an-alternative-one)
    - [Make a similarity table for each gene between two pannzer output](#make-a-similarity-table-for-each-gene-between-two-pannzer-output)
    - [Write exhaustive data table with metadata and similarity based on a multiple isoform annotation](#write-exhaustive-data-table-with-metadata-and-similarity-based-on-a-multiple-isoform-annotation)
    - [Perform GO enrichment analysis (GOEA) on previous results](#perform-go-enrichment-analysis-goea-on-previous-results)
    - [Evaluate isoforms diversity for each gene](#evaluate-isoforms-diversity-for-each-gene)

## Requirement

If you want to get input data by yourselves :

- [samtools](https://github.com/samtools/samtools)
- [AGAT](https://github.com/NBISweden/AGAT)

To compute GO semantic similarity :

- [GOGO](https://github.com/zwang-bioinformatics/GOGO) : to run this software, your current working directory needs to be the installation GOGO directory

To visualize the GOEA (GO Enrichment Analysis) :

- [Graphviz](https://www.graphviz.org/) : it's an open source graph visualization software.

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

### Computation of the number of GO terms with different annotation between files

To get an UpSetPlot of the number of genes which differ in term of GO annotation between input, you can run this command :

```sh
python3 src/number_genes_with_different_go_term_between_files.py -i data/pannzer_output/human.all.nr_off.out data/pannzer_output/human.long.nr_off.out data/pannzer_output/human.mane.nr_off.out -b data/pannzer_output/human.all.nr_off.out -o res/human_gene_count.pdf
```

- `-i` : pannzer output as input file (at least 2, or 1 with `-b` option)
- `-b` : path of a pannzer output where the best isoform will be kept (and this new element add to the plot)
- `-a` : infer GO term ancestry (longer)
- `-o` : name of the output file


python3 ./similarity_genes_between_files.py -i git_data/pannzer_output/human.all.nr_off.out git_data/pannzer_output/human.long.nr_off.out git_data/pannzer_output/human.mane.nr_off.out -g ~/software/GOGO/ -fb -o human_gene_sim

### Computation of GOGO similarity of genes between files

To get tables of mean gene similarity between files and a table of similarity between each gene for each pair, you can run this command :

```sh
python3 ./src/similarity_genes_between_files.py -i git_data/pannzer_output/human.all.nr_off.out git_data/pannzer_output/human.long.nr_off.out git_data/pannzer_output/human.mane.nr_off.out -g ~/software/GOGO/ -fb -o res/human_gene_sim
```

- `-i` : pannzer output as input file (at least 2, or 1 with `-b` option)
- `-g` : path of GOGO directory
- `-b` : path of a pannzer output where the best isoform will be kept (and this new element add to the plot)
- `-a` : infer GO term ancestry (longer)
- `-f` : filter gene_set to have only multiple-isoform gene used for similarity
- `-o` : name of the output file

Two files are generated :

- `output.txt` : a table with multiple column : `COMPARISON` is for the pair of annotation compared to get the similarity, `GENE_ID` is the gene id, `[BP/CC/MF]_SIMILARITY` is the similarity for each ontology (warning : NA is possible if GOGO wasn't able to compute similarity)
- `output.pdf` : 3 tables with mean similarity for each ontology and pair of input


### Print different plot to compare a reference annotation with an alternative one

This script prints plots to compare a reference annotation (usually a multi-isoform one) and alternative one which is generated.
There are two possibilities, the `longest` isoform (with the longest cds) and the `best` isoform (with the high number of GO term).
Plots generated include :
- Observed and expected similarities for each ontology (expected are calculated by randomly reassigned each transcript to a gene with conservation of the number of transcripts by gene)
- GOGO similarity by number of isoforms

```sh
python3 ./src/precise_analysis_of_one_multiisoform_annotation.py -i data/pannzer_output/human.all.nr_off.out -g ~/Software/GOGO/ -t long -f -o src/human_allVSlong_plots.pdf
```

- `-i` : reference pannzer output as input file
- `-g` : path of GOGO directory
- `-t` : `long` or `best` to choose if the alternative file keep only the longest or the best isoform
- `-f` : filter gene_set to have only multiple-isoform gene used for similarity
- `-o` : name of the output file


### Make a similarity table for each gene between two pannzer output

With two Pannzer output, one from a multiple-isoform annotation and the other from a single-isoform annotation, a table with the similarity for each gene between files is written by this script (with the number of coding isoform for the gene) :

```sh
python3 ./src/make_metadata_and_similarity_table.py -m data/pannzer_output/human.all.nr_off.out -s data/pannzer_output/human.long.nr_off.out -g ~/Software/GOGO/ -f -o human_allVSlong_sim
```

- `-m` : path to a pannzer output as input file (from a multiple-isoform annotation)
- `-s` : path to a pannzer output as input file (from a single-isoform annotation)
- `-g` : path of GOGO directory
- `-a` : infer GO term ancestry (longer)
- `-f` : filter gene_set to have only multiple-isoform gene used for similarity
- `-o` : name of the output file



### Write exhaustive data table with metadata and similarity based on a multiple-isoform annotation

From the Pannzer output of a multiple-isoform annotation's proteome, create two tables with metadata (parsing of the annotation) and optionally similarity table.

```sh
python3 ./src/description_table.py -m data/pannzer_output/human.all.nr_off.out -g ~/Software/GOGO/ -o res/exhaustive -lbc data/pannzer_output/human.mane.nr_off.out
```

- `-m` : path to a pannzer output as input file (from a multiple-isoform annotation)
- `-g` : path of GOGO directory
- `-a` : infer GO term ancestry (longer)
- `-f` : filter gene_set to have only multiple-isoform gene used for similarity
- `-o` : name of the output file
- `-c` : path to a pannzer output as input file (from a custom single-isoform annotation) to compute similarity with this one
- `-l` : to compare similarity with the longest isoform
- `-b` : to compare similarity with the beest isoform (isoform with the highest number of GO terms)

Two files will be generated :

- `output.metadata.txt` from the parsing of the multiple-isoform annotation (one line by **transcript**)
    - `chromosome` : chromosome of the gene
    - `gene_id` : ensembl gene id
    - `nb_isoform` : nb of coding isoform PRESENT in the pannzer output
    - `list_of_go_term_assigned_to_the_gene` : list of go term assigned to the gene by pannzer
    - `nb_go_term_assigned_to_the_gene` : nb of go term assigned to the gene by pannzer
    - `transcript_id` : transcript id of the transcript
    - `description` : predicted description of the transcript
    - `score` : ppv of this prediction
    - `list_of_go_term_assigned_to_the_transcript` : list of go term assigned to the transcript by pannzer
    - `nb_go_term_assigned_to_the_transcript` : nb of go term assigned to the transcript by pannzer


- `output.similarity.txt` if at least one option between `-c -b -l` is present, create a table with similarity (one line by **gene**)
    - `chromosome` : chromosome of the gene
    - `gene_id` : ensembl gene id
    - `nb_isoform` : nb of coding isoform PRESENT in the pannzer output
    - `list_of_go_term_assigned_to_the_gene` : list of go term assigned to the gene by pannzer
    - `nb_go_term_assigned_to_the_gene` : nb of go term assigned to the gene by pannzer

        **only if you use `-c`**
    - `transcript_id_custom` : transcript id of the chosen transcript for the custom single-isoform annotation
    - `description_custom` : predicted description of the chosen transcript
    - `score_custom` : ppv of this prediction
    - `list_of_go_term_assigned_to_the_transcript_custom` : list of go term assigned to the chosen transcript by pannzer
    - `nb_go_term_assigned_to_the_transcript_custom` : nb of go term assigned to the chosen transcript by pannzer
    - `BP_similarity_custom` : GOGO similarity of BP term between the set of go term assigned to the gene and to the custom transcript only
    - `CC_similarity_custom` : GOGO similarity of CC term between the set of go term assigned to the gene and to the custom transcript only
    - `MF_similarity_custom` : GOGO similarity of MF term between the set of go term assigned to the gene and to the custom transcript only

        **only if you use `-l`**

    - `transcript_id_longest` : transcript id of the longest transcript 
    - `description_longest` : predicted description of the longest transcript
    - `score_longest` : ppv of this prediction
    - `list_of_go_term_assigned_to_the_transcript_longest` : list of go term assigned to the chosen transcript by pannzer
    - `nb_go_term_assigned_to_the_transcript_longest` : nb of go term assigned to the chosen transcript by pannzer
    - `BP_similarity_longest` : GOGO similarity of BP term between the set of go term assigned to the gene and to the longest transcript only
    - `CC_similarity_longest` : GOGO similarity of CC term between the set of go term assigned to the gene and to the longest transcript only
    - `MF_similarity_longest` : GOGO similarity of MF term between the set of go term assigned to the gene and to the longest transcript only

        **only if you use `-b`**
 
    - `transcript_id_best` : transcript id of the best transcript (with the highest number of go term) 
    - `description_best` : predicted description of the best transcript
    - `score_best` : ppv of this prediction
    - `list_of_go_term_assigned_to_the_transcript_best` : list of go term assigned to the chosen transcript by pannzer
    - `nb_go_term_assigned_to_the_transcript_best` : nb of go term assigned to the chosen transcript by pannzer
    - `BP_similarity_best` : GOGO similarity of BP term between the set of go term assigned to the gene and to the best transcript only
    - `CC_similarity_best` : GOGO similarity of CC term between the set of go term assigned to the gene and to the best transcript only
    - `MF_similarity_best` : GOGO similarity of MF term between the set of go term assigned to the gene and to the best transcript only

**Note** : Results are not available on GitHub due to size issue.


### Perform GO Enrichment Analysis (GOEA) on previous results

To find if the subset of gene with low similarity is enriched in some terms, you can run :

```sh
python3 ./src/go_enrichment_analysis.py -i res/exhaustive.similarity.txt --bp_column BP_similarity_longest --cc_column CC_similarity_longest \
        --mf_column MF_similarity_longest --bp_thresold 0.9 --ensembl2ncbi data/gene2ensembl.gz --go_dag data/go-basic.obo --gene2go data/gene2go -o res/human_allVSlong
```

- `-i` : TSV table with at least a column for each BP, CC, MF similarity and one for gene ID
- `--bp_column` : Name of the column with BP similarity (default : BP_similarity)
- `--bp_thresold` : Maximum BP similarity to keep a gene list to perform the enrichment (default : 1)
- `--cc_column` : Name of the column with CC similarity (default : CC_similarity)
- `--cc_thresold` : Maximum CC similarity to keep a gene list to perform the enrichment (default : 1)
- `--mf_column` : Name of the column with MF similarity (default : MF_similarity)
- `--mf_thresold` : Maximum MF similarity to keep a gene list to perform the enrichment (default : 1)
- `--ensembl2ncbi` : path of gene2ensembl.gz file from NCBI
- `--go_dag` : path of go-basic.obo file
- `--gene2go` : path of gene2go UNZIP file from NCBI
- `--taxid` : NCBI TaxID of the studied species (default : 9606 for human)
- `-o` : Name of the output file

Four files are generated :

- `output.goea.tsv` : table with enriched go term (and p-value, FDR, etc)
- `output.{NS}_graph_of_significant_GO.png` `x3` : plot of a subgraph of the GO graph with GO term which enrich our gene set.


### Evaluate isoforms diversity for each gene

To evaluate isoforms diversity for each gene with different metrics, you should run :

```sh
python3 src/intragene_isoform_diversity.py -i data/pannzer_output/human.all.nr_off.out -o res/isoforms_diversity
```

#### Metrics computed in all isoform :

Let $G$, the set, a set of isoform defined as :
$$ G = \{ I_1, I_2, \ldots, I_n \} $$
where each $I_i$ is a set of multiple GO terms defined as :
$$ I_i = \{ T_{i1}, T_{i2}, \ldots, T_{im_i} \} $$


***Number of isoform : number of isoform***
$$n_{isoform} = |G|$$

***Standard deviation of the number of GO term***
$$\sigma_{m_i} = \sqrt{\frac{1}{n}\sum_{i=1}^{n}{(m_i-\bar{m_i} )^2}}$$

***Redudancy metric***\
This metrics was designed to have an idea of the number of times a GO term appear in the genes.
For each unique GO term, his number of reoccurence is count ($0$ if it appear only in $1$ isoform, and $n-1$ if it appear in all isoform). Then, the mean of this counting is done to have the mean count of reoccurence. Finally, this count is divided by the $n-1$.
If $r$ is close to $1$, that's means that all GO terms are present in all isoforms, and if $r$ is close to $0$, each isoform is different.

If there is 1 isoform, the redudancy metric is set to $1$.

Let $O$ be the set of all unique GO terms defined as :
$$O = \bigcup_{i=1}^n I_i = \{ T_{1}, T_{2}, \ldots, T_{n_o} \}$$
Let $count(T_i)$ the number of isoform where $T_i$ is present.
$$r = \frac{1}{n-1} \sum_{i=1}^ {n_o}(count(T_i)-1)$$


#### Metrics computed in all isoform :
Certain metrics were calculated for each pair of isoforms before all the values were averaged.

***Jaccard index***\
The Jaccard index is measure of similarity.\
If $I_1 \cup I_2  = \emptyset$, $J(I_1, I_2) = 1$, else

$$J(I_1, I_2) = \frac{|I_1 \cap I_2|}{|I_1 \cup I_2|}$$

***Dice coefficient***\
The Sørensen–Dice coefficient is measure of similarity.\
If $I_1 \cup I_2  = \emptyset$, $D(I_1, I_2) = 1$, else

$$D(I_1, I_2) = \frac{2|I_1 \cap I_2|}{|I_1| + |I_2|}$$

***Overlap coefficient***\
If $min(|I_1|,|I_2|) = 0$, then $overlap(I_1, I_2) = 1$, else

$$overlap(I_1, I_2) = \frac{|I_1 \cap I_2|}{min(|I_1|,|I_2|)}$$


Two files are generated :
- `output.data.tsv` : measurement for each gene
- `output.summary.tsv` : descriptive statistics for each measurement