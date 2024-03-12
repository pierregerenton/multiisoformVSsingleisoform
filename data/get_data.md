# Getting data

This file contains usefull information to regenerate input file for all scripts of this repository if you don't have them, as data are not provided.

## Human and mouse genome and annotation

Gencode assembly and annotation have been chosen. Data are available [here](https://ftp.ebi.ac.uk/pub/databases/gencode/) :

- Human  (release 45):
    - [Genome](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.p14.genome.fa.gz)
    - [Gene annotation](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gff3.gz)
- Mouse (release M34):
    - [Genome](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/GRCm39.genome.fa.gz) 
    - [Gene annotation](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.annotation.gff3.gz)

To get faster preliminary result, you can explore information of the chromosome 21 of the human genome (the smallest).

*Extracting the chr21 from the Gencode assembly*

There are multiple ways of doing that. An example with `samtools` :
```
samtools faidx GRCh38.p14.genome.fa.gz  # create the index
samtools faidx GRCh38.p14.genome.fa.gz chr21  # get the sequence
```

*Extracting the chr21 annotation*

An example with common `bash` command :
```
zcat gencode.v45.annotation.gff3.gz | grep -E '^chr21'
```

## Special annotation

To compare annotation, an another annotation was created by keeping the longest isoform of each gene. It was done with `AGAT` :

```
agat_sp_keep_longest_isoform.pl -gff [human/mouse].gff3 -o [human/mouse].longest_isoform.gff3
```

Also, the MANE annotation is use for the human genome. You can download it in the [NCBI ftp site](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/) or by clicking [here (release 1.2)](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.ensembl_genomic.gff.gz).

## Extracting proteome

From each annotation, a proteome can be extract. The proteome is a multi-fasta file where each sequence is extracted from the reference genome according to the coordinates of the CDS of a coding transcript in an annotation.

With `AGAT`, you can do that by running :
```
agat_sp_extract_sequences.pl -g [human/mouse].[all/long/mane].gff3 -f [human/mouse].fa -t cds -p -o protein.[human/mouse].[all/long/mane].fa
```

## GO term annotation

GO term annotation was done with the PANNZER web interface at http://ekhidna2.biocenter.helsinki.fi/sanspanz/.
We suggest you to uncheck **Remove redundant filter** (cf [`../src/choice_go_set`](../src/choice_go_set)), precise the species of your dataset (if known) and let other parameters by default.










----------------------------

Toronen P, Medlar A, Holm L (2018) PANNZER2: A rapid functional annotation webserver. Nucl. Acids Res. 46, W84-W88
