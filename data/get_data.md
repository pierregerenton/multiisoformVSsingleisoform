# Getting data

This file contains usefull information to regenerate input file for all scripts of this repository if you don't have them, as data are not provided.

## Human and mouse genome and annotation

Gencode assembly and annotation have been chosen. Data are available here :

- Human :
    - [Genome](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.p14.genome.fa.gz)
    - [Gene annotation](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gff3.gz)
- Mouse :
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


