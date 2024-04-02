## Old scripts

Some scripts were written / used but now aren't useful (too buggy, not use at all anymore, etc).
However, we decided to keep them because some old files were generated with them.


### make_metadata_and_similarity_table.py

Removed because it is redundant with description_table.py .

__Usage :__

With two Pannzer output, one from a multiple-isoform annotation and the other from a single-isoform annotation, a table with the similarity for each gene between files is written by this script (with the number of coding isoform for the gene) :

```sh
python3 ./src/make_metadata_and_similarity_table.py -m data/pannzer_output/human.all.nr_off.out -s data/pannzer_output/human.long.nr_off.out -g ~/Software/GOGO/ -f -o human_allVSlong_sim
```

- `-m` : path to a pannzer output as input file (from a multiple-isoform annotation) \[MANDATORY\]
- `-s` : path to a pannzer output as input file (from a single-isoform annotation) \[MANDATORY\]
- `-g` : path of GOGO directory \[MANDATORY\]
- `-a` : infer GO term ancestry (longer)
- `-f` : filter gene_set to have only multiple-isoform gene used for similarity
- `-o` : name of the output file \[MANDATORY\]
