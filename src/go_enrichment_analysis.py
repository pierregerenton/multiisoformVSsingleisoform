import argparse
import pandas as pd  # use data frame
import collections as cx  # use Counter
from goatools.obo_parser import GODag  # parse the .obo file (with GO ontology)
from goatools.anno.genetogo_reader import Gene2GoReader  # parse NCBI gene2go file
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS  # perform enrichment
from goatools.godag_plot import plot_results  # ploting


parser = argparse.ArgumentParser(
                    prog='go_enrichment_analysis',
                    description='Write a table with the result of the GOEA \
                        and some plot of the GO DAG with those enriched segment.',
                    epilog='For more information, contact pierre.gerenton@crg.eu')


parser.add_argument(
    '-i', '--input',
    type=str, required=True,
    help='TSV table with at least a column for each BP, CC, MF similarity and one for gene ID'
)

parser.add_argument(
    '--bp_column',
    type=str, default='BP_similarity', required=False,
    help='Name of the column with BP similarity'
)

parser.add_argument(
    '--bp_thresold',
    type=float, default=None, required=False,
    help='Maximum BP similarity to keep a gene list to perform the enrichment'
)

parser.add_argument(
    '--cc_column',
    type=str, default='CC_similarity', required=False,
    help='Name of the column with CC similarity'
)

parser.add_argument(
    '--cc_thresold',
    type=float, default=None, required=False,
    help='Maximum CC similarity to keep a gene list to perform the enrichment'
)

parser.add_argument(
    '--mf_column',
    type=str, default='MF_similarity', required=False,
    help='Name of the column with MF similarity'
)

parser.add_argument(
    '--mf_thresold',
    type=float, default=None, required=False,
    help='Maximum MF similarity to keep a gene list to perform the enrichment'
)

parser.add_argument(
    '--ensembl2ncbi',
    type=str, required=True,
    help="path of gene2ensembl.gz file from NCBI"
)

parser.add_argument(
    '--go_dag',
    type=str, required=True,
    help="path of go-basic.obo file"
)

parser.add_argument(
    '--gene2go',
    type=str, required=True,
    help="path of gene2go UNZIP file from NCBI"
)

parser.add_argument(
    '--taxid',
    type=int, required=False, default=9606,  #human
    help="NCBI TaxID of the studied species"
)

parser.add_argument(
    '-o', '--output-name',
    type=str, default='number_genes_with_different_go_term_between_files',
    help='Name of the output file'
)


args = parser.parse_args()

def main():
    
    print("Preparing background and gene list ...")
    
    similarity_table = pd.read_csv(args.input, sep='\t')
    ensembl2ncbi_df = pd.read_csv(args.ensembl2ncbi, sep='\t', compression='gzip')

    ensembl2ncbi_dict = dict(zip(ensembl2ncbi_df['Ensembl_gene_identifier'], ensembl2ncbi_df['GeneID']))

    background = similarity_table['gene_id'].tolist()
    nb_gene_available = len(background)
    background = [ensembl2ncbi_dict[ensembl_id] for ensembl_id in background if ensembl_id in ensembl2ncbi_dict]
    nb_gene_found = len(background)
    print(f"{nb_gene_found}/{nb_gene_available} used as background")

    if args.bp_thresold is not None:
        similarity_table = similarity_table[similarity_table[args.bp_column] <= args.bp_thresold]
    if args.cc_thresold is not None:
        similarity_table = similarity_table[similarity_table[args.cc_column] <= args.cc_thresold]
    if args.mf_thresold is not None:
        similarity_table = similarity_table[similarity_table[args.mf_column] <= args.mf_thresold]
    gene_list = similarity_table['gene_id'].tolist()
    nb_gene_available = len(gene_list)
    gene_list = [ensembl2ncbi_dict[ensembl_id] for ensembl_id in gene_list if ensembl_id in ensembl2ncbi_dict]
    nb_gene_found = len(gene_list)
    print(f"{nb_gene_found}/{nb_gene_available} used as gene_list")

    print("Done")
    print()

    print("Preparation of GOATOOLS ...")

    obo_dag = GODag(args.go_dag)  # no obsolete
    gene2go = Gene2GoReader(args.gene2go, taxids=[args.taxid])
    ns2assoc = gene2go.get_ns2assc()
    goea_obj = GOEnrichmentStudyNS(
        pop = background,
        ns2assoc = ns2assoc,
        godag = obo_dag,
        propagate_counts = False,
        alpha = 0.05,
        methods = ['fdr_bh']
    )
    print('Done !')
    print()

    print('Run analysis ...')

    goea_results = goea_obj.run_study(gene_list)
    goea_results_sig = [r for r in goea_results if r.p_fdr_bh < 0.05]

    print('Done !')
    print()

    print('Plot and write results ...')

    print('{N} of {M:,} results were significant'.format(
        N=len(goea_results_sig),
        M=len(goea_results)))
    print('Significant results: {E} enriched, {P} purified'.format(
        E=sum(1 for r in goea_results_sig if r.enrichment=='e'),
        P=sum(1 for r in goea_results_sig if r.enrichment=='p')))   
    ctr = cx.Counter([r.NS for r in goea_results_sig])
    print('Significant results[{TOTAL}] = {BP} BP + {MF} MF + {CC} CC'.format(
        TOTAL=len(goea_results_sig),
        BP=ctr['BP'],  # biological_process
        MF=ctr['MF'],  # molecular_function
        CC=ctr['CC'])) # cellular_component
    goea_obj.wr_tsv(args.output_name + ".goea_results.tsv", goea_results_sig)
    plot_results(args.output_name + ".{NS}_graph_of_significative_GO.png", goea_results_sig)


if __name__=="__main__":
    main()


