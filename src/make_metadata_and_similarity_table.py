import argparse
import pannzer_out_api as poa


parser = argparse.ArgumentParser(
                    prog='table_sim_metadata',
                    description='Write a table with BP, CC, MF similarity\
                        for each gene and the number of isoform',
                    epilog='For more information, contact pierre.gerenton@crg.eu')


parser.add_argument(
    '-m', '--multi-isoform-pannzer-output',
    type=str, required=True,
    help='Path of a panzzer output (from a multiple isoform annotation)'
)

parser.add_argument(
    '-s', '--single-isoform-pannzer-output',
    type=str, required=True,
    help='Path of a panzzer output (from a single isoform annotation)'
)

parser.add_argument(
    '-g', '--gogo-dir',
    type=str, required=True,
    help='Path of the GOGO directory'
)

parser.add_argument(
    '-a', '--ancestry',
    action='store_true',
    help='Will add all parents GO term of all GO term of a gene'
)

parser.add_argument(
    '-f', '--only-multiple-isoform',
    action='store_true', default=False,
    help='Compute similarity only for gene with multiple coding isoform (based on the first file after -i)'
)

parser.add_argument(
    '-o', '--output-name',
    type=str, default='number_genes_with_different_go_term_between_files',
    help='Name of the output file'
)


args = parser.parse_args()





def main():

    multi_isoform = poa.parse_pannzer_annotation(args.multi_isoform_pannzer_output, 'multi')
    single_isoform = poa.parse_pannzer_annotation(args.single_isoform_pannzer_output, 'single')


    print('Preparing gene set ...')
    gene_set = list(multi_isoform.genes)  # getting genes from multi isoform
    gene_set = [gene for gene in single_isoform.genes if gene in gene_set]  # Keeping genes if also present in singlecisoform
    if args.only_multiple_isoform==True:
        gene_set = [gene for gene in gene_set if len(multi_isoform.genes[gene].transcripts) > 1]
    print('Done')
    print()

    argument = {'gene_set':gene_set}
    fun = poa.gogo_similarity_between_annotation

    if args.ancestry: 
        print('Reading ancestry ...')
        go_graph = poa.read_go_graph()
        argument['go_graph'] = go_graph
        fun = poa.gogo_similarity_between_annotation_with_hierarchy
        print('Done')
        print()
    
    print('Computing values ...')

    with open(args.output_name + '.txt','w') as txt_output:
        txt_output.write(f"GENE_ID\tBP_SIMILARITY\tCC_SIMILARITY\tMF_SIMILARITY\tNB_ISOFORM\n")

    similarity = fun(multi_isoform, single_isoform, args.gogo_dir, **argument)
    print(poa.mean_similarity(similarity))

    for gene in similarity:
        similarity[gene]['nb_isoform'] = len(multi_isoform[gene].transcripts)
    
 
    print('Done !')
    print()

    print('Writing table ...')
      
    with open(args.output_name + '.txt','a') as txt_output:
        for gene in similarity:
            txt_output.write(f"{gene}\t{similarity[gene]['BP']}\t{similarity[gene]['CC']}\t{similarity[gene]['MF']}\t{similarity[gene]['nb_isoform']}\n")

    print('Done !')


if __name__=="__main__":
    main()


