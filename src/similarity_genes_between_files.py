import argparse
from matplotlib import pyplot
import pannzer_out_api as poa
from os.path import basename

parser = argparse.ArgumentParser(
                    prog='similarity',
                    description='Create 3 table for each ontology \
                        with the mean similarity for each pair of \
                            panzzer output in the input + \
                                write raw data',
                    epilog='For more information, contact pierre.gerenton@crg.eu')


parser.add_argument(
    '-i', '--pannzer-output',
    type=str, nargs='+', required=True,
    help='List of Pannzer Output file to compare (by pair)'
)

parser.add_argument(
    '-g', '--gogo-dir',
    type=str, required=True,
    help='Path of the GOGO directory'
)

parser.add_argument(
    '-b', '--add-best',
    action='store_true', default=None,
    help='Will add an element for the comparison that is the best isoform (isoform with the most isoform) from the first file after -i'
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


def check_number_args(args):
    c = len(args.pannzer_output)
    if args.add_best is not None:
        c += 1
    if c < 2:
        print('Need at least 2 element to compare')
        exit()

def read_pannzer_outputs(list_path):
    print("Reading Pannzer outputs ...")
    list_annotation = []
    for path in list_path:
        name = basename(path).strip('.out')
        annotation = poa.parse_pannzer_annotation(path, name)
        list_annotation.append(annotation)
    print("Done")
    return list_annotation



def main():

    check_number_args(args)

    pannzer_outputs = read_pannzer_outputs(args.pannzer_output)
    print()

    if args.add_best:
        print('Adding element with best isoform ...')
        pannzer_outputs.append(poa.make_best_single_isoform_annotation(pannzer_outputs[0]))
        print('Done')
        print()
    
    print('Preparing gene set ...')
    gene_set = list(pannzer_outputs[0].genes)  # getting genes from PO 1
    for annotation in pannzer_outputs[1:]:  # Keeping genes if also present in other PO
        gene_set = [gene for gene in annotation.genes if gene in gene_set]
    if args.only_multiple_isoform==True:
        gene_set = [gene for gene in gene_set if len(pannzer_outputs[0].genes[gene].transcripts) > 1]
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

    bp_table = []
    for i in range(len(pannzer_outputs)):
        line = []
        for j in range(len(pannzer_outputs)):
            line.append(' ')
        bp_table.append(line)

    cc_table = []
    for i in range(len(pannzer_outputs)):
        line = []
        for j in range(len(pannzer_outputs)):
            line.append(' ')
        cc_table.append(line)

    mf_table = []
    for i in range(len(pannzer_outputs)):
        line = []
        for j in range(len(pannzer_outputs)):
            line.append(' ')
        mf_table.append(line)

    with open(args.output_name + '.txt','w') as txt_output:
        txt_output.write(f"COMPARISON\tGENE_ID\tBP_SIMILARITY\tCC_SIMILARITY\tMF_SIMILARITY\n")

    for i in range(len(pannzer_outputs)-1):
        annotation_1 = pannzer_outputs[i]
        bp_table[i][i] = 1
        cc_table[i][i] = 1
        mf_table[i][i] = 1
        for j in range(i+1,len(pannzer_outputs)):
            annotation_2 = pannzer_outputs[j]
            comb = annotation_1.name + '_X_' + annotation_2.name
            print(f"Working on {comb}")
            similarity = fun(annotation_1, annotation_2, args.gogo_dir, **argument)
            with open(args.output_name + '.txt','a') as txt_output:
                for gene in similarity:
                    txt_output.write(f"{comb}\t{gene}\t{similarity[gene]['BP']}\t{similarity[gene]['CC']}\t{similarity[gene]['MF']}\n")
                msim = poa.mean_similarity(similarity)
                bp_table[i][j]=msim['BP']
                cc_table[i][j]=msim['CC']
                mf_table[i][j]=msim['MF']
    bp_table[j][j] = 1
    cc_table[j][j] = 1
    mf_table[j][j] = 1

    print('Done !')
    print()

    print('Ploting ...')

    rhead = chead = [annot.name for annot in pannzer_outputs]
    col_table = []
    for i in range(len(pannzer_outputs)):
        line = []
        for j in range(len(pannzer_outputs)):
            if j>i:
                line.append('lightgreen')
            else:
                line.append('lightgray')
        col_table.append(line)

    fig, axs = pyplot.subplots(3, constrained_layout=True)

    ax = axs[0]
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.axis('off')
    ax.set_title('BP Similarity')
    the_table = ax.table(bp_table, rowLabels=rhead, colLabels=chead, loc='center', cellColours=col_table)
    the_table.scale(1, 1.5)

    ax = axs[1]
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.axis('off')
    ax.set_title('CC Similarity')
    the_table = ax.table(cc_table, rowLabels=rhead, colLabels=chead, loc='center', cellColours=col_table)
    the_table.scale(1, 1.5)

    ax = axs[2]
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.axis('off')
    ax.set_title('MF Similarity')
    the_table = ax.table(mf_table, rowLabels=rhead, colLabels=chead, loc='center', cellColours=col_table)
    the_table.scale(1, 1.5)

    pyplot.suptitle('Mean gene GOGO Similarity for different ontology between files')
    pyplot.draw()
    pyplot.savefig(args.output_name + '.pdf', format='pdf')
    print('Done !')



if __name__=="__main__":
    main()


