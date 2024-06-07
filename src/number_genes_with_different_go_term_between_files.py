import argparse
import upsetplot
from matplotlib import pyplot
import pannzer_out_api as poa
from os.path import basename

__author__ = "Gérenton Pierre"
__credits__ = ["Gérenton Pierre", "Fabio Zanarello", "Roderic Guigó i Serra"]
__license__ = "CC0 1.0 Universal"

parser = argparse.ArgumentParser(
                    prog='number_genes_with_different_go_term_between_files',
                    description='For each pair of annotation (in input), will\
                          compute the number of genes where GO annotation \
                            differ and print an UpSetPlot',
                    epilog='For more information, contact fabio.zanarello@crg.eu')


parser.add_argument(
    '-i', '--pannzer-output',
    type=str, nargs='+', required=True,
    help='List of Pannzer Output file to compare (by pair)'
)

parser.add_argument(
    '-l', '--add-longest',
    type=str, metavar='PANZZER_OUTPUT', default=None,
    help='Will add an element for the comparison that is the longest isoform (isoform with the most isoform) from the file after -b'
)

parser.add_argument(
    '-b', '--add-best',
    type=str, metavar='PANZZER_OUTPUT', default=None,
    help='Will add an element for the comparison that is the best isoform (isoform with the most isoform) from the file after -b'
)

parser.add_argument(
    '-a', '--ancestry',
    action='store_true',
    help='Will add all parents GO term of all GO term of a gene'
)

parser.add_argument(
    '-o', '--output-name',
    type=str, default='number_genes_with_different_go_term_between_files',
    help='Name of the output file'
)



args = parser.parse_args()


def check_number_args(args):
    c = len(args.pannzer_output)
    if args.add_longest is not None:
        c += 1    
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

def read_longest(path):
    print("Adding longest isoform annotation ...")
    name = basename(path).strip('.out')
    annotation = poa.parse_pannzer_annotation(path, name)
    longest_annotation = poa.make_longest_single_isoform_annotation(annotation)
    print("Done")
    return longest_annotation

def read_best(path):
    print("Adding best isoform annotation ...")
    name = basename(path).strip('.out')
    annotation = poa.parse_pannzer_annotation(path, name)
    best_annotation = poa.make_best_single_isoform_annotation(annotation)
    print("Done")
    return best_annotation


def main():

    check_number_args(args)

    pannzer_outputs = read_pannzer_outputs(args.pannzer_output)
    print()

    if args.add_longest is not None:
        pannzer_outputs.append(read_longest(args.add_longest))
        print()

    if args.add_best is not None:
        pannzer_outputs.append(read_best(args.add_best))
        print()
    
    print('Preparing gene set ...')
    gene_set = list(pannzer_outputs[0].genes)  # getting genes from PO 1
    for annotation in pannzer_outputs[1:]:  # Keeping genes if also present in other PO
        gene_set = [gene for gene in annotation.genes if gene in gene_set]
    print('Done')
    print()

    argument = {'gene_set':gene_set}
    fun = poa.genes_with_diff_go_terms

    if args.ancestry: 
        print('Reading ancestry ...')
        go_graph = poa.read_go_graph()
        argument['go_graph'] = go_graph
        fun = poa.genes_with_diff_go_terms_with_hierarchy
        print('Done')
        print()
    
    print('Computing values ...')
    comb_name = []
    comb_val = []
    for i in range(len(pannzer_outputs)-1):
        annotation_1 = pannzer_outputs[i]
        for annotation_2 in pannzer_outputs[i+1:]:
            comb_name.append([annotation_1.name, annotation_2.name])
            comb_val.append(len(fun(annotation_1, annotation_2, **argument)))
    print('Done !')
    print()

    print('Ploting ...')
    data = upsetplot.from_memberships(comb_name, data = comb_val)
    upsetplot.plot(data, sort_by='cardinality', sort_categories_by='input', totals_plot_elements=0 )
    pyplot.title(f'Number of genes with \na different GO annotation\nbetween files\nTotal number of genes : {len(gene_set)}')
    pyplot.ylabel('Number of genes')
    pyplot.savefig(args.output_name, format='pdf',bbox_inches="tight")
    print('Done !')



if __name__=="__main__":
    main()


