import pannzer_out_api as poa
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import seaborn
import argparse
from os.path import basename



parser = argparse.ArgumentParser(
                    prog='precise_analysis_of_one_multiisoform_annotation',
                    description='Print a lot of plot when comparing a reference annotation\
                        with an alternative single-isoform annotation',
                    epilog='For more information, contact pierre.gerenton@crg.eu')


parser.add_argument(
    '-i', '--pannzer-output',
    type=str, required=True,
    help='Reference multiple isoform annotation'
)

parser.add_argument(
    '-g', '--gogo-dir',
    type=str, required=True,
    help='Path of the GOGO directory'
)

parser.add_argument(
    '-t', '--type',
    choices=['long','best'], default=None,
    help='Compare reference to longest OR best isoform'
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

    print('Reading input ...')
    reference_annotation = poa.parse_pannzer_annotation(args.pannzer_output)
    shuffle_annotation = poa.make_shuffle_copy(reference_annotation)
    
    if args.type == 'long':
        alternative_annotation = poa.make_longest_single_isoform_annotation(reference_annotation)
        alternative_shuffle_annotation = poa.make_longest_single_isoform_annotation(shuffle_annotation)
    elif args.type == 'best':
        alternative_annotation = poa.make_best_single_isoform_annotation(reference_annotation)
        alternative_shuffle_annotation = poa.make_best_single_isoform_annotation(shuffle_annotation)
    else:
        raise ValueError('unknown single isoform type')
    


    print('Done\n')



    print('Preparing gene set ...')
    gene_set = list(reference_annotation.genes)
    if args.only_multiple_isoform==True:
        gene_set = [gene for gene in gene_set if len(reference_annotation.genes[gene].transcripts) > 1]

    print('Done')
    print()


    print('Computing similarity ...')
    similarity_observed = poa.gogo_similarity_between_annotation(reference_annotation, alternative_annotation, args.gogo_dir, gene_set)
    BP_observed = [ float(similarity_observed[gene]['BP']) for gene in similarity_observed if similarity_observed[gene]['BP'] != 'NA']
    CC_observed = [ float(similarity_observed[gene]['CC']) for gene in similarity_observed if similarity_observed[gene]['CC'] != 'NA']
    MF_observed = [ float(similarity_observed[gene]['MF']) for gene in similarity_observed if similarity_observed[gene]['MF'] != 'NA']
    print('Mean observed similarity :', poa.mean_similarity(similarity_observed))
    similarity_expected = poa.gogo_similarity_between_annotation(shuffle_annotation, alternative_shuffle_annotation, args.gogo_dir, gene_set)
    BP_expected = [ float(similarity_expected[gene]['BP']) for gene in similarity_expected if similarity_expected[gene]['BP'] != 'NA']
    CC_expected = [ float(similarity_expected[gene]['CC']) for gene in similarity_expected if similarity_expected[gene]['CC'] != 'NA']
    MF_expected = [ float(similarity_expected[gene]['MF']) for gene in similarity_expected if similarity_expected[gene]['MF'] != 'NA'] 
    print('Mean expected similarity :', poa.mean_similarity(similarity_expected))
    print('Done')
    print()    

    print('Ploting ...')

    bins_number = int(len(gene_set)/100)
    figs = []

    # Observed distribution 
    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.histplot(BP_observed, bins=bins_number,element='poly')
    seaborn.histplot(CC_observed, bins=bins_number,element="poly")
    seaborn.histplot(MF_observed, bins=bins_number,element="poly")
    pyplot.legend(['BP', 'CC', 'MF'])
    pyplot.title("Observed distribution of similarity between coding genes")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of genes")
    pyplot.ylim(0,250)
    pyplot.draw()

    # Expected distribution 
    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.histplot(BP_expected, bins=bins_number,element='poly')
    seaborn.histplot(CC_expected, bins=bins_number,element="poly")
    seaborn.histplot(MF_expected, bins=bins_number,element="poly")
    pyplot.legend(['BP', 'CC', 'MF'])
    pyplot.title("Expected distribution of similarity between coding genes")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of genes")
    pyplot.ylim(0,250)
    pyplot.draw()

    # BP
    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.histplot(BP_expected, bins=bins_number, element="poly")
    seaborn.histplot(BP_observed, bins=bins_number, element='poly')
    pyplot.legend(['BP : Expected Distribution', 'BP : Observed Distribution'])
    pyplot.title("Observed and expected distribution of the BP similarity between a gene and its longest CDS")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of genes")
    pyplot.ylim(0,250)


    x_similarity = []
    y_nb_isoform = []
    for gene in gene_set:
        if similarity_expected[gene]['BP'] != 'NA':
            x_similarity.append(float(similarity_expected[gene]['BP']))
            y_nb_isoform.append(len(shuffle_annotation.genes[gene].transcripts))

    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.scatterplot(x = x_similarity, y = y_nb_isoform)
    pyplot.title("Number of isoform by BP similarity of genes (expected)")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of isoforms")

    x_similarity = []
    y_nb_isoform = []
    for gene in gene_set:
        if similarity_observed[gene]['BP'] != 'NA':
            x_similarity.append(float(similarity_observed[gene]['BP']))
            y_nb_isoform.append(len(reference_annotation.genes[gene].transcripts))

    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.scatterplot(x = x_similarity, y = y_nb_isoform)
    pyplot.title("Number of isoform by BP similarity of genes (observed)")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of isoforms")

    # CC
    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.histplot(CC_expected, bins=bins_number, element="poly")
    seaborn.histplot(CC_observed, bins=bins_number, element='poly')
    pyplot.legend(['CC : Expected Distribution', 'CC : Observed Distribution'])
    pyplot.title("Observed and expected distribution of the CC similarity between a gene and its longest CDS")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of genes")
    pyplot.ylim(0,250)


    x_similarity = []
    y_nb_isoform = []
    for gene in gene_set:
        if similarity_expected[gene]['CC'] != 'NA':
            x_similarity.append(float(similarity_expected[gene]['CC']))
            y_nb_isoform.append(len(shuffle_annotation.genes[gene].transcripts))

    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.scatterplot(x = x_similarity, y = y_nb_isoform)
    pyplot.title("Number of isoform by CC similarity of genes (expected)")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of isoforms")

    x_similarity = []
    y_nb_isoform = []
    for gene in gene_set:
        if similarity_observed[gene]['CC'] != 'NA':
            x_similarity.append(float(similarity_observed[gene]['CC']))
            y_nb_isoform.append(len(reference_annotation.genes[gene].transcripts))

    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.scatterplot(x = x_similarity, y = y_nb_isoform)
    pyplot.title("Number of isoform by CC similarity of genes (observed)")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of isoforms")

    # MF

    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.histplot(MF_expected, bins=bins_number, element="poly")
    seaborn.histplot(MF_observed, bins=bins_number, element='poly')
    pyplot.legend(['MF : Expected Distribution', 'MF : Observed Distribution'])
    pyplot.title("Observed and expected distribution of the MF similarity between a gene and its longest CDS")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of genes")
    pyplot.ylim(0,250)


    x_similarity = []
    y_nb_isoform = []
    for gene in gene_set:
        if similarity_expected[gene]['MF'] != 'NA':
            x_similarity.append(float(similarity_expected[gene]['MF']))
            y_nb_isoform.append(len(shuffle_annotation.genes[gene].transcripts))

    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.scatterplot(x = x_similarity, y = y_nb_isoform)
    pyplot.title("Number of isoform by MF similarity of genes (expected)")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of isoforms")

    x_similarity = []
    y_nb_isoform = []
    for gene in gene_set:
        if similarity_observed[gene]['MF'] != 'NA':
            x_similarity.append(float(similarity_observed[gene]['MF']))
            y_nb_isoform.append(len(reference_annotation.genes[gene].transcripts))

    figs.append(pyplot.figure(figsize=(15, 6)))
    seaborn.scatterplot(x = x_similarity, y = y_nb_isoform)
    pyplot.title("Number of isoform by MF similarity of genes (observed)")
    pyplot.xlabel("GOGO Similarity")
    pyplot.ylabel("Number of isoforms")


    with PdfPages(args.output_name) as pdf:
        for fig in figs:
            pdf.savefig(fig, bbox_inches='tight') 

    print('Done')
    print()


if __name__=="__main__":
    main()


