import pannzer_out_api as poa
import statistics as stats
import pandas as pd
import argparse
from os.path import basename


parser = argparse.ArgumentParser(
                    prog='intragene_isoform_diversity',
                    description='Compute diversity between all isoform of a gene',
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
    '-o', '--output-name',
    type=str, default='number_genes_with_different_go_term_between_files',
    help='Name of the output file'
)


args = parser.parse_args()


def precise_quantile(sample, n : int) -> float:
    return stats.quantiles(sample, n = 100)[n-1]



def main():

    print('Reading input ...')

    annotation : poa.Annotation = poa.parse_pannzer_annotation(args.pannzer_output)
    shuffle : poa.Annotation = poa.make_shuffle_copy(annotation)

    print('Done\n')

    print('Computing diversity for given data ...')

    data = pd.DataFrame()
    data['Gene'] = annotation.genes.keys()
    data['Type'] = "Observed"
    data['Number of isoform'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.number_of_isoforms)
    data['Jaccard Index'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.diversity_by_pair, similarity_function=poa.jaccard_index)
    data['Dice coefficient'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.diversity_by_pair, similarity_function=poa.dice_coefficient)
    data['Overlap coefficient'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.diversity_by_pair, similarity_function=poa.overlap_coefficient)
    data['Redundance metric'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.go_redundance_metric)
    data['Stdev number GO term'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.stdev_number_of_go_by_isoform)
    similarity = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.intra_gene_gogo_similarity, gogo_dir=args.gogo_dir, gogo_file='iig.observed.'+basename(args.pannzer_output).strip('.out'))
    data['BP similarity'] = similarity.apply(lambda x : x[0])
    data['CC similarity'] = similarity.apply(lambda x : x[1])
    data['MF similarity'] = similarity.apply(lambda x : x[2])

    print('Done\n')

    print('Computing diversity for shuffle data ...')

    reassignment = pd.DataFrame()
    reassignment['Gene'] = shuffle.genes.keys()
    reassignment['Type'] = "Expected"
    reassignment['Number of isoform'] = reassignment['Gene'].apply(shuffle.get_gene).apply(poa.Gene.number_of_isoforms)
    reassignment['Jaccard Index'] = reassignment['Gene'].apply(shuffle.get_gene).apply(poa.Gene.diversity_by_pair, similarity_function=poa.jaccard_index)
    reassignment['Dice coefficient'] = reassignment['Gene'].apply(shuffle.get_gene).apply(poa.Gene.diversity_by_pair, similarity_function=poa.dice_coefficient)
    reassignment['Overlap coefficient'] = reassignment['Gene'].apply(shuffle.get_gene).apply(poa.Gene.diversity_by_pair, similarity_function=poa.overlap_coefficient)
    reassignment['Redundance metric'] = reassignment['Gene'].apply(shuffle.get_gene).apply(poa.Gene.go_redundance_metric)
    reassignment['Stdev number GO term'] = reassignment['Gene'].apply(shuffle.get_gene).apply(poa.Gene.stdev_number_of_go_by_isoform)
    similarity = reassignment['Gene'].apply(shuffle.get_gene).apply(poa.Gene.intra_gene_gogo_similarity, gogo_dir=args.gogo_dir, gogo_file='iig.expected.'+basename(args.pannzer_output).strip('.out'))
    reassignment['BP similarity'] = similarity.apply(lambda x : x[0])
    reassignment['CC similarity'] = similarity.apply(lambda x : x[1])
    reassignment['MF similarity'] = similarity.apply(lambda x : x[2])

    print('Done\n')

    summary = pd.DataFrame()
    summary['Metrics'] = data.columns[2:]
    summary['Mean (observed)'] = summary['Metrics'].apply(data.get).apply(stats.fmean, axis = 1)
    summary['Mean (expected)'] = summary['Metrics'].apply(reassignment.get).apply(stats.fmean, axis = 1)
    summary['Harmonic Mean (observed)'] = summary['Metrics'].apply(data.get).apply(stats.harmonic_mean, axis = 1)
    summary['Harmonic Mean (expected)'] = summary['Metrics'].apply(reassignment.get).apply(stats.harmonic_mean, axis = 1)
    summary['Median (observed)'] = summary['Metrics'].apply(data.get).apply(stats.median, axis = 1)
    summary['Median (expected)'] = summary['Metrics'].apply(reassignment.get).apply(stats.median, axis = 1)
    summary['Q25 (observed)'] = summary['Metrics'].apply(data.get).apply(precise_quantile, n=25, axis = 1)
    summary['Q25 (expected)'] = summary['Metrics'].apply(reassignment.get).apply(precise_quantile, n=25, axis = 1)
    summary['Q75 (observed)'] = summary['Metrics'].apply(data.get).apply(precise_quantile, n=75, axis = 1)
    summary['Q75 (expected)'] = summary['Metrics'].apply(reassignment.get).apply(precise_quantile, n=75, axis = 1)
    summary['Sample Variance (observed)'] = summary['Metrics'].apply(data.get).apply(stats.variance, axis = 1)
    summary['Sample Variance (expected)'] = summary['Metrics'].apply(reassignment.get).apply(stats.variance, axis = 1)
    summary['Sample Standard Deviation (observed)'] = summary['Metrics'].apply(data.get).apply(stats.stdev, axis = 1)
    summary['Sample Standard Deviation (expected)'] = summary['Metrics'].apply(reassignment.get).apply(stats.stdev, axis = 1)
    summary['Population Variance (observed)'] = summary['Metrics'].apply(data.get).apply(stats.pvariance, axis = 1)
    summary['Population Variance (expected)'] = summary['Metrics'].apply(reassignment.get).apply(stats.pvariance, axis = 1)
    summary['Population Standard Deviation (observed)'] = summary['Metrics'].apply(data.get).apply(stats.pstdev, axis = 1)
    summary['Population Standard Deviation (expected)'] = summary['Metrics'].apply(reassignment.get).apply(stats.pstdev, axis = 1)



    print('Done')
    print()

    print('Writing output ...')

    full_data = pd.concat((data, reassignment))

    full_data.to_csv(args.output_name + '.data.tsv', sep='\t', index=False)
    summary.to_csv(args.output_name + '.summary.tsv', sep='\t', index=False)

    print('Done')


if __name__=="__main__":
    main()


