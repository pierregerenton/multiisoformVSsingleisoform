import pannzer_out_api as poa
import statistics as stats
import pandas as pd
import argparse


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

    print('Done\n')

    print('Computing diversity ...')

    data = pd.DataFrame()
    data['Gene'] = annotation.genes.keys()
    data['Number of isoform'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.number_of_isoforms)
    data['Jaccard Index'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.diversity_by_pair, similarity_function=poa.jaccard_index)
    data['Dice coefficient'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.diversity_by_pair, similarity_function=poa.dice_coefficient)
    data['Overlap coefficient'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.diversity_by_pair, similarity_function=poa.overlap_coefficient)
    data['Redundance metric'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.go_redundance_metric)
    data['Stdev number GO term'] = data['Gene'].apply(annotation.get_gene).apply(poa.Gene.stdev_number_of_go_by_isoform)

    summary = pd.DataFrame()
    summary['Metrics'] = data.columns[1:]
    summary['Mean'] = summary['Metrics'].apply(data.get).apply(stats.fmean, axis = 1)
    summary['Harmonic Mean'] = summary['Metrics'].apply(data.get).apply(stats.harmonic_mean, axis = 1)
    summary['Median'] = summary['Metrics'].apply(data.get).apply(stats.median, axis = 1)
    summary['Q25'] = summary['Metrics'].apply(data.get).apply(precise_quantile, n=25, axis = 1)
    summary['Q75'] = summary['Metrics'].apply(data.get).apply(precise_quantile, n=75, axis = 1)
    summary['Sample Variance'] = summary['Metrics'].apply(data.get).apply(stats.variance, axis = 1)
    summary['Sample Standard Deviation'] = summary['Metrics'].apply(data.get).apply(stats.stdev, axis = 1)
    summary['Population Variance'] = summary['Metrics'].apply(data.get).apply(stats.pvariance, axis = 1)
    summary['Population Standard Deviation'] = summary['Metrics'].apply(data.get).apply(stats.pstdev, axis = 1)

    print('Done')
    print()

    print('Writing output ...')

    data.to_csv(args.output_name + '.data.tsv', sep='\t', index=False)
    summary.to_csv(args.output_name + '.summary.tsv', sep='\t', index=False)

    print('Done')


if __name__=="__main__":
    main()


