import pandas as pd
import argparse


parser : argparse.ArgumentParser = argparse.ArgumentParser(
                    prog='intragene_isoform_diversity',
                    description='Compute diversity between all isoform of a gene',
                    epilog='For more information, contact pierre.gerenton@crg.eu')


parser.add_argument(
    '-i', '--input',
    type=str, required=True,
    help='Input GFF'
)


parser.add_argument(
    '-o', '--output',
    type=str, required=True,
    help='Output GFF'
)

args : argparse.Namespace = parser.parse_args()

def infogff2dict(x) -> dict[str, str] :
    return dict( (k,v) for k,v in [element.split('=') for element in x.split(';')])

def keeprow(row, selected_gene : set[str], selected_transcript : set[str]):
    if row[2] == 'gene':
        return infogff2dict(row[8])['gene_id'] in selected_gene
    else:
        return infogff2dict(row[8])['transcript_id'] in selected_transcript

def main():

    # Read input 

    gff : pd.DataFrame = pd.read_csv(args.input , sep='\t', header=None, comment='#')

    # Removing useless features

    gff = gff[gff[2].isin(['gene', 'transcript', 'exon', 'CDS'])]


    # Working with subdataframe for gene

    gene : pd.DataFrame = gff[gff[2].isin(['gene'])]

    only_protein_coding = gene[8].apply(lambda x : infogff2dict(x)["gene_type"] == "protein_coding")
    gene = gene[only_protein_coding]
    
    with_readthrough = gene[8].apply(lambda x : "readthrough" in x)
    gene = gene[~with_readthrough]
    
    gene['ID'] = gene[8].apply(lambda x : x.split(';')[0].lstrip('ID='))
    selected_gene : set = set(gene['ID'].tolist())

    print(f"Number of non-readthrough coding gene : {len(selected_gene)}")

    # Working with subdataframe for transcript

    transcript = gff[gff[2].isin(['transcript'])]

    only_protein_coding = transcript[8].apply(lambda x : infogff2dict(x)['transcript_type'] == 'protein_coding')
    transcript = transcript[only_protein_coding]
    transcript['ID'] = transcript[8].apply(lambda x : x.split(';')[0].lstrip('ID='))
    selected_transcript : set = set(transcript['ID'].tolist())
    
    print(f"Number of coding transcript : {len(selected_transcript)}")

    transcript_from_selected_gene = transcript[8].apply(lambda x : infogff2dict(x)['gene_id'] in selected_gene)
    transcript = transcript[transcript_from_selected_gene]
    selected_transcript : set = set(transcript['ID'].tolist()) 

    print(f"Number of coding transcript IN a non-readthrough coding gene : {len(selected_transcript)}")

    # Working with subdataframe for CDS (to remove redundant transcript -> transcript with exactly the same CDS coordinates)

    cds : pd.DataFrame = gff[gff[2].isin(['CDS'])]

    cds_from_selected_transcript = cds[8].apply(lambda x : infogff2dict(x)['transcript_id'] in selected_transcript)
    cds = cds[cds_from_selected_transcript]

    # Making a dictionnary with gene -> transcript -> cds structure
    repartition : dict[str, dict[str, set[tuple[str, int, int]]]] = dict() 
    for index, row in cds.iterrows():
        row : list = list(row)
        chr : str = row[0]
        start : int = int(row[3])
        end : int = int(row[4])
        info = infogff2dict(row[8])
        gene_id : str = info['gene_id']
        transcript_id : str = info['transcript_id']
        if gene_id in repartition:
            if transcript_id in repartition[gene_id]:
                repartition[gene_id][transcript_id].add((chr, start, end))
            else:
                repartition[gene_id][transcript_id] = {(chr, start, end)}
        else:
            repartition[gene_id] = {transcript_id : {(chr, start, end)}}

    # Keeping only one transcript if the cds set is the same between two transcript
    repartition_no_dup : dict[str, dict[str, set[tuple[str, int, int]]]] = dict() 
    for gene in repartition:
        repartition_no_dup[gene] = dict()
        for transcript in repartition[gene]:
            if repartition[gene][transcript] not in repartition_no_dup[gene].values():
                repartition_no_dup[gene][transcript] = repartition[gene][transcript]

    selected_transcript : set = set()
    for gene_str in repartition_no_dup.values():
        selected_transcript.update(gene_str.keys())
    
    print(f"Number of non-redundant protein_coding transcript IN a non-readthrough coding gene : {len(selected_transcript)}")

    row_kept = gff.apply(keeprow, axis=1, selected_gene=selected_gene, selected_transcript=selected_transcript)
    gff = gff[row_kept]
    gff.to_csv(args.output, header=None, sep='\t', index=False)


if __name__ == '__main__':
    main()