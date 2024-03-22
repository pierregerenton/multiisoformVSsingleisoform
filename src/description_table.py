import argparse
import pannzer_out_api as poa


parser = argparse.ArgumentParser(
                    prog='description_table',
                    description='Write a table with isoform information (reformating of Pannzer Output)\
                        and a table of similarity with metadata for a multi-isoform annotation vs singles',
                    epilog='For more information, contact pierre.gerenton@crg.eu')


parser.add_argument(
    '-m', '--multi-isoform',
    type=str, required=True,
    help='Path of a panzzer output (from a multiple isoform annotation)'
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

parser.add_argument(
    '-c', '--custom-single-isoform',
    type=str, required=False, default=None,
    help='Path of a panzzer output (from a single isoform annotation)'
)

parser.add_argument(
    '-l', '--longest-single-isoform',
    required=False, default=False, action="store_true",
    help='Computing similarity against the longest isoform'
)

parser.add_argument(
    '-b', '--best-single-isoform',
    required=False, default=False, action="store_true",
    help='Computing similarity against the best isoform'
)

args = parser.parse_args()



def is_similarity_table_needed(args) -> bool:
    c = 0
    if args.custom_single_isoform is not None:
        c += 1
    if args.longest_single_isoform:
        c += 1
    if args.best_single_isoform:
        c += 1
    if c == 0:
        return False
    return True

def main():

    similarity_table = is_similarity_table_needed(args)  # checking if an argument is available for gogo similarity

    multi_isoform = poa.parse_pannzer_annotation(args.multi_isoform, 'multi')

    print('Reading alternative single isoform annotation...')
    if args.custom_single_isoform is not None:
        custom_isoform = poa.parse_pannzer_annotation(args.custom_single_isoform, 'custom')
    if args.longest_single_isoform:
        longest_isoform = poa.make_longest_single_isoform_annotation(multi_isoform)
    if args.best_single_isoform:
        best_isoform = poa.make_best_single_isoform_annotation(multi_isoform)
    print('Done')
    print()

    if similarity_table:
        print('Preparing gene set ...')
        gene_set= list(multi_isoform.genes)  # getting genes from multi isoform

        if args.custom_single_isoform is not None: # keeping genes if present in our custom annotation
            gene_set = [gene for gene in custom_isoform.genes if gene in gene_set]

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
    
        print('Computing similarity ...')

        gogo_name = args.multi_isoform.split('.')[0]
        if args.custom_single_isoform is not None:
            sim_custom = fun(multi_isoform, custom_isoform, args.gogo_dir, gogo_file=gogo_name + ".custom", **argument)
        if args.longest_single_isoform:
            sim_longest = fun(multi_isoform, longest_isoform, args.gogo_dir, gogo_file=gogo_name + ".longest", **argument)
        if args.best_single_isoform:
            sim_best = fun(multi_isoform, best_isoform, args.gogo_dir, gogo_file=gogo_name + ".best", **argument)

        print('Done')
        print()

        print('Writing similarity table ...')
        with open(args.output_name + '.similarity.txt','w') as sim_output:
          
            header = "chromosome\tgene_id\tnb_isoform\tlist_of_go_term_assigned_to_the_gene\tnb_go_term_assigned_to_the_gene"
            if args.custom_single_isoform is not None:
                header += "\ttranscript_id_custom\tdescription_custom\tscore_custom\tlist_of_go_term_assigned_to_the_transcript_custom"\
                          "\tnb_go_term_assigned_to_the_transcript_custom\tBP_similarity_custom\tCC_similarity_custom\tMF_similarity_custom"
            if args.longest_single_isoform:
                header += "\ttranscript_id_longest\tdescription_longest\tscore_longest\tlist_of_go_term_assigned_to_the_transcript_longest"\
                          "\tnb_go_term_assigned_to_the_transcript_longest\tBP_similarity_longest\tCC_similarity_longest\tMF_similarity_longest"                
            if args.best_single_isoform:
                header += "\ttranscript_id_best\tdescription_best\tscore_best\tlist_of_go_term_assigned_to_the_transcript_best"\
                          "\tnb_go_term_assigned_to_the_transcript_best\tBP_similarity_best\tCC_similarity_best\tMF_similarity_best"                
            header += '\n'
            sim_output.write(header)

            for gene in gene_set:
                gene = multi_isoform.get_gene(gene)
                chromosome = gene.chromosome
                gene_id = gene.id
                nb_isoform = len(gene.transcripts)
                list_of_go_term_assigned_to_the_gene = gene.get_go_term_id()
                nb_go_term_assigned_to_the_gene = len(list_of_go_term_assigned_to_the_gene)
                list_of_go_term_assigned_to_the_gene = ';'.join(list_of_go_term_assigned_to_the_gene)
                line = f"{chromosome}\t{gene_id}\t{nb_isoform}\t{list_of_go_term_assigned_to_the_gene}\t{nb_go_term_assigned_to_the_gene}"
                if args.custom_single_isoform is not None:
                    transcript_id_custom = custom_isoform.genes[gene_id].get_transcripts_id()[0]
                    transcript = custom_isoform.genes[gene_id].get_transcript(transcript_id_custom)
                    description_custom = transcript.desc
                    score_custom = transcript.ppv
                    list_of_go_term_assigned_to_the_transcript_custom = transcript.get_go_term_id()
                    nb_go_term_assigned_to_the_transcript_custom = len(list_of_go_term_assigned_to_the_transcript_custom)
                    list_of_go_term_assigned_to_the_transcript_custom = ';'.join(list_of_go_term_assigned_to_the_transcript_custom)
                    BP_similarity_custom = sim_custom[gene_id]['BP']
                    CC_similarity_custom = sim_custom[gene_id]['CC']
                    MF_similarity_custom = sim_custom[gene_id]['MF']
                    line += f"\t{transcript_id_custom}\t{description_custom}\t{score_custom}\t{list_of_go_term_assigned_to_the_transcript_custom}"\
                          f"\t{nb_go_term_assigned_to_the_transcript_custom}\t{BP_similarity_custom}\t{CC_similarity_custom}\t{MF_similarity_custom}"
                if args.longest_single_isoform:
                    transcript_id_longest = longest_isoform.genes[gene_id].get_transcripts_id()[0]
                    transcript = longest_isoform.genes[gene_id].get_transcript(transcript_id_longest)
                    description_longest = transcript.desc
                    score_longest = transcript.ppv
                    list_of_go_term_assigned_to_the_transcript_longest = transcript.get_go_term_id()
                    nb_go_term_assigned_to_the_transcript_longest = len(list_of_go_term_assigned_to_the_transcript_longest)
                    list_of_go_term_assigned_to_the_transcript_longest = ';'.join(list_of_go_term_assigned_to_the_transcript_longest)
                    BP_similarity_longest = sim_longest[gene_id]['BP']
                    CC_similarity_longest = sim_longest[gene_id]['CC']
                    MF_similarity_longest = sim_longest[gene_id]['MF']
                    line += f"\t{transcript_id_longest}\t{description_longest}\t{score_longest}\t{list_of_go_term_assigned_to_the_transcript_longest}"\
                          f"\t{nb_go_term_assigned_to_the_transcript_longest}\t{BP_similarity_longest}\t{CC_similarity_longest}\t{MF_similarity_longest}"
                if args.best_single_isoform:
                    transcript_id_best = best_isoform.genes[gene_id].get_transcripts_id()[0]
                    transcript = best_isoform.genes[gene_id].get_transcript(transcript_id_best)
                    description_best = transcript.desc
                    score_best = transcript.ppv
                    list_of_go_term_assigned_to_the_transcript_best = transcript.get_go_term_id()
                    nb_go_term_assigned_to_the_transcript_best = len(list_of_go_term_assigned_to_the_transcript_best)
                    list_of_go_term_assigned_to_the_transcript_best = ';'.join(list_of_go_term_assigned_to_the_transcript_best)
                    BP_similarity_best = sim_best[gene_id]['BP']
                    CC_similarity_best = sim_best[gene_id]['CC']
                    MF_similarity_best = sim_best[gene_id]['MF']
                    line += f"\t{transcript_id_best}\t{description_best}\t{score_best}\t{list_of_go_term_assigned_to_the_transcript_best}"\
                          f"\t{nb_go_term_assigned_to_the_transcript_best}\t{BP_similarity_best}\t{CC_similarity_best}\t{MF_similarity_best}"
                line += '\n'
                sim_output.write(line)

    
 
        print('Done !')
        print()

    print('Writing metadata table ...')
      
    with open(args.output_name + '.metadata.txt','w') as txt_output:
        header = "chromosome\tgene_id\tnb_isoform\tlist_of_go_term_assigned_to_the_gene\tnb_go_term_assigned_to_the_gene\ttranscript_id\t"\
                 "description\tscore\tlist_of_go_term_assigned_to_the_transcript\tnb_go_term_assigned_to_the_transcript\n"
        txt_output.write(header)
        for gene in multi_isoform.genes:
            gene = multi_isoform.get_gene(gene)
            chromosome = gene.chromosome
            gene_id = gene.id
            nb_isoform = len(gene.transcripts)
            list_of_go_term_assigned_to_the_gene = gene.get_go_term_id()
            nb_go_term_assigned_to_the_gene = len(list_of_go_term_assigned_to_the_gene)
            list_of_go_term_assigned_to_the_gene = ';'.join(list_of_go_term_assigned_to_the_gene)
            for transcript in gene.transcripts.values():
                transcript_id = transcript.id
                description = transcript.desc
                score = transcript.ppv
                list_of_go_term_assigned_to_the_transcript = transcript.get_go_term_id()
                nb_go_term_assigned_to_the_transcript = len(list_of_go_term_assigned_to_the_transcript)
                list_of_go_term_assigned_to_the_transcript = ';'.join(list_of_go_term_assigned_to_the_transcript)
                txt_output.write(
                    f"{chromosome}\t{gene_id}\t{nb_isoform}\t{list_of_go_term_assigned_to_the_gene}\t{nb_go_term_assigned_to_the_gene}\t{transcript_id}\t"\
                    f"{description}\t{score}\t{list_of_go_term_assigned_to_the_transcript}\t{nb_go_term_assigned_to_the_transcript}\n"
                )

    print('Done !')


if __name__=="__main__":
    main()


