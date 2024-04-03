"""API to process Panzzer2 output"""

import copy
from obonet import read_obo
import os
import networkx
from copy import deepcopy
from random import shuffle
from itertools import combinations
from statistics import pstdev


class Annotation:
    """Represent a set of genes."""

    def __init__(self, name) -> None:
        self.name : str = name
        self.genes : dict[str, Gene] = dict()

    def __getitem__(self, gene : str):
        return self.genes[gene]
    
    def __str__(self) -> str:
        return f"Annotation '{self.name}' with {len(self.genes)} genes"
    
    def __repr__(self) -> str:
        return f"<Annotation '{self.name}'>"
    
    def __iter__(self):  # delegate iteration to the dict
        return iter(self.genes)
    
    def __len__(self) -> int:
        return len(self.genes)


    # Setters
        
    def add_gene(self, gene):
        """gene is type Gene"""
        if gene.id in self.genes:
            raise ValueError("Trying to add a duplicate gene")
        else:
            self.genes[gene.id] = gene

    def add_transcript(self, gene, transcript):
        """gene is type string and transcript is type transcript"""
        if gene not in self.genes:
            raise ValueError("Gene not present in the annotation")
        else:
            self[gene].add_transcript(transcript)     

    # Getters

    def get_gene(self, gene : str) :
        return self[gene]
    

    # Methods

    def number_of_genes(self) -> int:
        return len(self.genes)
    
    def number_of_isoforms(self) -> int:
        return sum(x.number_of_isoforms() for x in self.genes.values())



class Gene:
    """Represent a gene in an annotation."""

    def __init__(self, id, chromosome) -> None:
        self.id : str = id
        self.transcripts : dict[str, Transcript] = dict()
        self.chromosome : str = chromosome


    def __getitem__(self, transcript : str):
        return self.transcripts[transcript]

    def __str__(self) -> str:
        return f"Gene '{self.id}' with {len(self.transcripts)} transcripts"
    def __repr__(self) -> str:
        return f"<Gene '{self.id}'>"

    # Setters
        
    def add_transcript(self, transcript):
        if transcript.id in self.transcripts:
            raise ValueError(f"Transcript {transcript.__repr__()} already in the gene")
        else:
            self.transcripts[transcript.id] = transcript

    # Getters
            
    def get_go_term_id(self, threesold=0) -> list[str]:
        ids = set()
        for transcripts in self.transcripts.values():
            ids.update(transcripts.get_go_term_id(threesold=threesold))
        return list(ids)
    
    def get_transcript(self, transcript_id : str):
        return self.transcripts[transcript_id]
    
    def get_transcripts_id(self) -> list[str]:
        return list(self.transcripts.keys())

    # Methods

    def number_of_isoforms(self) -> int:
        """Return the number of isoforms of the gene."""
        return len(self.transcripts)
    
    def stdev_number_of_go_by_isoform(self) -> float:
        """Return the POPULATION sd of the number of GO term
        by isoform."""
        data : list[int] = [iso.number_of_go_terms() for iso in self.transcripts.values()]
        return pstdev(data)
    
    def go_redundance_metric(self) -> float:
        """Custom metrics between 0 and 1 that look at GO redundance.
        If RM = 0, each GO term is present in only 1 isoform, if RM = 1, 
        each GO term is present in all isoform == all isoform are the same."""
        if self.number_of_isoforms() < 2:
            return 1.0
        if len(self.get_go_term_id()) == 0:
            return 1.0
        count_dict = dict()
        for go_term in self.get_go_term_id():
            count = -1  # init count at -1 to not count first occurence
            for isoform in self.transcripts.values():
                if go_term in isoform.get_go_term_id():
                    count += 1
            count_dict[go_term] = count
        sum_count = sum(count_dict.values())
        mean_count = sum_count / len(count_dict)
        rm = mean_count / (self.number_of_isoforms() - 1)
        return rm

    def diversity_by_pair(self, similarity_function, default_value : float = 1.0) -> float:
        """For each pair of isoform compute distance.
        Then compute the mean of all of them.
        Need a distance function
        Default value should be 1.0 if we compute a similarity and 0.0 if we compute a distance."""
        if self.number_of_isoforms() < 2:
            return default_value
        else:
            nb_combinaison = 0
            sum_distance = 0
            for isoform1, isoform2 in combinations(self.transcripts, 2):
                isoform1 = self.transcripts[isoform1]
                isoform2 = self.transcripts[isoform2]

                sum_distance += similarity_function(isoform1, isoform2)
                nb_combinaison += 1
            return sum_distance / nb_combinaison       

    def jaccard_diversity(self) -> float:
        """Mean of Jaccard index between all isoforms.
        Jaccard distance is 1 - intersection/union (using GO term)."""
        return self.diversity_by_pair(jaccard_index)
    
    def intra_gene_gogo_similarity(self, gogo_dir : str, gogo_file : str = "gogo") -> tuple[float, float, float]:
        """
        Return a dict of key = transcript_id_1_X_transcript_id_2 and value = dict(ontology, similary)
        For example, sim['ENST0001_X_ENST0002']['BP'] is the similarity of the BP GO term between both ENST0001 and ENST0002 transcript(for the two input annotation)
        """
        # check number of isoform

        if len(self.transcripts) < 2:
            return (1, 1, 1)
        
        bp_sum : float = 0
        cc_sum : float = 0
        mf_sum : float = 0
        nb_pair : int = 0
        with open(f"{gogo_file}.gogo_input.txt", 'w') as input_file:
            for transcript1, transcript2 in combinations(self.transcripts, 2):
                nb_pair += 1
                line_p1 = transcript1 + f"{nb_pair} "
                line_p2 = "; " + transcript2 + f"{nb_pair} "
                bp_t1 = self[transcript1].get_go_term_id('BP')
                bp_t2 = self[transcript2].get_go_term_id('BP')
                cc_t1 = self[transcript1].get_go_term_id('CC')
                cc_t2 = self[transcript2].get_go_term_id('CC')
                mf_t1 = self[transcript1].get_go_term_id('MF')
                mf_t2 = self[transcript2].get_go_term_id('MF')

                if set(bp_t1) == set(bp_t2):  # include both length == 0
                    bp_sum += 1
                elif len(bp_t1) == 0 or len(bp_t2) == 0:
                    bp_sum += 0
                else:
                    line_p1 += ' '.join(bp_t1) + ' '
                    line_p2 += ' '.join(bp_t2) + ' '
                
                if set(cc_t1) == set(cc_t2):  # include both length == 0
                    cc_sum += 1
                elif len(cc_t1) == 0 or len(cc_t2) == 0:
                    cc_sum += 0
                else:
                    line_p1 += ' '.join(cc_t1) + ' '
                    line_p2 += ' '.join(cc_t2) + ' '

                if set(mf_t1) == set(mf_t2):  # include both length == 0
                    mf_sum += 1
                elif len(mf_t1) == 0 or len(mf_t2) == 0:
                    mf_sum += 0
                else:
                    line_p1 += ' '.join(mf_t1) + ' '
                    line_p2 += ' '.join(mf_t2) + ' '
                
                input_file.write(line_p1 + line_p2 +';\n')

        run_gogo_script(gogo_dir, gogo_file)

        with open(f"{gogo_file}.gogo_output.txt") as output_file:
            for line in output_file:
                line = line.strip().split(' ')
                bp_sum += 0 if line[-5] == 'NA' else float(line[-5])
                cc_sum += 0 if line[-3] == 'NA' else float(line[-3])
                mf_sum += 0 if line[-1] == 'NA' else float(line[-1])
        return bp_sum/nb_pair, cc_sum/nb_pair, mf_sum/nb_pair

class Transcript:
    """Represent a transcripts."""

    def __init__(self, id) -> None:
        self.id : str = id
        self.gos : list[GO] = list()
        self.keggs : list[KEGG] = list()
        self.ecs : list[EC] = list()
        self.seq : str = None
        self.desc : str = None  # description
        self.ppv = None

    def __str__(self) -> str:
        return f"Transcript '{self.id}' with {len(self.gos)} GO terms"
    
    def __repr__(self) -> str:
        return f"<Transcript '{self.id}'>"
    
    # Getters

    def get_go_term_id(self, namespace :str = 'all', threesold : float = 0) -> list[str]:
        """Return a list of ids of the GO term assigned to this transcript."""
        if namespace == 'all':
            return [go.id for go in self.gos if go.ppv > threesold]
        else:
            return [go.id for go in self.gos if (go.ppv > threesold and go.namespace == namespace)]

    # Setters
        
    def add_go(self, go):
        """go is type GO"""
        self.gos.append(go)

    def add_ec(self, ec):
        """ec is type EC"""
        self.ecs.append(ec)

    def add_kegg(self, kegg):
        """kegg is type KEGG"""
        self.keggs.append(kegg)

    def set_seq(self, seq : str):
        self.seq = seq

    def set_desc(self, desc : str, ppv):
        self.desc = desc
        self.ppv = float(ppv)

    # Methods
        
    def number_of_go_terms(self) -> int:
        return len(self.gos)


class Feature:
    """Represent a feature."""

    def __init__(self, id, desc, ppv) -> None:
        self.id = id
        self.desc = desc
        self.ppv = float(ppv)

    def __str__(self):
        return self.id

    def __repr__(self):
        return self.id   

class GO(Feature):

    def __init__(self, id, desc, ppv, namespace) -> None:
        """namespace is BP, MF, CC"""
        self.namespace = namespace
        super().__init__(id, desc, ppv)

class KEGG(Feature):
    pass

class EC(Feature):
    pass
    # def relation(self, ec):
    #     """ec is another EC
    #     s
    #     return 'same' is its the same EC
    #     return 'self' if the self EC is more precise
    #     return 'arg' if the other EC is more precise
    #     return 'diff' if the two are different and don't include each other"""
    #     if self.id == ec.id:
    #         return 'same'
    #     else:
    #         self_id = self.id.split('.')
    #         alter_id = ec.id.split('.')
    #         for s,a in zip(self_id, alter_id):
    #             if s == a:
    #                 pass
    #             if s == '-':
    #                 return 'arg'
    #             if a == '-':
    #                 return 'self'
    #             if a != s:
    #                 return 'diff'



def parse_pannzer_annotation(path, name = 'undefined_name'):
    """
    Return an Annotation object from a PANZZER output
    """
    annotation = Annotation(name)
    with open(path) as file:
        for line in file:
            line = line.strip().split('\t')
            match line[1]:  # line[1] is type of line
                case "original_DE":
                    # get gene info from sequence description
                    sequence_description = line[5].split(' ')
                    gene_id = sequence_description[0].split('=')[1].split('.')[0]
                    chr = sequence_description[1].split('=')[1]

                    # creating gene and adding him to the annotation if not present
                    if gene_id not in annotation.genes:
                        gene = Gene(gene_id, chr)
                        annotation.add_gene(gene)

                    # creating current transcript
                    transcript_id = line[0]
                    transcript = Transcript(transcript_id)
                    annotation.add_transcript(gene_id, transcript)
                case "qseq":
                    # sequence available
                    transcript.set_seq(line[-1])
                case "DE":
                    # prediction of the description of the transcript
                    transcript.set_desc(line[-1], line[3])
                case "BP_ARGOT":
                    id = 'GO:' + line[-2]
                    desc = line[-1]
                    ppv = float(line[-3])
                    feature = GO(id, desc, ppv, 'BP')
                    transcript.add_go(feature)
                case "CC_ARGOT":
                    id = 'GO:' + line[-2]
                    desc = line[-1]
                    ppv = float(line[-3])
                    feature = GO(id, desc, ppv, 'CC')
                    transcript.add_go(feature)
                case "MF_ARGOT":
                    id = 'GO:' + line[-2]
                    desc = line[-1]
                    ppv = float(line[-3])
                    feature = GO(id, desc, ppv, 'MF')
                    transcript.add_go(feature)
                case "EC_ARGOT":
                    id = line[-2]
                    desc = line[-1]
                    ppv = float(line[-3])
                    feature = EC(id, desc, ppv)
                    transcript.add_ec(feature) 
                case "KEGG_ARGOT":  
                    id = line[-2]
                    desc = line[-1]
                    ppv = float(line[-3])
                    feature = KEGG(id, desc, ppv)
                    transcript.add_kegg(feature)

    return annotation


def make_best_single_isoform_annotation(annotation : Annotation):
    """Return a single-isoform annotation from a multiple one by taking, for each gene,
    the isoform with the greater number of GO term"""
    best_annotation = Annotation(annotation.name + "_best")
    for gene in annotation.genes:
        id = annotation[gene].id
        chr = annotation[gene].chromosome
        new_gene = Gene(id, chr)
        best_isoform = max(annotation[gene].transcripts.values(), key= lambda x : len(x.gos))  # select isoform with highest number of GO terms
        new_gene.add_transcript(copy.deepcopy(best_isoform))
        best_annotation.add_gene(new_gene)
    return best_annotation

def jaccard_index(isoform1 : Transcript, isoform2 : Transcript) -> float:
    """ Return Jaccard distance between 2 isoform.
    Jaccard distance is intersection/union (using GO term)"""
    go_term1 = set(isoform1.get_go_term_id())
    go_term2 = set(isoform2.get_go_term_id())
    intersection = set.intersection(go_term1, go_term2)
    union = set.union(go_term1, go_term2)
    if len(union) == 0:
        return 1.0
    else:
        return len(intersection) / len(union)
    
def dice_coefficient(isoform1 : Transcript, isoform2 : Transcript) -> float:
    """Return Sorensen-Dice coefficient"""
    go_term1 = set(isoform1.get_go_term_id())
    go_term2 = set(isoform2.get_go_term_id())
    intersection = set.intersection(go_term1, go_term2)
    if len(go_term1) + len(go_term2) == 0:
        return 1.0
    else:
        return (2 * len(intersection)) / (len(go_term1) + len(go_term2))

def overlap_coefficient(isoform1 : Transcript, isoform2 : Transcript) -> float:
    """Return Sorensen-Dice coefficient"""
    go_term1 = set(isoform1.get_go_term_id())
    go_term2 = set(isoform2.get_go_term_id())
    intersection = set.intersection(go_term1, go_term2)
    if min(len(go_term1), len(go_term2)) == 0:
        return 1.0
    else:
        return len(intersection) / min(len(go_term1), len(go_term2))

 
def read_go_graph():
    url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    go_graph = read_obo(url, ignore_obsolete=False)
    go_graph.remove_edges_from([edge for edge in go_graph.edges if (edge[2] != 'is_a' and edge[2] != 'part_of')])
    go_graph.remove_nodes_from([node for node in go_graph.nodes if node in ['GO:0005575', 'GO:0008150', 'GO:0003674']])
    return go_graph


def go_list_to_go_ancestry(go_list, go_graph):
    """
        go_list is a list(str) with each string formatted GO:XXXXXX
        go_graph is a networkx graph from the .obo gene ontology graph
        return a list with all go term and the parent recursively
    """
    go_set = set(go_list)
    for id in go_list:
        go_set.update(networkx.descendants(go_graph,id))
    return list(go_set)

def genes_with_diff_go_terms(annotation_1, annotation_2, gene_set = None):
    """
        Output the list of gene with a different number of GO term between annotation
    """
    if gene_set is None:  # if no geneset, create a geneset of the common gene between annotation
        gene_set = [gene for gene in annotation_1.genes if gene in annotation_2.genes]
    gene_list = []
    for gene in gene_set:
        goterm_1 = set(annotation_1[gene].get_go_term_id())
        goterm_2 = set(annotation_2[gene].get_go_term_id())
        if goterm_1 !=  goterm_2:
            gene_list.append(gene)
    return gene_list

def genes_with_diff_go_terms_with_hierarchy(annotation_1, annotation_2, go_graph, which = (True, True), gene_set = None):
    """
        Output the list of gene with a different number of GO term between annotation, considering hierarchy
        go_graph is a networkx graph of go term
        which is a tuple of boolean specifing if the hierarchy need to be compute for annotation1 and annotation2 respectively
    """
    if gene_set is None:  # if no geneset, create a geneset of the common gene between annotation
        gene_set = [gene for gene in annotation_1.genes if gene in annotation_2.genes]
    gene_list = []
    for gene in gene_set:
        goterm_l1 = annotation_1[gene].get_go_term_id()
        if which[0]:
            goterm_l1 = go_list_to_go_ancestry(goterm_l1, go_graph)
        goterm_l2 = annotation_2[gene].get_go_term_id()
        if which[1]:
            goterm_l2 = go_list_to_go_ancestry(goterm_l2, go_graph)
        if set(goterm_l1) != set(goterm_l2):
            gene_list.append(gene)
    return gene_list

def run_gogo_script(gogo_dir : str, gogo_file : str) -> None:
    os.system(f"""
        CURDIR=$(pwd)
        cd {gogo_dir}
        perl gene_pair_comb.pl $CURDIR/{gogo_file}.gogo_input.txt $CURDIR/{gogo_file}.gogo_output.txt 
        cd $CURDIR
        """)   

def gogo_similarity_between_transcript(transcript_1 : Transcript, transcript_2 : Transcript, gogo_dir : str, gogo_file : str = "gogo") -> tuple[float, float, float]:
    bp_sim = gogo_similarity_between_set_from_the_same_namespace(set(transcript_1.get_go_term_id(namespace='BP')), set(transcript_2.get_go_term_id(namespace='BP')), gogo_dir, gogo_file)
    cc_sim = gogo_similarity_between_set_from_the_same_namespace(set(transcript_1.get_go_term_id(namespace='CC')), set(transcript_2.get_go_term_id(namespace='CC')), gogo_dir, gogo_file)
    mf_sim = gogo_similarity_between_set_from_the_same_namespace(set(transcript_1.get_go_term_id(namespace='MF')), set(transcript_2.get_go_term_id(namespace='MF')), gogo_dir, gogo_file)
    return bp_sim, cc_sim, mf_sim

def gogo_similarity_between_set_from_the_same_namespace(set_1 : set[str], set_2 : set[str], gogo_dir : str, gogo_file : str = "gogo") -> float:
    """
    Given 2 set of GO term, compute the similarity with GOGO
    BOTH SET NEED TO BE FROM THE SAME NAMESPACE (BP, CC or MF)
    """
    if set_1 == set_2:  # include len(set_1) == 0 and len(set_2) == 0
        return 1
    if len(set_1) == 0 or len(set_2) == 0:
        return 0
    
    # write input file for the bash script
    with open(f"{gogo_file}.gogo_input.txt", 'w') as input_file:
        line = 'set_1 '
        line += ' '.join(list(set_1)) + "; "
        line += 'set2 '
        line += ' '.join(list(set_2)) + ";"
        input_file.write(line) 

    # running bash script
    run_gogo_script(gogo_dir, gogo_file)
    
    # parsing output script
    with open(f"{gogo_file}.gogo_output.txt") as output_file:
        line = output_file.readline()
        similarity = line.split(';')[-1].replace('BPO', '').replace('CCO', '').replace('MFO', '').replace('NA', '').strip()

    return float(similarity)


def gogo_similarity_between_annotation(annotatation_1 : Annotation, annotation_2 : Annotation, gogo_dir:str, gene_set = None, gogo_file : str = "gogo") -> dict[str, dict[str, str]]:
    """
        Return a dict of key = gene_id and value = dict(ontology, similary)
        For example, sim['ENSG0001']['BP'] is the similarity of the BP GO term for the gene ENSG0001 (for the two input annotation)
    """
    if gene_set is None: # if no geneset, create a geneset of the common gene between annotation
        gene_set = [gene for gene in annotatation_1.genes if gene in annotation_2.genes]

    # write input file for the bash script
    with open(f"{gogo_file}.gogo_input.txt", 'w') as input_file:
        for gene in gene_set:
            id = annotatation_1[gene].id
            line = id + '-' + annotatation_1.name + ' '
            line = line + ' '.join(annotatation_1[gene].get_go_term_id()) + "; "
            line = line + id + '-' + annotation_2.name + ' '
            line = line + ' '.join(annotation_2[gene].get_go_term_id()) + "; \n"
            input_file.write(line)  # gene1 GO:0001 GO:0002; gene2 GO:0001 GO:0002

    # running bash script
    run_gogo_script(gogo_dir, gogo_file)
    
    # parsing output script
    similarity = dict()
    with open(f"{gogo_file}.gogo_output.txt") as output_file:
        for line, gene in zip(output_file, gene_set):
            # line = line.replace('NA','1.000') # do not do that, if geneA has no go term but geneB have, NA will be displayed
            line = line.strip().split(' ')
            similarity[gene] = dict()
            similarity[gene]['BP'] = line[-5]
            similarity[gene]['CC'] = line[-3]
            similarity[gene]['MF'] = line[-1]
    return similarity

# def gogo_write_input(annotatation_1 : Annotation, annotation_2 : Annotation, path: str, gene_set = None):
#     if gene_set is None: # if no geneset, create a geneset of the common gene between annotation
#         gene_set = [gene for gene in annotatation_1.genes if gene in annotation_2.genes]

#     # write input file for the bash script
#     with open(path, 'w') as input_file:
#         for gene in gene_set:
#             id = annotatation_1[gene].id
#             line = id + '-' + annotatation_1.name + ' GO:'
#             line = line + ' GO:'.join(annotatation_1[gene].get_go_term_id()) + "; "
#             line = line + id + '-' + annotation_2.name + ' GO:'
#             line = line + ' GO:'.join(annotation_2[gene].get_go_term_id()) + "; \n"
#             input_file.write(line)  # gene1 GO:0001 GO:0002; gene2 GO:0001 GO:0002


# def gogo_read_output(output_path : str, gene_set):
#     similarity = dict()
#     with open(output_path) as output_file:
#         for line, gene in zip(output_file, gene_set):
#             # line = line.replace('NA','1.000') # do not do that, if geneA has no go term but geneB have, NA will be displayed
#             line = line.strip().split(' ')
#             similarity[gene] = dict()
#             similarity[gene]['BP'] = line[-5]
#             similarity[gene]['CC'] = line[-3]
#             similarity[gene]['MF'] = line[-1]
#     return similarity

def gogo_similarity_between_annotation_with_hierarchy(annotatation_1 : Annotation, annotation_2 : Annotation, gogo_dir:str, go_graph, gogo_file : str = "gogo", which = (True, True), gene_set = None)  -> dict[str, dict[str, str]]:
    """
        Return a dict of key = gene_id and value = dict(ontology, similary)
        For example, sim['ENSG0001']['BP'] is the similarity of the BP GO term for the gene ENSG0001 (for the two input annotation)
    """
    if gene_set is None: # if no geneset, create a geneset of the common gene between annotation
        gene_set = [gene for gene in annotatation_1.genes if gene in annotation_2.genes]

    # write input file for the bash script
    with open(f"{gogo_file}.gogo_input.txt", 'w') as input_file:
        for gene in gene_set:
            id = annotatation_1[gene].id
            line = id + '-' + annotatation_1.name + '_A1 '
            go_term_annotation_1 = annotatation_1[gene].get_go_term_id()
            if which[0]:
                go_term_annotation_1 = go_list_to_go_ancestry(go_term_annotation_1, go_graph)
            line = line + ' '.join(go_term_annotation_1) + "; "
            line = line + id + '-' + annotation_2.name + '_A2 '
            go_term_annotation_2 = annotation_2[gene].get_go_term_id()
            if which[1]:
                go_term_annotation_2 = go_list_to_go_ancestry(go_term_annotation_2, go_graph)
            line = line + ' '.join(go_term_annotation_2) + "; \n"
            input_file.write(line)  # gene1 GO:0001 GO:0002; gene2 GO:0001 GO:0002

    # running bash script
    run_gogo_script(gogo_dir, gogo_file)
    
    # parsing output script
    similarity = dict()
    with open(f"{gogo_file}.gogo_output.txt") as output_file:
        for line, gene in zip(output_file, gene_set):
            # line = line.replace('NA','1.000') # do not do that, if geneA has no go term but geneB have, NA will be displayed
            line = line.strip().split(' ')
            similarity[gene] = dict()
            similarity[gene]['BP'] = line[-5]
            similarity[gene]['CC'] = line[-3]
            similarity[gene]['MF'] = line[-1]
    return similarity

def mean_similarity(similarity):
    result = dict()
    sBP = nBP = 0
    sCC = nCC = 0
    sMF = nMF = 0
    for gene in similarity:
        if similarity[gene]["BP"] != 'NA':
            sBP += float(similarity[gene]["BP"])
            nBP += 1
        if similarity[gene]["CC"] != 'NA':
            sCC += float(similarity[gene]["CC"])
            nCC += 1
        if similarity[gene]["MF"] != 'NA':
            sMF += float(similarity[gene]["MF"])
            nMF += 1
    sBP /= nBP
    result['BP'] = sBP
    sCC /= nCC
    result['CC'] = sCC
    sMF /= nMF
    result['MF'] = sMF
    return result


def make_shuffle_copy(annotation : Annotation):
    copied = deepcopy(annotation)
    transcript_list = []
    for gene in copied.genes.values():
    #  creating a list of transcripts
        transcript_list.extend(gene.transcripts.values())

    shuffle(transcript_list)  # shuffle the list

    index = 0
    for gene in copied.genes.values():
        n_transcript = len(gene.transcripts)
        gene.transcripts.clear()  # removing current transcripts
        for transcript in transcript_list[index:index+n_transcript]:
            gene.add_transcript(transcript)  # adding random transcripts
        index += n_transcript
    copied.name = annotation.name + "_shuffle"  # important when computing similarity
    return copied


def make_longest_single_isoform_annotation(annotation : Annotation):
    """Return a single-isoform annotation from a multiple one by taking, for each gene,
    the isoform with the longest sequence"""
    long_annotation = Annotation(annotation.name + "_long")
    for gene in annotation.genes:
        id = annotation[gene].id
        chr = annotation[gene].chromosome
        new_gene = Gene(id, chr)
        long_isoform = max(annotation[gene].transcripts.values(), key= lambda x : len(x.seq))  # select isoform with highest number of GO terms
        new_gene.add_transcript(deepcopy(long_isoform))
        long_annotation.add_gene(new_gene)
    return long_annotation

