import os
import json
from tqdm import tqdm
import pandas as pd
from operator import itemgetter
import obonet
import GOA
import networkx as nx


def file_length(filename):
    return sum(1 for line in open(filename, 'rb'))


def filter_evidence(annot_file, save_location, highthr=True):

    exp_evidence = ['EXP', 'IPI', 'IDA', 'IMP', 'IGI', 'IEP']
    inferred_evidence = ['TAS', 'IC']
    h_evidence = ['HTP', 'HDA', 'HMP', 'HGI', 'HEP']
    ok_evidence = exp_evidence + inferred_evidence
    if highthr:
        ok_evidence += h_evidence
    
    gaf_filelen = file_length(annot_file)

    with open(annot_file) as handle, open(save_location, 'a', newline='') as f:
        pbar = tqdm(GOA.gafiterator(handle), total=gaf_filelen, unit_scale=True)
        for rec in pbar:
            if rec['Evidence'] in ok_evidence and rec['DB']=='UniProtKB':
                # save full records
                #with open(save_location, 'a', newline='') as f:
                json.dump(rec, f)
                f.write(os.linesep)

    return pd.read_json(save_location, lines=True)


def clean_annotations(filtered_annotations):
    """ because some entried have multiple values of the column DB:Reference,
        some terms (rows) are duplicated. We will get rid of these duplicates
        and other columns we don't need, in particular columns that have lists
        as values, which will prevent pandas from finding duplicates. 
        We will also get rid of any negative labels ('NOT' in the qualifier)"""

    annotations_df = filtered_annotations.drop(['DB:Reference', 'Synonym', 'With',
                                                'Annotation_Extension', 'Gene_Product_Form_ID'], axis=1)

    # The Qualifier column is a list that can contain two entries if the terms is NOT something
    annotations_df['not_Qualifier'] = annotations_df.Qualifier.apply(lambda x: 'NOT' in x)
    annotations_df['qualifier'] = annotations_df.Qualifier.apply(lambda x: x[-1])

    # process taxonID to get just the number, eg '[taxon:9606]' becomes '9606'
    annotations_df['species']=annotations_df.Taxon_ID.apply(lambda x: x[0].split(':')[-1])

    annotations_df = annotations_df.drop(['Qualifier', 'Taxon_ID'],axis=1)
    annotations_df = annotations_df.drop_duplicates()

    # Exclude negative lables 
    annotations_df = annotations_df[~annotations_df.not_Qualifier]

    return annotations_df 


def propagate_terms(df, obo_file):
    obsolete_replace={'GO:0006975':'GO:0042770',
                     'GO:0031617':'GO:0000776',
                     'GO:1901720':'GO:1905560',
                     'GO:0034291':'GO:0140911',
                     'GO:0034290':'GO:0140911',
                     'GO:0034292':'GO:0140911',
                     'GO:0004147':'GO:0043754',
                     'GO:0044662':'GO:0051673',
                     'GO:0044649':'GO:0051715',
                     'GO:0050828':'GO:0043129',
                     'GO:1990142':'GO:0044179',
                     'GO:0102430':'GO:0102431',
                     'GO:0006295': 'GO:0006289',
                     'GO:0006296': 'GO:0006289',
                     'GO:0006875': 'GO:0030003',
                     'GO:0008022': 'GO:0005515',
                     'GO:0008852': 'GO:0008310',
                     'GO:0008853': 'GO:0008311',
                     'GO:0030004': 'GO:0030003',
                     'GO:0030320': 'GO:0030002',
                     'GO:0031997': 'GO:0005515',
                     'GO:0032199': 'GO:0032197',
                     'GO:0033683': 'GO:0006289',
                     'GO:0046916': 'GO:0030003',
                     'GO:0047485': 'GO:0005515',
                     'GO:0072507': 'GO:0055080',
                     'GO:0097056': 'GO:0001717',
                     'GO:1904608': 'GO:1902065',
                     'GO:1990731': 'GO:0070914'}

    # load graph
    obonet_graph = obonet.read_obo(obo_file)
    # keep only "is_a" and "part_of" edges
    remove_edges = [(i, j, k) for i, j, k in obonet_graph.edges if not(k=="is_a" or k=="part_of")]
    obonet_graph.remove_edges_from(remove_edges)

    # replace obsolete
    df.GO_ID = df.GO_ID.replace(obsolete_replace)
    # remove remaining obsolete
    obsolete_remove = set(df.GO_ID.values).difference(set(obonet_graph.nodes))

    # look up ancestors once to save time
    ancestor_lookup = {t: nx.descendants(obonet_graph,t) for t in df.GO_ID.unique()
                       if t in obonet_graph}
    obsolete_lookup = {t: set() for t in obsolete_remove}
    ancestor_lookup.update(obsolete_lookup)

    propagated_terms = []
    for aspect, subontology in zip(['P','C','F'], ['BPO', 'CCO', 'MFO']):

        print(f'Processing {subontology} annotations')

        for gene, annotations in df[df.Aspect==aspect].groupby('DB_Object_ID'):
            terms = annotations.GO_ID.values
            if len(terms) > 1:
                gene_terms = set(terms).union(*itemgetter(*set(terms))(ancestor_lookup))
            else:
                gene_terms = set([terms[0]]).union(ancestor_lookup[terms[0]])

            gene_terms -= obsolete_remove
            if len(gene_terms) >= 1:  # perhaps gene only had obsolete terms
                propagated_terms.append({'EntryID': gene, 'term': list(gene_terms), 'aspect':subontology})

    return pd.DataFrame(propagated_terms)
