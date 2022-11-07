import os
import sys
import json   
import GOA
from tqdm import tqdm
import urllib.request
import gzip
import shutil
import argparse

import pandas as pd
import networkx as nx
import obonet


def file_length(filename):
    return sum(1 for line in open(filename, 'rb'))


def download_gofile(source_path, save_path):
    if not os.path.exists(save_path):
        urllib.request.urlretrieve(source_path, filename=os.path.join(save_path))
    
    out_file = save_path.split('.gz')[0]
    if not os.path.exists(out_file):
        with gzip.open(save_path, 'rb') as f_in:
            with open(out_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    return out_file


def filter_evidence(annot_file, save_location):
    
    count = 0
    exp_evidence = ['EXP', 'IPI', 'IDA', 'IMP', 'IGI', 'IEP']
    inferred_evidence = ['TAS', 'IC']
    h_evidence = ['HTP', 'HDA', 'HMP', 'HGI', 'HEP']
    ok_evicende = exp_evidence + inferred_evidence + h_evidence
    gaf_filelen = file_length(annot_file)

    with open(annot_file) as handle:
        pbar = tqdm(GOA.gafiterator(handle), total=gaf_filelen, unit_scale=True)
        for rec in pbar:
            pbar.set_description("Matches Found {}".format(count))
            if rec['Evidence'] in ok_evidence:
                count += 1
                rec_info = {'GeneID': rec['DB_Object_ID'], 'term': rec['GO_ID'], 'evidence': rec['Evidence'],
                            'aspect': rec['Aspect'], 'synonym': rec['Synonym'], 'symbol': rec['DB_Object_Symbol'],
                            'taxon': rec['Taxon_ID'], 'date': rec['Date']}
                with open(save_location, 'a', newline='') as f:
                    json.dump(rec_info, f)
                    f.write(os.linesep)     


def propogate_terms(annotation_file, save_location, db_name):

    # a couple known term replacements
    obsolete_lookup={'GO:0006975':'GO:0042770', 
                     'GO:0031617':'GO:0000776',
                     'GO:1901720':'GO:1905560', 
                     'GO:0034291':'GO:0140911', 'GO:0034290':'GO:0140911', 'GO:0034292':'GO:0140911',
                     'GO:0004147':'GO:0043754'}

    def get_subontology(aspect):
        """aspect should be 'BPO', 'CCO', and 'MFO' respetively"""

        obonet_graph = obonet.read_obo('http://purl.obolibrary.org/obo/go/go-basic.obo')

        subontology_roots = {'BPO':'GO:0008150',
                             'CCO':'GO:0005575',
                             'MFO':'GO:0003674'}
        return obonet_graph.subgraph(nx.ancestors(obonet_graph, subontology_roots[aspect]))

    
    df = pd.read_json(annotation_file, lines = True)

    for aspect, subontology in zip(['P','C','F'], ['BPO', 'CCO', 'MFO']):
        print(f'Processing {subontology} annotations')

        db_file = os.path.join(save_location, f'{db_name}_{subontology}.json')
  
        subontology = get_subontology(subontology)

        for gene, annotations in tqdm(df[df.aspect==aspect].groupby('GeneID')):
            gene_terms = set()
            for term in annotations.term.values:
                if term not in subontology:
                    if term in obsolete_lookup:
                        term = obsolete_lookup[term]
                    else:
                        print(f'{term} Term not found')
                        continue
                gene_terms = gene_terms.union(nx.descendants(subontology, term))
            
            with open(db_file, "a") as f:
                f.write(json.dumps({'GeneID': gene, 'terms': list(gene_terms)})+'\n')
        


if __name__ == '__main__':
    """
    Process GOA annotations from Swiss-Prot (or Tremble+Swiss-Prot) and extract annotations with 
    acceptable evidence codes. Specify data location in first input argument and GAF file name in the second argument.
    """

    parser = argparse.ArgumentParser(description='Download Gene Ontology annotations with experimental evidence codes and propogate labels')
    parser.add_argument('--species', '-s', default='all', help='either "human" or "all" to specify which file to download')
    parser.add_argument('--dest', '-d', default='GOAdata', help='Path to save raw and processed files. Will be created if does not exist.')    
    args = parser.parse_args()
    
    # set up save path
    if args.dest is None: 
        data_location=path.join(os.getcwd(),'data')
    else:
        data_location=args.dest

    os.makedirs(data_location, exist_ok=True)

    # get raw annotations
    if args.species=='all':
        goa_file = 'goa_uniprot_gcrp.gaf.gz'
        #goa_file = 'goa_uniprot_all.gaf.gz'
    elif args.species=='human':    
        goa_file = 'goa_human.gaf.gz'

    file_source='http://geneontology.org/gene-associations/'


    # download and decompress file
    annot_file = os.path.join(data_location, goa_file)
    
    if not os.path.exists(annot_file):
        print(f'Downloading GO Annotation File (GAF) file {goa_file} from {file_source}')    
        annot_file = download_gofile(os.path.join(file_source,goa_file), annot_file)
        print(f'Decompressed file location: {annot_file}')
 

    filtered_file=os.path.join(data_location, goa_file.split('.')[0]+'_exp.json')
    
    if not os.path.exists(filtered_file):
        print('GO Annotation File ready. Extracting annotations with experiment evidence codes')
        print(f'Saving to {filtered_file}')
        filter_evidence(annot_file, filtered_file) 


    propogate_terms(filtered_file, data_location, args.species)


