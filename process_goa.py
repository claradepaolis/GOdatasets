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
    ok_evidence = exp_evidence + h_evidence + inferred_evidence 
    gaf_filelen = file_length(annot_file)

    with open(annot_file) as handle:
        pbar = tqdm(GOA.gafiterator(handle), total=gaf_filelen, unit_scale=True)
        for rec in pbar:
            pbar.set_description("Matches Found {}".format(count))
            if rec['Evidence'] in ok_evidence and rec['DB']=='UniProtKB':
                count += 1
                #rec_info = {'GeneID': rec['DB_Object_ID'], 'term': rec['GO_ID'], 'evidence': rec['Evidence'],
                #            'aspect': rec['Aspect'], 'synonym': rec['Synonym'], 'symbol': rec['DB_Object_Symbol'],
                #            'taxon': rec['Taxon_ID'], 'date': rec['Date']}
                #with open(save_location, 'a', newline='') as f:
                #    json.dump(rec_info, f)
                #    f.write(os.linesep)     
               
                # save full records
                with open(save_location, 'a', newline='') as f:
                    json.dump(rec, f)
                    f.write(os.linesep) 

def propogate_terms(annotation_file, save_location, db_name, obo_file=None):

    # a couple known term replacements
    obsolete_lookup={'GO:0006975':'GO:0042770', 
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
                    
    def get_subontology(aspect, obo):
        """aspect should be 'BPO', 'CCO', and 'MFO' respetively"""
        
        if obo is None:
             obo = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
        obonet_graph = obonet.read_obo(obo)
        # keep only "is_a" and "part_of" edges
        remove_edges = [(i, j, k) for i, j, k in obonet_graph.edges if not(k=="is_a" or k=="part_of")]
        obonet_graph.remove_edges_from(remove_edges)


        subontology_roots = {'BPO':'GO:0008150',
                             'CCO':'GO:0005575',
                             'MFO':'GO:0003674'}
        return obonet_graph.subgraph(nx.ancestors(obonet_graph, subontology_roots[aspect]))

    
    df = pd.read_json(annotation_file, lines = True)

    for aspect, subontology in zip(['P','C','F'], ['BPO', 'CCO', 'MFO']):
        print(f'Processing {subontology} annotations')

        db_file = os.path.join(save_location, f'{db_name}_{subontology}.json')
  
        subontology = get_subontology(subontology, obo_file)

        for gene, annotations in tqdm(df[df.Aspect==aspect].groupby('DB_Object_ID')):
            gene_terms = set()
            for term in annotations.GO_ID.values:
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
    parser.add_argument('--gaf', '-g', default=None, help='Path to GAF file with annotations. If none provided, will download from ebi.ac.uk')
    parser.add_argument('--obo', '-o', default=None, help='Path to OBO graph file if local. If empty (default) current OBO structure at run-time will be downloaded from http://purl.obolibrary.org/obo/go/go-basic.obo')
    args = parser.parse_args()
    
    # set up save path
    if args.dest is None: 
        data_location=path.join(os.getcwd(),'data')
    else:
        data_location=args.dest
    
    os.makedirs(data_location, exist_ok=True)

    # get ontology structure
    if args.obo is not None:
        obo_file = args.obo
    else:
        obo_file = None

    # get raw annotations
    if args.gaf is not None:
        goa_file = args.gaf
    else:
        if args.species=='all':
            goa_file = 'goa_uniprot_all.gaf'
            file_source = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/'

        elif args.species=='human':    
            goa_file = 'goa_human.gaf'
            file_source ='https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/'

        # download and decompress file
        annot_file = os.path.join(data_location, goa_file)
    
        if not os.path.exists(annot_file):
            print(f'Downloading GO Annotation File (GAF) file {goa_file} from {file_source}')    
            annot_file = download_gofile(os.path.join(file_source,goa_file), annot_file)
            print(f'Decompressed file location: {annot_file}')
 
    
    #output file
    filtered_file=os.path.join(data_location, os.path.split(goa_file)[-1].split('.')[0]+'_exp.json')
    
    # save filtered evidence to file
    if not os.path.exists(filtered_file):
        print('GO Annotation File ready. Extracting annotations with experiment evidence codes')
        print(f'Saving to {filtered_file}')
        filter_evidence(goa_file, filtered_file) 


    propogate_terms(filtered_file, data_location, args.species, obo_file)


