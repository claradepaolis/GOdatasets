import os
import sys
import json   
import time
from tqdm import tqdm
import urllib.request
import gzip
import shutil
import argparse

import pandas as pd
import networkx as nx
import obonet
from Bio import SeqIO
import GOA
from parsers.fasta_utils import read_fasta_sql
from parsers.goa_utils import filter_evidence, clean_annotations, propagate_terms


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


def get_all_taxonomies(annotations_df):

    # get species for each protein
    taxonomy_df = annotations_df[['DB_Object_ID', 'species']].drop_duplicates()
    taxonomy_df = taxonomy_df.rename({'DB_Object_ID':'EntryID', 'species':'taxonomyID'}, axis=1)
    # make sure there aren't repeated protein entries and all proteins are in list
    assert len(taxonomy_df)==len(taxonomy_df.EntryID.unique()), 'Some taxonomy ID is duplicated'
  
    return taxonomy_df


def find_sequences(terms_df, swissprot, trembl, trembl_index):
    # load sequences
    record_dict = SeqIO.to_dict(SeqIO.parse(swissprot, "fasta"))

    # SwissProt FASTA file uses keys with format "sp|uniprot_acession_id|gene_name". 
    # We want to use on the accession id, so we map the orignal keys to this id

    # first make sure each UniProt accession ID maps to one key in the record dictionary
    assert len(record_dict)==len(set(orig_key.split('|')[1] for orig_key in record_dict.keys()))

    accession_id_map = {orig_key.split('|')[1]:orig_key for orig_key in record_dict.keys()}

    # find all proteins that are in the SwissProt DB
    proteins_with_terms = set(terms_df.EntryID.unique().tolist())
    proteins_in_sp = set(accession_id_map.keys()).intersection(proteins_with_terms)

    # Get sequences for all proteins in dataset
    # Set ID in FASTA records to UniProt accession ID
    # First get sequences that are in SwissProt
    for accession in proteins_in_sp:
        record_dict[accession_id_map[accession]].id = accession

    records = [record_dict[accession_id_map[protein]] for protein in proteins_in_sp] 

    if len(proteins_with_terms.difference(proteins_in_sp)) > 0:
        print('Processing trEMBL sequence file.')
        if os.path.exists(trembl_index):
            tr_records = read_fasta_sql(trembl_index, trembl, set(proteins_with_terms).difference(proteins_in_sp))
        else:
            print('Processing may take about 30 minues if TrEMBL file has not been previously indexed') 
            start_time = time.time()
            def get_accession(name):
                parts = name.split("|")
                return parts[1]
            indices_sql = SeqIO.index_db(trembl_index, trembl, "fasta",
                                         key_function=get_accession)
            tr_records = []
            for accession in set(proteins_with_terms).difference(proteins_in_sp):
                seqrecord = indices_sql[accession]
                seqrecord.id = accession  # change the ID to just the accession number
                tr_records.append(seqrecord)

            print(f'Time elapsed: {(time.time()-start_time)/60}')

        records += tr_records

    return records
  


if __name__ == '__main__':
    """
    Process GOA annotations from Swiss-Prot (or Tremble+Swiss-Prot) and extract annotations with 
    acceptable evidence codes. Specify data location in first input argument and GAF file name in the second argument.
    """

    parser = argparse.ArgumentParser(
        description='Download Gene Ontology annotations with experimental evidence codes and propogate labels')
    parser.add_argument('--raw', default='.', 
                        help='Path to raw files for annotations, OBO, and fasta files')
    parser.add_argument('--dest', '-d', 
                        help='Path to save processed files. Will be created if does not exist.')    
    parser.add_argument('--gaf', '-g', default=None, 
                        help='Path to GAF file with annotations. If none provided, will download from ebi.ac.uk')
    parser.add_argument('--obo', '-o', default=None, 
                        help='Path to OBO graph file if local. If empty (default) current OBO structure at run-time will be downloaded from http://purl.obolibrary.org/obo/go/go-basic.obo')
    parser.add_argument('--swiss', '-s',
                        help='Path to swissprot fasta file')
    parser.add_argument('--trembl', '-t',
                        help='Path to trembl fasta file. Index file should have the name and in the same path but with .idx file extension')

    args = parser.parse_args()
    
    # set up save path
    data_location = args.raw
 
    if args.dest is None: 
        save_location=path.join(os.getcwd(),'data')
    else:
        save_location=args.dest
    
    os.makedirs(data_location, exist_ok=True) # create source dir in case we need to download any raw files
    os.makedirs(save_location, exist_ok=True)

    # get ontology structure
    if args.obo is not None:
        obo_file = args.obo
    else:
        obo_file = None

    # get raw annotations
    if args.gaf is not None:
        goa_file = args.gaf
    else:
        goa_file = 'goa_uniprot_all.gaf'
        file_source = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/'

        # download and decompress file
        annot_file = os.path.join(data_location, goa_file)
    
        if not os.path.exists(annot_file):
            print(f'Downloading GO Annotation File (GAF) file {goa_file} from {file_source}')    
            annot_file = download_gofile(os.path.join(file_source,goa_file), annot_file)
            print(f'Decompressed file location: {annot_file}')
 
    swissprot_file = args.swiss
    trembl_file = args.trembl
    trembl_index_file = args.trembl.split('.fasta')[0]+'.idx'

    # output file to save
    filtered_file=os.path.join(data_location, os.path.split(goa_file)[-1].split('.')[0]+'_evidence.json')
    
    # save filtered evidence to file
    if not os.path.exists(filtered_file):
        print('Extracting annotations with experiment evidence codes')
        print(f'Saving to {filtered_file}')
        filtered_annotations = filter_evidence(goa_file, filtered_file) 
    else:
        print(f'Filtered evidence GO annotation file exists. Loading annotations from file {filtered_file}')
        filtered_annotations = pd.read_json(filtered_file, lines=True)

    # Remove any duplicates or negative labels 
    print('Removing duplicates and negative labels')
    annotations_df = clean_annotations(filtered_annotations)

    # Propagate labels to root and save to file
    print('Propagating annotations to ontology roots')
    terms_df = propagate_terms(annotations_df, obo_file) 
    # In some cases, only obsolete terms are annotated for some protein
    # so we get rid of any proteins without terms
    annotations_df = annotations_df[annotations_df.DB_Object_ID.isin(terms_df.EntryID)]
    # expand list of terms to that each one is in one row
    all_terms = terms_df.explode('term').reset_index(drop=True) 
    terms_file = os.path.join(save_location, 'terms.tsv')
    print(f'Saving terms to file {terms_file}')
    all_terms.to_csv(terms_file, index=False, sep='\t') 
    

    # Find taxonomies and save to file
    taxdf = get_all_taxonomies(annotations_df)
    taxonomy_file = os.path.join(save_location, 'taxonomy.tsv')
    print(f'Saving terms to file {taxonomy_file}')
    taxdf.to_csv(taxonomy_file, index=False, sep='\t') 

    
    # Find sequences and save to file
    print('Finding sequences for annotated proteins')
    records = find_sequences(all_terms, swissprot=swissprot_file, trembl=trembl_file, trembl_index=trembl_index_file)
    seq_file_path = os.path.join(save_location,'sequences.fasta')
    print('Saving sequences to file {}'.format(seq_file_path))
    for r in records:
        with open(seq_file_path, 'a') as handle:
            SeqIO.write(r, handle, 'fasta')


