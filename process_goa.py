import os
import argparse
import pandas as pd
import urllib
from goa_utils import filter_evidence, clean_annotations, propagate_terms


def download_file(source_path, save_path):
    if not os.path.exists(save_path):
        urllib.request.urlretrieve(source_path, filename=os.path.join(save_path))
    
    out_file = save_path.split('.gz')[0]
    if not os.path.exists(out_file):
        with gzip.open(save_path, 'rb') as f_in:
            with open(out_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    return out_file



if __name__ == '__main__':
    """
    Process GOA annotations from UniProt extract annotations with acceptable evidence codes.
    Specify data location in first input argument and GAF file name in the second argument.
    """

    parser = argparse.ArgumentParser(description='Download Gene Ontology annotations with experimental evidence codes and propogate labels')
    parser.add_argument('--dest', '-d', default=None,
                        help='Path to save processed files. Will be created if does not exist.')
    parser.add_argument('--gaf', '-g',
                        help='Path to GAF file with annotations.')
    parser.add_argument('--obo', '-o', default=None,
                        help='Path to OBO graph file if local. If empty (default) current OBO structure at run-time will be downloaded from http://purl.obolibrary.org/obo/go/go-basic.obo')
    parser.add_argument('--high', action='store_true',
                        help='Flag to include high-throughput evidence codes')    
    args = parser.parse_args() 
    
    # get raw annotations
    if args.gaf is not None:
        goa_file = args.gaf
        assert os.path.exists(goa_file), "GAF file does not exist"

    # get destination location for output files
    if args.dest is not None:
        data_location = args.dest
    else:
        data_location = os.path.split(goa_file)[0]
    
    # get ontology structure
    if args.obo is not None:
        obo_file = args.obo
    else:
        print('Downloading OBO file from http://purl.obolibrary.org/obo/go/go-basic.obo')
        obo_file = download_file('http://purl.obolibrary.org/obo/go/go-basic.obo', 
                                  os.path.join(data_location, 'go-basic.obo'))
    
    # high-throughput evidence codes
    include_highthr = args.high

    # output file to save
    filtered_file=os.path.join(data_location, os.path.split(goa_file)[-1].split('.')[0]+'_evidence.json')

    # save filtered evidence to file
    if not os.path.exists(filtered_file):
        print('Extracting annotations with experiment evidence codes')
        print(f'Saving to {filtered_file}')
        filtered_annotations = filter_evidence(goa_file, filtered_file, highthr=include_highthr)
    else:
        print(f'Filtered evidence GO annotation file exists. Loading annotations from file {filtered_file}')
        filtered_annotations = pd.read_json(filtered_file, lines=True)

    # Remove any duplicates or negative labels 
    print('Removing duplicates and negative labels')
    annotations_df = clean_annotations(filtered_annotations)

    # Propagate labels to root and save to file
    print('Propagating annotations to ontology roots and removing obsolete labels')
    terms_df = propagate_terms(annotations_df, obo_file)
    # In some cases, only obsolete terms are annotated for some protein
    # so we get rid of any proteins without terms
    annotations_df = annotations_df[annotations_df.DB_Object_ID.isin(terms_df.EntryID)]

    # expand list of terms to that each one is in one row
    all_terms = terms_df.explode('term').reset_index(drop=True)
    terms_file = os.path.join(data_location, 'terms.tsv')
    print(f'Saving terms to file {terms_file}')
    all_terms.to_csv(terms_file, index=False, sep='\t')


