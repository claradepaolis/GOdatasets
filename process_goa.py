import os
import sys
import json   
import GOA
from tqdm import tqdm


def file_length(filename):
    return sum(1 for line in open(filename))

if __name__ == '__main__':

    """
    Process GOA annotations from Swiss-Prot (or Tremble+Swiss-Prot) and extract annotations with 
    acceptable evidence codes. Specify data location in first input argument and GAF file name in the second argument.
    """

    data_location = sys.argv[1]
    annotations_dataset = os.path.join(data_location, sys.argv[2])
    save_location = os.path.join(data_location, 'goa_all_full_exp_evidence.json')

    count = 0
    ok_evidence = ['TAS', 'EXP', 'IC', 'IPI', 'IDA', 'IMP', 'IGI', 'IEP']
    gaf_filelen = file_length(annotations_dataset)

    with open(annotations_dataset) as handle:
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

