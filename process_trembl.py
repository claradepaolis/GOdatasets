import os
import sys
import time
from Bio import SeqIO

if __name__ == '__main__':
    """
    Create SQL db for fast access to TrEMBL sequence data.
    Specify data location in first input argument 
    """

    data_location = sys.argv[1]
    save_location = sys.argv[2]
    save_file = os.path.join(save_location, 'uniprot_trembl.idx')
    source_file = os.path.join(data_location, 'uniprot_trembl.fasta')
    
    print("Creating SQLite index file for fasta file {}. Will be saved to {}".format(source_file, save_file))
    start_time = time.time()

    def get_accession(name):
        parts = name.split("|")
        return parts[1]
    indices_sql = SeqIO.index_db(save_file, source_file, "fasta", 
                                 key_function=get_accession)


    #indices = SeqIO.index_db(save_file,
    #                         os.path.join(data_location, 'uniprot_trembl.fasta'), "fasta")
    duration = time.time()-start_time

    print('Elapsed time: {}s'.format(duration))
