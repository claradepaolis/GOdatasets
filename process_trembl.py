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
    save_file = 'uniprot_trembl.idx'
    print("Creating SQLite index file. Will be saved to {}/{}".format(data_location, save_file))
    start_time = time.time()
    indices = SeqIO.index_db(os.path.join(data_location, save_file),
                             os.path.join(data_location, 'uniprot_trembl.fasta'), "fasta")
    duration = time.time()-start_time

    print('Elapsed time: {}s'.format(duration))
