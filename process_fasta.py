import time
import sqlite3
from tqdm import tqdm
from Bio import SeqIO, SeqRecord, Seq
import argparse


def create_sql_db(input_file):
    
    save_file = input_file + '.idx'
    
    print("Creating SQLite index file for fasta file {}. Will be saved to {}".format(input_file, save_file))
    start_time = time.time()

    def get_accession(name):
        parts = name.split("|")
        return parts[1]
    
    indices_sql = SeqIO.index_db(save_file, input_file, 'fasta', key_function=get_accession)
    
    duration = time.time()-start_time

    print(f'{duration:.2f} s elapsed')
    return indices_sql


def read_fasta_sql(index_file, sequence_file, proteins):
    """
    Find a set of proteins in a fasta sequence file that has been indexed into a
    sqlite file by biopython's Bio.SeqIO.index_db method
    :param index_file: SQLite index file for sequence_file
    :param sequence_file: fasta format file with sequences indexed by index_file
    :param proteins: list or set of proteins to be matched. These must match the key function for the database
    :return: list of matching records. List items are a biopython SeqRecord object
    """
    con = sqlite3.connect(index_file)

    matching_records = [] 
    with open(sequence_file, 'r') as filehandle:
        for gene in tqdm(set(proteins)):
            sql = f'select * from offset_data where key = "{gene}"'
            cursor = con.execute(sql)
            record = cursor.fetchone()
            if record is not None:
                offset = record[-2]
                length = record[-1]
                filehandle.seek(offset, 0)
                seq_record = filehandle.read(length)
                sequence = ''.join(seq_record.split('\n')[1:-1])
                description = seq_record.split('\n')[0].split('>')[1]
                rec=SeqRecord.SeqRecord(Seq.Seq(sequence), description=description, 
                                        id=gene, name=description.split(' ')[0])
                matching_records.append(rec)

    return matching_records


if __name__ == '__main__':
    """
    Create SQL index/db of Uniprot fasta protein sequence files.
    """
    
    parser = argparse.ArgumentParser(
        description='Create or read SQL index/db of Uniprot fasta protein sequence files.')
    parser.add_argument('--fasta', help='Path to raw fasta file')
    args = parser.parse_args()    

    create_sql_db(args.fasta) 
