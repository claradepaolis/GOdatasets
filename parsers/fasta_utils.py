import sqlite3
from Bio import SeqIO, SeqRecord, Seq


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
        for gene in set(proteins):
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
