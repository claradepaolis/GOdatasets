#!/bin/bash

DATA_DIR="/data/$USER/dataset/raw"  #specify where to save data files
CODE_DIR="." # specify code directory
cd $DATA_DIR

# Sequence data
cd $DATA_DIR
echo "============================================================"
echo "Downloading UniProtKB/Swiss-Prot fasta file: approx 86M"
echo "============================================================"
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gzip -d uniprot_sprot.fasta.gz
echo "============================================================"
echo "Downloading UniProtKB/TrEMBL fasta file: approx 47G"
echo "============================================================"
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gzip -d uniprot_trembl.fasta.gz
echo "============================================================"
echo "Sequence files ready. Running 'python process_trembl.py' to create SQLite lookup table"
echo "============================================================"
cd $CODE_DIR
python create_dataset/process_trembl.py $DATA_DIR
