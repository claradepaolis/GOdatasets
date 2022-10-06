#!/bin/bash

DATA_DIR="/data/$USER/dataset/raw"  # specify where to save data
CODE_DIR="."  # specify code location
cd $DATA_DIR

# GO Annotation File (GAF) file (GCRP)
echo "============================================================"
echo "Downloading GO Annotation File (GAF) file (GCRP subset): approx 2.7GB"
echo "============================================================"
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_gcrp.gaf.gz
gzip -d goa_uniprot_gcrp.gaf.gz
echo "============================================================"
echo "GO Annotation File ready. Running 'python process_goa.py' to extract experiment evidence codes"
echo "============================================================"
cd $CODE_DIR
python process_goa.py $DATA_DIR goa_uniprot_gcrp.gaf
