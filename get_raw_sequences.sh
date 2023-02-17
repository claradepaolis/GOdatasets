#!/bin/bash

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Download decompress and protein sequence databases from Swiss-Prot and (optionally) TrEMBL.
         If TrEMBL database is downloaded, a SQLite lookup table will be created to allow for fast look-ups by (Uniprot) Gene ID.
         All files will be saved in a directory ./GOdataset/raw, which will be created if it does not exist"
   echo
   echo "Syntax: get_raw_sequences [-t]"
   echo "options:"
   echo "-h     Print this Help."
   echo "-t     Flag to download TrEMBL sequence database (large)."
   echo
}
############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
dl_TrEMBL=false
while getopts ":ht" option; do
   case $option in
      h) # display Help
         Help
         exit;;
     t) # TrEMBL option
         dl_TrEMBL=true;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done


DATA_DIR="./GOdataset/raw"  #specify where to save data files
CODE_DIR=$(pwd) # specify code directory
mkdir -p $DATA_DIR

# Sequence data
cd $DATA_DIR
echo "============================================================"
echo "Downloading UniProtKB/Swiss-Prot fasta file: approx 86M"
echo "============================================================"
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gzip -d uniprot_sprot.fasta.gz

if $dl_TrEMBL; then
 echo "============================================================"
 echo "Downloading UniProtKB/TrEMBL fasta file: approx 47G"
 echo "============================================================"
 wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
 gzip -d uniprot_trembl.fasta.gz
 echo "============================================================"
 echo "Sequence files ready. Running 'python process_trembl.py' to create SQLite lookup table"
 echo "============================================================"
 cd $CODE_DIR
 python process_trembl.py $DATA_DIR
fi
