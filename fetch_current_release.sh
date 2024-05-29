#!/bin/bash

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Download UniProt and Gene Ontology dataset raw files. 
         Protein sequence databases in fasta format from Swiss-Prot and (optionally with -t flag) TrEMBL will be downloaded.
         The Gene Ontology annotations can also be downloaded with the -g flag. (About 16 GB)
         The full UniProt ID mapping file can be downloaded with the -m flag, but it is 20+ GB. Avoid downloading if not needed. 

         It's suggested to create a SQLite lookup table (with the -i flag) to allow for fast indexing into the fasta files by UniProt ID.
         All files will be saved in a directory ./YYYY-XX to match the release date and number (XX) of the current release.
         The directory will be created if it does not already exist"
   echo
   echo "Syntax: get_raw_sequences [-t|m|i|g]"
   echo "options:"
   echo "-h     Print this Help."
   echo "-t     Flag to download TrEMBL sequence database (large)."
   echo "-m     Flag to download ID mapping file"
   echo "-i     Flag to create index file for fast UniProt ID lookups"
   echo "-g     Flag to download GO annotations"
   echo
}
############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
dl_TrEMBL=false
dl_mapping=false
create_index=false
get_goa=false
while getopts ":htmig:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      t) # TrEMBL option
         dl_TrEMBL=true;;
      m) # ID mapping file option
         dl_mapping=true;;
      i) 
         create_index=true;;
      g)
         get_goa=true;;
      \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

# Get current release number
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt
relname=$(cat relnotes.txt | head -n 1 | cut -d " " -f 3)
echo "Downloading UniProt Release " $relname

DATA_DIR="./"$relname
mkdir -p $DATA_DIR
mv relnotes.txt $DATA_DIR
cd $DATA_DIR

# Fetch support files: ID mapping, READMEs
if $dl_mapping; then
 if [ ! -f idmapping.dat.gz ]; then
  wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
 fi
 if [ ! -f idmapping_README ]; then
  wget -O idmapping_README ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
 fi
fi

if [ ! -f README ]; then
 wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/README
fi

if [ ! -f reldate.txt ]; then
 wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt
fi

# Sequence data
echo "create index" $create_index
echo "trembl" $dl_TrEMBL
echo "get goa" $get_goa

if [[ -f "uniprot_sprot.fasta.gz" ]] || [[ -f "uniprot_sprot.fasta" ]]; then 
 echo "release" $relname "file uniprot_sprot.fasta already exists"
else
 echo "============================================================"
 echo "Downloading UniProtKB/Swiss-Prot fasta file: approx 86M"
 echo "============================================================"
 wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
fi

if $dl_TrEMBL; then
 if [[ -f "uniprot_trembl.fasta.gz" ]] || [[ -f "uniprot_trembl.fasta" ]]; then
  echo "release" $relname "file uniprot_trembl.fasta already exists"
 else
  echo "============================================================"
  echo "Downloading UniProtKB/TrEMBL fasta file: approx 47G"
  echo "============================================================"
  wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
 fi
fi

if $create_index; then
 if [[ -f "uniprot_sprot.fasta.gz" ]] || [[ -f "uniprot_sprot.fasta" ]]; then
  if [[ -f "uniprot_sprot.fasta.idx" ]]; then
   echo "============================================================"
   echo "UniProtKB/Swiss-Prot index file uniprot_sprot.fasta.idx already exists"
   echo "============================================================"
  else 
   echo "============================================================"
   echo "UniProtKB/Swiss-Prot fasta file downloaded. Processing to create SQLite lookup table"
   echo "============================================================"
   gzip -d uniprot_sprot.fasta.gz
   python ../process_fasta.py --fasta uniprot_sprot.fasta 
  fi
 fi

 if [[ -f "uniprot_trembl.fasta.gz" ]] || [[ -f "uniprot_trembl.fasta" ]]; then
  if [[ -f "uniprot_trembl.fasta.idx" ]]; then
   echo "============================================================"
   echo "UniProtKB/TrEMBL index file uniprot_trembl.fasta.idx already exists"
   echo "============================================================"
  else
   echo "============================================================"
   echo "UniProtKB/TrEMBL fasta file downloaded. Processing to create SQLite lookup table"
   echo "============================================================"
   gzip -d uniprot_trembl.fasta.gz
   python ../process_fasta.py --fasta uniprot_trembl.fasta
  fi
 fi
fi


# GO annotations
if $get_goa; then
 wget -O goa_release_numbers.txt https://ftp.ebi.ac.uk/pub/databases/GO/goa/current_release_numbers.txt 
 wget https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz

 # Get the release date of the ontology and use it to name the OBO files 
 wget -O go_reldate.json  http://current.geneontology.org/metadata/release-date.json
 godate="$(cat go_reldate.json | grep -o ....-..-..)"
  
 wget -O go-basic-"$godate".obo http://current.geneontology.org/ontology/go-basic.obo
 wget -O go-"$godate".obo http://current.geneontology.org/ontology/go.obo
fi

