#!/bin/bash
# M. Clara De Paolis Kaluza, 2022. GNU GPL3.0 license
# Help and argument handling code inspired by https://www.redhat.com/sysadmin/arguments-options-bash-scripts

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Download decompress and Gene Ontology Uniprot annotations (in GAF format) 
         File is downloaded from ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/.
         By default, the complete set goa_uniprot_all.gaf will be downloaded. To download annotatation for only a single species. use the -s argument."
   echo
   echo "Syntax: get_raw_sequences [-h|-s]"
   echo "options:"
   echo "-h     Print this Help."
   echo "-s     Species name if downloading annotations for only a single species"
   echo
}
############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
species="all"
while getopts ":hs:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      s) # Species option
         species="$OPTARG";;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done
echo "Species: $species"
DATA_DIR="./raw_data"  # specify where to save data
CODE_DIR=$(pwd)  # specify code location
mkdir -p $DATA_DIR

species_file=goa_${species}.gaf.gz
all_file=goa_uniprot_all.gaf.gz
if [ "$species" = "all" ]; then
 gaf_file="UNIPROT/"${all_file}
else
 prefix=`echo ${species} | tr a-z A-Z`
 gaf_file=${prefix}"/"${species_file}
fi

# GO Annotation File (GAF) file (GCRP)
echo "============================================================"
echo "Downloading GO Annotation File (GAF) file (GCRP subset): approx 2.7GB"
echo "============================================================"
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/$gaf_file -P $DATA_DIR

cd $DATA_DIR
basefile=$(basename ${gaf_file} ".gz")
gzip -d $basefile".gz"
echo "============================================================"
echo "GO Annotation File ready." 
echo "To filter annotations by experiment evidence codes and propagate labels, run"
echo "python process_goa.py --gaf $DATA_DIR/$basefile"
echo "============================================================"
