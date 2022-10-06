# GOdatasets
Process Gene Ontology annotation data


## Usage
### Set up
After cloning this repo locally, set up dependecies:  
`conda install tqdm`  
We also need support for the GAF 2.2 file format, not supported in the current 
release of [Biopython (1.79)](https://github.com/biopython/biopython).  
Instead, we just copy the file that contains GAF parsing code from the commit that added support: 

`wget https://raw.githubusercontent.com/biopython/biopython/ba5dfd472862c9efe797cdd6d5fe011e8cf96f0e/Bio/UniProt/GOA.py`


### Getting Annotation Data
To get GO annotations, `get_raw_annotations.sh` will:
1. download GO Annotations from `ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`,
2. filter by experimental evidence codes,
3. save records to json file.
Run time is approx 1 hour on my server. The resulting json file is approximately 200MB

### Getting Sequences
To download protein sequences, `get_raw_sequences.sh` will:
1. download Swiss-Prot and TrEMBL sequence data
2. generate a SQL database for fast access to TrEMBL sequences
Run time is approx 3.5 hours on my server.

## Background
[The Gene Ontology](http://geneontology.org/docs/ontology-documentation/) provides a set of ontologies, a set of classes (terms) and relations between them, that describes the functions of genes. The ontology is comprised of three subontologies: Biological Process, Cellular Component, and Molecular Function.

### Annotations and Ground Truth 
If a term appears as annotated for a gene, it means that this gene is believed to have this function. By processing these annotated terms, we can generate a dataset of genes and their ground truth labels for each term. The absence of a term annotation does not necessarily mean a gene _does not have this function_, only that this annotation does not exist (yet) in the GO.  

### Evidence
Different types of evidence may be used to determine a term annotation for a gene. More information about evidence can be found here: http://geneontology.org/docs/guide-go-evidence-codes/. In this repo, we filter only terms with experimental evidence codes, but this can be modified within the file `process_goa.py`.   

### Other available information
The annotation data contains several fields in the GO Annotation File (GAF) 2.2 format. More information on these fields can be found here http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/.

## Alternatives and More info
The GO annotations used here, and others, can be found here: https://www.uniprot.org/downloads.  
You can find more information and READMEs for database files here and its subdirectories: https://ftp.uniprot.org/pub/databases/uniprot/. 
