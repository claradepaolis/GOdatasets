# GOdatasets
Process Gene Ontology annotation data


## Usage
### Set up
After cloning this repo locally, set up dependecies:  
```
conda install tqdm
conda install pandas
conda install -c biobuilds obonet
```  

We also need support for the GAF 2.2 file format, not supported in the current 
release of [Biopython (1.79)](https://github.com/biopython/biopython).  
Instead, we just copy the file that contains GAF parsing code from the commit that added support: 

`wget https://raw.githubusercontent.com/biopython/biopython/ba5dfd472862c9efe797cdd6d5fe011e8cf96f0e/Bio/UniProt/GOA.py`


### Getting Annotation Data
To get GO annotations,  
`python process_goa.py --species human` or `python process_goa.py --species all`

This will  
1. download GO Annotations from [http://geneontology.org/gene-associations](http://geneontology.org/gene-associations),  
2. filter by experimental evidence codes (TAS, EXP, IC, IPI, IDA, IMP, IGI, IEP),  
3. save records to json file.  
4. propogate labeled annotations according to the onotology graph at [http://purl.obolibrary.org/obo/go/go-basic.obo](http://purl.obolibrary.org/obo/go/go-basic.obo) 
5. save propogated annotaions to files, one for each subontology (BPO, CCO, and MFO) in json format


This code creates a directory `GOAdata` under the source directory and will save all files there.
To specify a different save location, use the `--dest` flag. 


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
