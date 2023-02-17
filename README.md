# GOdatasets
Process Gene Ontology annotation data

This code will download Gene Ontology annotations and filter them by evidence code and source database.  
Once genes and annotations are processed (filtered), annotations will propogated for each gene so that a full labeling in the ontology is created.  
Finally, sequences for the set proteins with annotations will be found in SwissProt (and optionally TrEMBL).


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

If also processing sequences, the complete Biopython package will be needed:
`conda install -c conda-forge biopython`

### Annotation Data
To get GO annotations,  

`get_raw_annotations.sh` will download Uniprot GO annotations from [https://ftp.ebi.ac.uk/pub/databases/GO/goa/](https://ftp.ebi.ac.uk/pub/databases/GO/goa/)

`./get_raw_annotations.sh -s human` to download only human proteins.  Species name can also be any of the stand-alone species provided at the link above, eg. chicken, cow, mouse, etc.
`./get_raw_annotations.sh` to get all proteins with annotations.  


#### Processing Annotations
Once annotation file has been downloaded, process annotations with `process_goa.py`. This file will:  
1. filter by experimental, implied, and high-throughput evidence codes (TAS, EXP, IC, IPI, IDA, IMP, IGI, IEP),  
2. save records to json file.  
3. propogate labeled annotations according to the onotology graph at [http://purl.obolibrary.org/obo/go/go-basic.obo](http://purl.obolibrary.org/obo/go/go-basic.obo) 

To run, specify the location of the downloaded gaf file with the `--gaf` flag, for example    
`python process_goa.py --gaf raw_data/goa_human.gaf`


### Getting Sequences
To download protein sequences, `get_raw_sequences.sh` will:
1. download Swiss-Prot and optionally TrEMBL sequence data (use the `-t` flag to enable. This is a large file that can take a long time to download and process)
2. If TrEMBL data was downloaded, this will also generate a SQL database for fast access to TrEMBL sequences

## Background
[The Gene Ontology](http://geneontology.org/docs/ontology-documentation/) provides a set of ontologies, a set of classes (terms) and relations between them, that describes the functions of genes. The ontology is comprised of three subontologies: Biological Process, Cellular Component, and Molecular Function.

### Annotations and Ground Truth 
If a term appears as annotated for a gene, it means that this gene is believed to have this function. By processing these annotated terms, we can generate a dataset of genes and their ground truth labels for each term. The absence of a term annotation does not necessarily mean a gene _does not have this function_, only that this annotation does not exist (yet) in the GO.  

### Evidence
Different types of evidence may be used to determine a term annotation for a gene. More information about evidence can be found here: http://geneontology.org/docs/guide-go-evidence-codes/. In this repo, we filter for terms with experimental, implied, and high-throughput evidence codes ('EXP', 'IPI', 'IDA', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'), but this can be modified within the file `process_goa.py`.   

### Obsolete Annotated Terms 
In some cases, the GO annotations contain terms that have been obsoleted and do not appear in the OBO graph. We have manually checked all obsolete terms that appear in the filtered set of annotations against the obsoleted comments in [AmiGO](http://amigo.geneontology.org/amigo/). When comments indicate a one-to-one replacement with another term, ie "GO:XXXXXXX replaced by GO:YYYYYYY, we make this substitution in the code and replace the obsolete term with its replacement. In many cases, a direct replacement is not available and we remove the obsolete term from the dataset. The details for each obsolete term are provided in the file `obsolete_terms.tsv`. Only those that include "replaced by" are replaced. All others are ignored. 

### Other available information
The annotation data contains several fields in the GO Annotation File (GAF) 2.2 format. More information on these fields can be found here http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/.

## Alternatives and More info
The GO annotations used here, and others, can be found here: https://www.uniprot.org/downloads.  
You can find more information and READMEs for database files here and its subdirectories: https://ftp.uniprot.org/pub/databases/uniprot/. 

