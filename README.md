# GOdatasets
Process Gene Ontology annotation data and corresponding protein sequences.

This code will download Gene Ontology annotations and filter them by evidence code and source database (UniProt).  

Once genes and annotations are processed (filtered), annotations will propogated for each gene so that a full labeling in the ontology is created.  

Finally, sequences for the set proteins with annotations will be found in SwissProt (and optionally TrEMBL).


## Usage
### Set up
After cloning this repo locally, set up dependecies (see requirements.txt).
 
To process annotations, you will need the following packages: pandas, networkx, obonet, and tqdm.  
If using conda to install dependencies, obonet can be found in the biobuilds channel (`conda install -c biobuilds obonet`)

If also processing sequences, the Biopython package will be needed. If using conda, find it in the conda-forge channel:  
`conda install -c conda-forge biopython`

### Getting Raw Data
To download protein sequences, `fetch_current_release.sh` will:
1. Create a directory with the current release name, eg 2023_01, and
2. download Swiss-Prot and optionally TrEMBL sequence data (use the `-t` flag to enable. This is a large file that can take a long time to download and process).
3. You can also generate a SQL index file for fast access sequences using UniProt ID (use the -i flag)

The full UniProt ID mapping file can also be downloaded with the -m flag, but it is 20+ GB. Avoid downloading if not needed.  
At this point you can also download the ful set of Gene Ontology annotations with the -g flag. 
This will also download the current GO ontology structure file, go-basic.obo.

If you only need a single species subset of the GO, use `get_raw_annotations.sh` instead
 

### Annotation Data 
To get GO annotations you can run the above when downloading sequences. Alternatively,  

`get_raw_annotations.sh` will download Uniprot GO annotations from [https://ftp.ebi.ac.uk/pub/databases/GO/goa/](https://ftp.ebi.ac.uk/pub/databases/GO/goa/)

`./get_raw_annotations.sh -s human` to download only human proteins.  Species name can also be any of the stand-alone species provided at the link above, eg. chicken, cow, mouse, etc.
`./get_raw_annotations.sh` to get all proteins with annotations. This is the same file that will download when running `fetch_current_release.sh -g`


#### Processing Annotations
Once annotation file has been downloaded, process annotations with `process_goa.py`. This file will:  
1. filter by experimental(EXP, IPI, IDA, IMP, IGI, IEP), inferred (IC), Traceable Author Statement (TAS), and (optionally) high-throughput evidence codes (HTP, HDA, HMP, HGI, HEP),  
2. save records to json file,  
3. propogate labeled annotations according to the onotology graph at [http://purl.obolibrary.org/obo/go/go-basic.obo](http://purl.obolibrary.org/obo/go/go-basic.obo) or local file.

To run, specify the location of the downloaded gaf file with the `--gaf` flag, for example    
`python process_goa.py --gaf raw_data/goa_human.gaf`  

To include high-throughput evidence codes, use the `--high` flag, for example  
`python process_goa.py --gaf raw_data/goa_human.gaf --high`

If you have a local ontology graph file (.obo), specify its path with the `--obo` flag:  
`python process_goa.py --gaf raw_data/goa_human.gaf --obo raw_data/go.obo`

## Background
[The Gene Ontology](http://geneontology.org/docs/ontology-documentation/) provides a set of ontologies, a set of classes (terms) and relations between them, that describes the functions of genes. The ontology is comprised of three subontologies: Biological Process, Cellular Component, and Molecular Function.

### Annotations and Ground Truth 
If a term appears as annotated for a gene, it means that this gene is believed to have this function. By processing these annotated terms, we can generate a dataset of genes and their ground truth labels for each term. The absence of a term annotation does not necessarily mean a gene _does not have this function_, only that this annotation does not exist (yet) in the GO.  

### Evidence
Different types of evidence may be used to determine a term annotation for a gene. More information about evidence can be found here: http://geneontology.org/docs/guide-go-evidence-codes/. In this repo, we filter for terms with experimental, implied, and high-throughput evidence codes (EXP, IPI, IDA, IMP, IGI, IEP, TAS, IC, HTP, HDA, HMP, HGI, HEP), but this can be modified within the file `process_goa.py`.   

### Obsolete Annotated Terms 
In some cases, the GO annotations contain terms that have been obsoleted and do not appear in the OBO graph. We have manually checked all obsolete terms that appear in the filtered set of annotations against the obsoleted comments in [AmiGO](http://amigo.geneontology.org/amigo/). When comments indicate a one-to-one replacement with another term, ie "GO:XXXXXXX replaced by GO:YYYYYYY, we make this substitution in the code and replace the obsolete term with its replacement. In many cases, a direct replacement is not available and we remove the obsolete term from the dataset. The details for each obsolete term are provided in the file `obsolete_terms.tsv`. Only those that include "replaced by" are replaced. All others are ignored. 

### Other available information
The annotation data contains several fields in the GO Annotation File (GAF) 2.2 format. More information on these fields can be found here http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/.

## Alternatives and More info
The GO annotations used here, and others, can be found here: https://www.uniprot.org/downloads.  
You can find more information and READMEs for database files here and its subdirectories: https://ftp.uniprot.org/pub/databases/uniprot/. 

