import os
import obonet


if __name__ == "__main__":

    obonet_graph = obonet.read_obo('http://purl.obolibrary.org/obo/go/go-basic.obo')
    # ontology roots
    bpo_root ='GO:0008150' 
    cco_root ='GO:0005575'
    mfo_root ='GO:0003674'
    
    raw_location = os.path.join('/data/clara/GOdatasets/GOdataset', 'raw')
    
    annotation_file = os.path.join(raw_location, 'goa_fullexpevidence.json')
    annotations_df, all_labeled_genes = annotation_data(annotation_file, obonet_graph)

    
