# ensembl ncbi mapping
Improve multiple IDs mapping to one Ensembl gene ID for mygene.info website

Functions are run in order as shown below using three input files.
Mygene-py client is required for mygene.info lookups. 

## Functions in script
```python
multi_mapping_dict = find_multiple_mappings_from_entrezgene_file("gene_ensembl__xref_entrezgene__dm.txt")
ensembl_dict = create_ensembl_gene_id_dict("gene_ensembl__gene__main.txt", multi_mapping_dict)
ensembl_dict = find_ncbi_ids_from_gene2ensembl(ensembl_dict, "gene2ensembl")
write_mapping_ids_to_file(ensembl_dict)
```
