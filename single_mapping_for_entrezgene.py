
import ast
import collections

import mygene

from collections import defaultdict
from collections import Counter


def find_multiple_mappings_from_entrezgene_file(gene_ensembl_entrezgene_dm_file):
    """Input gene_ensembl_entrezgene_dm_file, and identify how many NCBI gene IDs there are for
    each ensembl gene ID. Lines in input file are:

    'gene_ensembl__xref_entrezgene__dm.txt' (useful columns in input_file):

    col1: Ensembl gene ID
    col2: NCBI gene ID

    If there is > 1 NCBI gene ID, we need to process further.
    """
    ensembl_dict_from_entrez = defaultdict(list)

    with open(gene_ensembl_entrezgene_dm_file) as file_in:
        next(file_in)
        for line in file_in:
            split_line = line.split("\t")
            ensembl_gene_id_from_entrez = split_line[1].strip()
            ncbi_gene_id_from_entrez = split_line[2].strip()
            ensembl_dict_from_entrez[ensembl_gene_id_from_entrez].append(ncbi_gene_id_from_entrez)

    multi_mapping_dict = {}
    for key in ensembl_dict_from_entrez:
        if len(ensembl_dict_from_entrez[key]) > 1:
            multi_mapping_dict[key] = ensembl_dict_from_entrez[key]

    return multi_mapping_dict

def create_ensembl_gene_id_dict(gene_ensembl_main_file, multi_mapping_dict):
    """Using gene_ensembl_main_file, identify correct ensembl symbol for each
    ensembl gene ID. Add this information to a new dictionary.

    'gene_ensembl__gene__main.txt' (useful columns in input file):

    col1: Ensembl gene ID
    col2: Ensembl symbol
    """
    ensembl_dict = defaultdict(list)
    symbol_file_list = []
    with open(gene_ensembl_main_file) as file_in:
        next(file_in)
        for line in file_in:
            split_line = line.split("\t")
            ensembl_id_from_main = str(split_line[1].strip())
            ensembl_symbol_from_main = split_line[2].strip()

            if ensembl_id_from_main in multi_mapping_dict:
                ensembl_id_dict = {}
                ensembl_id_dict['data'] = {'ncbi_list': multi_mapping_dict[ensembl_id_from_main],
                                           'symbol': ensembl_symbol_from_main,
                                           'gene2ensembl': [] }
                ensembl_dict[ensembl_id_from_main] = ensembl_id_dict

    return ensembl_dict


def find_ncbi_ids_from_gene2ensembl(ensembl_dict, gene2ensembl_file):
    """Input is gene2ensembl_file; maps one NCBI gene ID to one Ensembl gene ID.

    'gene2ensembl' (useful columns in input file):

    col1: NCBI gene ID
    col2: Ensembl gene ID
    """
    with open(gene2ensembl_file) as file_in:
        next(file_in)

        for line in file_in:
            split_line = line.split("\t")
            ensembl_gene_id_from_gene2ensembl = split_line[2].strip()
            ncbi_gene_id_from_gene2ensembl = split_line[1].strip()

            if ensembl_gene_id_from_gene2ensembl in ensembl_dict:
                ensembl_dict[ensembl_gene_id_from_gene2ensembl]['data']['gene2ensembl'].append(ncbi_gene_id_from_gene2ensembl)

    return ensembl_dict


def write_mapping_ids_to_file(ensembl_dict):
    """First use gene2ensembl as single match NCBI gene ID (if == 1 match).
    Next, if no gene2ensembl match, then look at mygene.info to find which NCBI
    ID from the NCBI multi mapping list returns the same ensembl symbol as the
    ensembl main file, and use corresponding NCBI gene ID as single match.
    """
    final_mapping_file = open("final_mapping_file.txt", "w")
    ncbi_list_for_mygene_querymany = []
    ncbi_list_for_mygene_querymany2 = []

    for key in ensembl_dict:
        ncbi_list = ensembl_dict[key]['data']['ncbi_list']
        ensembl_symbol = ensembl_dict[key]['data']['symbol'].upper()
        gene2ensembl_ncbi_gene_id_match_list = ensembl_dict[key]['data']['gene2ensembl']
        ncbi_list_for_mygene_querymany.append(ncbi_list)
        if len(gene2ensembl_ncbi_gene_id_match_list) == 1:
            final_mapping_file.write(key + '\t')
            final_mapping_file.write(gene2ensembl_ncbi_gene_id_match_list[0] + '\n')
        else:
            # only append list if need to query mygene.info
            ncbi_list_for_mygene_querymany.append(ncbi_list)

    ncbi_list_for_mygene_querymany = [item for sublist in ncbi_list_for_mygene_querymany for item in sublist]

    # no need to have duplicates because mygene.info returns one symbol
    ncbi_list_for_mygene_querymany = list(set(ncbi_list_for_mygene_querymany))

    mg = mygene.MyGeneInfo()
    ensembl_symbol_list_from_mygene = mg.querymany(ncbi_list_for_mygene_querymany, scopes='entrezgene', species="all", fields="symbol")

    mygene_website_dict = {}
    for dic in ensembl_symbol_list_from_mygene:
        try:
            mygene_website_dict[dic['query']] = dic['symbol']
        except:
            #if there is no symbol for the query, just don't add it to mygene_website_dict
            pass

    for key in ensembl_dict:
        ncbi_list = ensembl_dict[key]['data']['ncbi_list']
        ensembl_symbol = ensembl_dict[key]['data']['symbol'].upper()
        gene2ensembl_ncbi_gene_id_match_list = ensembl_dict[key]['data']['gene2ensembl']
        if len(gene2ensembl_ncbi_gene_id_match_list) != 1:

            ensembl_symbol_list_from_mygene = []
            for ncbi_id in ncbi_list:
                try:
                    ensembl_symbol_list_from_mygene.append(mygene_website_dict[ncbi_id].upper())
                except:
                    # need this here; keeps list size/order will never match with ensembl_symbol)
                    ensembl_symbol_list_from_mygene.append('symbol_not_found')

            if ensembl_symbol in ensembl_symbol_list_from_mygene:
                if ensembl_symbol_list_from_mygene.count(ensembl_symbol) == 1:
                    final_mapping_file.write(key + '\t')
                    ncbi_idx = ensembl_symbol_list_from_mygene.index(ensembl_symbol)
                    final_mapping_file.write('\t' + ncbi_list[ncbi_idx] + '\n')
    final_mapping_file.close()


# Call all the functions above in order:
multi_mapping_dict = find_multiple_mappings_from_entrezgene_file("gene_ensembl__xref_entrezgene__dm.txt")
ensembl_dict = create_ensembl_gene_id_dict("gene_ensembl__gene__main.txt", multi_mapping_dict)
ensembl_dict_appended = find_ncbi_ids_from_gene2ensembl(ensembl_dict, "gene2ensembl")
write_mapping_ids_to_file(ensembl_dict_appended)
