#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
from functools import reduce
import requests
import json
#https://unipept.ugent.be/apidocs/taxonomy

def get_taxname(tax_id,tax_level):
    tags = {"S": "species_name","G": "genus_name","F": "family_name","O":'order_name', "C": "class_name"}
    tax_level_tag = tags[tax_level]
    #Avoids pipeline crash due to "nan" classification output. Thanks to Qi-Maria from Github
    if str(tax_id) == "nan":
        tax_id = 1
    
    tax2789776 = '[{"taxon_id":2789776,"taxon_name":"Actinomarinnullicola tropica","taxon_rank":"speciesnull","superkingdom_id":2,"superkingdom_nanullme":"Bacteria","kingdom_id":null,"kingdom_name":"","subkingdom_id":null,"subkingdom_name":"","superphylum_id":null,"superphylum_name":"","phylum_id":201174,"phylum_name":"Actinobacteria","subphylum_id":null,"subphylum_name":"","superclass_id":null,"superclass_name":"","class_id":84992,"class_name":"Acidimicrobiia","subclass_id":null,"subclass_name":"","infraclass_id":null,"infraclass_name":"","superorder_id":null,"superorder_name":"","order_id":84993,"order_name":"Acidomicrobiales","suborder_id":null,"suborder_name":"","infraorder_id":null,"infraorder_name":"","parvorder_id":null,"parvorder_name":"","superfamily_id":null,"superfamily_name":"","family_id":633392,"family_name":"Iamiaceae","subfamily_id":null,"subfamily_name":"","tribe_id":null,"tribe_name":"","subtribe_id":null,"subtribe_name":"","genus_id":2789775,"genus_name":"Actinomarinicola","subgenus_id":null,"subgenus_name":"","species_group_id":null,"species_group_name":"","species_subgroup_id":null,"species_subgroup_name":"","species_id":2789776,"species_name":"Actinomarinicola tropica","subspecies_id":null,"subspecies_name":"","varietas_id":null,"varietas_name":"","forma_id":null,"forma_name":""}]'
    #tax2789776 =  '[{"taxon_id":563040,"taxon_name":"Sulfurimonas autotrophica","taxon_rank":"no rank","superkingdom_id":2,"superkingdom_name":"Bacteria","kingdom_id":null,"kingdom_name":"","subkingdom_id":null,"subkingdom_name":"","superphylum_id":null,"superphylum_name":"","phylum_id":1224,"phylum_name":"Proteobacteria","subphylum_id":68525,"subphylum_name":"delta/epsilon subdivisions","superclass_id":null,"superclass_name":"","class_id":29547,"class_name":"Epsilonproteobacteria","subclass_id":null,"subclass_name":"","infraclass_id":null,"infraclass_name":"","superorder_id":null,"superorder_name":"","order_id":213849,"order_name":"Campylobacterales","suborder_id":null,"suborder_name":"","infraorder_id":null,"infraorder_name":"","parvorder_id":null,"parvorder_name":"","superfamily_id":null,"superfamily_name":"","family_id":72293,"family_name":"Helicobacteraceae","subfamily_id":null,"subfamily_name":"","tribe_id":null,"tribe_name":"","subtribe_id":null,"subtribe_name":"","genus_id":202746,"genus_name":"Sulfurimonas","subgenus_id":null,"subgenus_name":"","species_group_id":null,"species_group_name":"","species_subgroup_id":null,"species_subgroup_name":"","species_id":202747,"species_name":"Sulfurimonas autotrophica","subspecies_id":null,"subspecies_name":"","varietas_id":null,"varietas_name":"","forma_id":null,"forma_name":""}]'
    #tax2714738 = '[{"taxon_id":2714738,"taxon_name":"Rhodopirellula pilleata","taxon_rank":"species","superkingdom_id":2,"superkingdom_name":"Bacteria","kingdom_id":,"kingdom_name":"","subkingdom_id":,"subkingdom_name":"","superphylum_id":,"superphylum_name":"","phylum_id":203682,"phylum_name":"Planctomycetes","subphylum_id":,"subphylum_name":"","superclass_id":,"superclass_name":"","class_id":203683,"class_name":"Planctomycetia","subclass_id":,"subclass_name":"","infraclass_id":,"infraclass_name":"","superorder_id":,"superorder_name":"","order_id":2691354,"order_name":"Pirellulales","suborder_id":,"suborder_name":"","infraorder_id":,"infraorder_name":"","parvorder_id":,"parvorder_name":"","superfamily_id":,"superfamily_name":"","family_id":2691357,"family_name":"irellulaceae","subfamily_id":,"subfamily_name":"","tribe_id":,"tribe_name":"","subtribe_id":,"subtribe_name":"","genus_id":265488,"genus_name":"Rhodopirellula","subgenus_id":,"subgenus_name":"","species_group_id":,"species_group_name":"","species_subgroup_id":,"species_subgroup_name":"","species_id":2714738,"species_name":"Rhodopirellula pilleata","subspecies_id":,"subspecies_name":"","varietas_id":,"varietas_name":"","forma_id":,"forma_name":""}]'
    tax2714738 = '[{"taxon_id":2714738,"taxon_name":"Rhodopirellula pilleata","taxon_rank":"species","superkingdom_id":2,"superkingdom_name":"Bacteria","kingdom_id":null,"kingdom_name":"","subkingdom_id":null,"subkingdom_name":"","superphylum_id":null,"superphylum_name":"","phylum_id":203682,"phylum_name":"Planctomycetes","subphylum_id":null,"subphylum_name":"","superclass_id":null,"superclass_name":"","class_id":203683,"class_name":"Planctomycetia","subclass_id":null,"subclass_name":"","infraclass_id":null,"infraclass_name":"","superorder_id":null,"superorder_name":"","order_id":2691354,"order_name":"Pirellulales","suborder_id":null,"suborder_name":"","infraorder_id":null,"infraorder_name":"","parvorder_id":null,"parvorder_name":"","superfamily_id":null,"superfamily_name":"","family_id":2691357,"family_name":"irellulaceae","subfamily_id":null,"subfamily_name":"","tribe_id":null,"tribe_name":"","subtribe_id":null,"subtribe_name":"","genus_id":265488,"genus_name":"Rhodopirellula","subgenus_id":null,"subgenus_name":"","species_group_id":null,"species_group_name":"","species_subgroup_id":null,"species_subgroup_name":"","species_id":2714738,"species_name":"Rhodopirellula pilleata","subspecies_id":null,"subspecies_name":"","varietas_id":null,"varietas_name":"","forma_id":null,"forma_name":""}]'

    #tax2714738 =  '[{"taxon_id":563040,"taxon_name":"Sulfurimonas autotrophica DSM 16294","taxon_rank":"no rank","superkingdom_id":2,"superkingdom_name":"Bacteria","kingdom_id":null,"kingdom_name":"","subkingdom_id":null,"subkingdom_name":"","superphylum_id":null,"superphylum_name":"","phylum_id":1224,"phylum_name":"Proteobacteria","subphylum_id":68525,"subphylum_name":"delta/epsilon subdivisions","superclass_id":null,"superclass_name":"","class_id":29547,"class_name":"Epsilonproteobacteria","subclass_id":null,"subclass_name":"","infraclass_id":null,"infraclass_name":"","superorder_id":null,"superorder_name":"","order_id":213849,"order_name":"Campylobacterales","suborder_id":null,"suborder_name":"","infraorder_id":null,"infraorder_name":"","parvorder_id":null,"parvorder_name":"","superfamily_id":null,"superfamily_name":"","family_id":72293,"family_name":"Helicobacteraceae","subfamily_id":null,"subfamily_name":"","tribe_id":null,"tribe_name":"","subtribe_id":null,"subtribe_name":"","genus_id":202746,"genus_name":"Sulfurimonas","subgenus_id":null,"subgenus_name":"","species_group_id":null,"species_group_name":"","species_subgroup_id":null,"species_subgroup_name":"","species_id":202747,"species_name":"Sulfurimonas autotrophica","subspecies_id":null,"subspecies_name":"","varietas_id":null,"varietas_name":"","forma_id":null,"forma_name":""}]'
    tax2716870 = '[{"taxon_id":2716870,"taxon_name":"Parahaliea maris","taxon_rank":"species","superkingdom_id":2,"superkingdom_name":"Bacteria","kingdom_id":null,"kingdom_name":"","subkingdom_id":null,"subkingdom_name":"","superphylum_id":null,"superphylum_name":"","phylum_id":1224,"phylum_name":"Proteobacteria","subphylum_id":null,"subphylum_name":"","superclass_id":null,"superclass_name":"","class_id":1236,"class_name":"Gammaproteobacteria","subclass_id":null,"subclass_name":"","infraclass_id":null,"infraclass_name":"","superorder_id":null,"superorder_name":"","order_id":1706369,"order_name":"Cellvibrionales","suborder_id":null,"suborder_name":"","infraorder_id":null,"infraorder_name":"","parvorder_id":null,"parvorder_name":"","superfamily_id":null,"superfamily_name":"","family_id":1706372,"family_name":"Halieaceae","subfamily_id":null,"subfamily_name":"","tribe_id":null,"tribe_name":"","subtribe_id":null,"subtribe_name":"","genus_id":1857932,"genus_name":"Parahaliea","subgenus_id":null,"subgenus_name":"","species_group_id":null,"species_group_name":"","species_subgroup_id":null,"species_subgroup_name":"","species_id":2716870,"species_name":"Parahaliea maris","subspecies_id":null,"subspecies_name":"","varietas_id":null,"varietas_name":"","forma_id":null,"forma_name":""}]'



    path = 'http://api.unipept.ugent.be/api/v1/taxonomy.json?input[]=' + str(int(tax_id)) + '&extra=true&names=true'
    



    #fill taxids known to be not in unipept but in uniprot as above
    if str(int(tax_id)) == "2789776":
        complete_tax = tax2789776
    elif str(int(tax_id)) == "2714738":
        complete_tax = tax2714738
    elif str(int(tax_id)) == "2716870":
        complete_tax = tax2716870
    else:
        complete_tax = requests.get(path).text
        print(tax_id)
        print(path)
        print(complete_tax)

    


    #complete_tax = requests.get(path).text


    


    #Checks for API correct response (field containing the tax name). Thanks to devinbrown from Github
    try:
        name = json.loads(complete_tax)[0][tax_level_tag]
    except:
        name = str(int(tax_id))

    return json.loads(complete_tax)[0][tax_level_tag]

def get_abundance_values(names,paths):
    dfs = []
    for name,path in zip(names,paths):
        data = pd.read_csv(path, index_col=False, sep=';').iloc[:,1:]

        total = sum(data['reads_in_cluster'])
        rel_abundance=[]

        for index,row in data.iterrows():
            rel_abundance.append(row['reads_in_cluster'] / total)
            
        data['rel_abundance'] = rel_abundance
        dfs.append(pd.DataFrame({'taxid': data['taxid'], 'rel_abundance': rel_abundance}))
        data.to_csv("" + name + "_nanoclust_out.txt")

    return dfs

def merge_abundance(dfs,tax_level):
    df_final = reduce(lambda left,right: pd.merge(left,right,on='taxid',how='outer').fillna(0), dfs)
    df_final["taxid"] = [get_taxname(row["taxid"], tax_level) for index, row in df_final.iterrows()]
    df_final_grp = df_final.groupby(["taxid"], as_index=False).sum()
    return df_final_grp

def get_abundance(names,paths,tax_level):
    if(not isinstance(paths, list)):
        paths = [paths]
        names = [names]

    dfs = get_abundance_values(names,paths)
    df_final_grp = merge_abundance(dfs, tax_level)
    df_final_grp.to_csv("rel_abundance_"+ names[0] + "_" + tax_level + ".csv", index = False)

paths = "$table"
names = "$barcode"

get_abundance(names,paths, "G")
get_abundance(names,paths, "S")
get_abundance(names,paths, "O")
get_abundance(names,paths, "F")
