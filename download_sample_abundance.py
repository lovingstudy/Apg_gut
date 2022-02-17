# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 14:22:04 2022

@author: xingjin1
"""
import requests
import json
import pandas as pd

meta_table = pd.read_csv('all_runs_associated_with_D006262.csv', index_col=1)
outdir = './Healthy/'

meta_table = meta_table[(meta_table['QC status'] == 1) & (meta_table['recent antibiotics use'] != 'Y')]

url = 'https://gmrepo.humangut.info/api/getFullTaxonomicProfileByRunID'
for num, rid in enumerate(meta_table.index):
    query = {'run_id': rid}
    try:
        data = requests.post(url, data=json.dumps(query)).json()
        species = pd.DataFrame(data.get("species"))
    except Exception as e:
        print(rid, 'failed.', num)
        print(e)
        continue
    if len(species) == 0: continue
    species.to_csv(outdir + '%s_species_abundance.csv'%rid)
    
    if num % 500 == 0: print(num)
    
    