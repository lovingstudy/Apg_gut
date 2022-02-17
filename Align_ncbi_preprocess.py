# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 18:57:04 2021

@author: xingjin1
"""
import pandas as pd

def load_fasta(filename):
    buff = []
    with open(filename) as f:
        for line in f: buff.append(line.strip())
    seqs, key = dict(), ''
    for line in buff:
        if line.startswith('>'):
            key = line
            seqs[key] = []
        else:
            seqs[key].append(line)
    del buff
    return seqs


fasta_inp = 'Mah1_blast_ncbi_nr_idt80_small_align.fasta'
meta_file = 'Mah1_blast_ncbi_nr_idt80_small_sele.csv'
seqs = load_fasta(fasta_inp)
meta = pd.read_csv(meta_file, index_col='Accession')
new_seqs = dict()
for k, v in seqs.items():
    if 'Mah1' in k:
        k = '>Mah1'
    else:
        tmp = k.split(' ')
        k = '>' + meta.loc[tmp[1], 'Description']
        k = k.strip().replace(' ', '_').replace('[', '').replace(']', '').replace(':', '')
    new_seqs[k] = v

with open(fasta_inp.replace('.fasta', '_chkd.fasta'), 'w') as f:
    for k, v in new_seqs.items():
        f.write('\n'.join([k] + v) + '\n')

