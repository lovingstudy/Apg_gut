# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 16:27:28 2021

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


#query_seq_file = 'SEQ.fasta'
blast_seq_file = 'Mah1_blast_ncbi_nr_5k.fasta'
blast_hit_file = 'Mah1_blast_ncbi_nr_5k_chkd.csv'

hit_table = pd.read_csv(blast_hit_file)
hit_table = hit_table[hit_table['Per. ident'] >= 80]
hit_table = hit_table.sort_values(by='Per. ident', ascending=False)
#hit_table = hit_table[~hit_table['Description'].str.contains('unclassified')]
#hit_table = hit_table[~hit_table['Description'].str.contains('MULTISPECIES')]
hit_table = hit_table[~hit_table['Description'].duplicated()]
hit_table['Genus'] = [s.split(' ')[0] for s in hit_table['Scientific Name']]
#hit_table['Sub_Name'] = [s.replace(t, '').strip() for s, t in zip(hit_table['Scientific Name'], hit_table['TYPE'])]
hit_table['Protein_Name'] = [s.split('[')[0].strip() + ' ' + t + ' ' + a for s, t, a in zip(hit_table['Description'], hit_table['Genus'], hit_table['Accession'])]
#viz_table = pd.concat([hit_table[hit_table['Genus'] == ti].iloc[:5] for ti in hit_table['Genus'].value_counts().index])
viz_table = hit_table.copy()
viz_table['Protein_Name'] = [s.replace('maltodextrin glucosidase', 'Apg') for s in viz_table['Protein_Name']]
viz_table.to_csv('Apg_blast_ncbi_nr_idt80_full.csv')

blast_seqs = load_fasta(blast_seq_file)
sele_set = set(viz_table['Accession'].values)
sele_seqs = [(k, v) for k, v in blast_seqs.items() if any([a in k for a in sele_set])]

tmppool = set()
with open('Apg_blast_ncbi_nr_idt80_full.fasta', 'w') as f:
    #tmp = [l for l in open(query_seq_file)]
    #f.write(''.join(tmp))
    for k, v in sele_seqs:
        k1 = k.replace('maltodextrin glucosidase', 'Apg').replace('MULTISPECIES: ', '')
        #k1 = '>' + ' '.join(k1.split(' ')[1:])
        if k1 in tmppool: continue
        f.write('\n'.join([k1] + v) + '\n')
        tmppool.add(k1)





