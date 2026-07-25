#!/usr/bin/env python3
"""Detailed WebMM vs RDKit type diff with neighbor context."""
import json
from collections import Counter

rd = json.load(open('scripts/rdkit_atom_types.json'))
us = json.load(open('scripts/webmm_atom_types.json'))

# MMFF numeric type -> name (common ones)
NAMES = {
    1:'C_3',2:'C_2',3:'C_1',4:'C_AR',5:'H',6:'O_3',7:'O_2',8:'N_2',9:'N_1',
    10:'N_AM',11:'N_AR',12:'N_PL3',13:'S? ',15:'S?',16:'S_O2',19:'Si',20:'C_CAT',
    21:'H_ONC',23:'H_N3',24:'H_COOH',25:'?',27:'H_NIM',28:'H_NAM',29:'H_OAR',
    30:'C?',32:'O_CO2',37:'C_AR',38:'N_AR',39:'NPYL',40:'N_PL3',41:'C_CO2',
    43:'N_SO2',45:'N_NO2',63:'C5A',64:'C5B',65:'N5A',66:'N5B',69:'N_POX',
    70:'OH2',71:'H?',89:'F_M',90:'CL_M',
}

def n(t):
    return f'{t}({NAMES.get(t,"?")})'

for name in sorted(set(rd) & set(us)):
    r_entry, u_entry = rd[name], us[name]
    if not isinstance(r_entry, list) or not isinstance(u_entry, list): continue
    for mi,(rm,um) in enumerate(zip(r_entry,u_entry)):
        if not isinstance(rm,dict) or not isinstance(um,dict): continue
        if 'error' in rm or 'error' in um: continue
        rt,ut = rm.get('types',[]), um.get('types',[])
        rsyms = rm.get('syms',[])
        bonds = rm.get('bonds',[])
        # build neighbor map
        nbrs = {i:[] for i in range(len(rsyms))}
        for a,b,o in bonds:
            nbrs[a].append((b,rsyms[b],o)); nbrs[b].append((a,rsyms[a],o))
        for i,(r,u) in enumerate(zip(rt,ut)):
            if int(r)!=int(u):
                nb = ' '.join(f'{s}{j}(b{o})' for j,s,o in nbrs[i])
                print(f'{name}#{mi} atom{i} {rsyms[i]}: rdkit={n(int(r))}  webmm={n(int(u))}  nbrs:[{nb}]')
