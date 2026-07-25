#!/usr/bin/env python3
"""Detailed WebMM vs RDKit type diff for the validation set, with neighbor context."""
import json
from collections import Counter

rd = json.load(open('scripts/val_set/rdkit_ref.json'))
us = json.load(open('scripts/val_set/webmm_ref.json'))

NAMES = {1:'C_3',2:'C_VIN',3:'C_2',4:'C_1',5:'H',6:'O_3',7:'O_2',8:'N_3',9:'N_2',
         10:'N_AM',16:'S_2',20:'CR4R',21:'H_ONC',22:'CR3R',23:'H_N3',24:'H_COOH',
         25:'P_4',26:'P_3',28:'H_NAM',30:'CE4R',32:'O_CO2',33:'HOS',34:'N_4',
         36:'HNR+',37:'C_AR',38:'N_AR',39:'NPYL',40:'N_PL3',44:'S_AR',59:'OFUR',
         63:'C5A',64:'C5B',65:'N5A',66:'N5B',70:'OH2',71:'HS',72:'S2CM'}

# we don't have bonds in rdkit_ref; load SDFs for neighbor info
import os
def neighbors(name, sdf_dir='scripts/val_set'):
    p = f'{sdf_dir}/{name}.sdf'
    if not os.path.exists(p):
        return {}, []
    txt = open(p).read()
    lines = txt.split('\n')
    try:
        counts = lines[3]
        na = int(counts[0:3]); nb = int(counts[3:6])
    except Exception:
        return {}, []
    syms = [lines[4+i][31:34].strip() for i in range(na)]
    nbrs = {i: [] for i in range(na)}
    for k in range(nb):
        b = lines[4+na+k]
        a = int(b[0:3])-1; c = int(b[3:6])-1
        nbrs.setdefault(a, []).append((c, syms[c] if c < len(syms) else '?'))
        nbrs.setdefault(c, []).append((a, syms[a] if a < len(syms) else '?'))
    return nbrs, syms

def n(t): return f'{t}({NAMES.get(t,"?")})'

for name in sorted(set(rd) & set(us)):
    r, u = rd[name], us[name]
    if 'error' in r or 'error' in u:
        continue
    rt, ut = r['types'], u['types']
    nbrs, syms = neighbors(name)
    for i, (a, b) in enumerate(zip(rt, ut)):
        if int(a) != int(b):
            nb = ','.join(f'{s}{j}' for j, s in nbrs.get(i, []))
            print(f'{name} atom{i} {syms[i] if i<len(syms) else "?"}: rdkit={n(int(a))} webmm={n(int(b))} nbrs=[{nb}]')
