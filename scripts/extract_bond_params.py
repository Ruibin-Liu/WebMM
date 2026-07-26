"""Extract MMFF bond params (r0, kb) from RDKit by finite-differencing val_set SDFs."""
import math
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def bond_params(sdf, i, j):
    def make():
        m = Chem.MolFromMolBlock(sdf, sanitize=False)
        m.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        Chem.SetAromaticity(m)
        return m
    m0 = make(); conf = m0.GetConformer()
    pi = conf.GetAtomPosition(i); pj = conf.GetAtomPosition(j)
    dx, dy, dz = pj.x-pi.x, pj.y-pi.y, pj.z-pi.z
    n = math.sqrt(dx*dx+dy*dy+dz*dz); dx/=n; dy/=n; dz/=n
    def E_at_r(r):
        m = make(); mp = AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94s')
        m.GetConformer().SetAtomPosition(j, (pi.x+dx*r, pi.y+dy*r, pi.z+dz*r))
        return AllChem.MMFFGetMoleculeForceField(m, mp, confId=0).CalcEnergy()
    rs = [1.2 + 0.02*k for k in range(40)]
    es = [E_at_r(r) for r in rs]
    emin = min(es); r0 = rs[es.index(emin)]
    rs2 = [r0-0.01+0.001*k for k in range(21)]
    es2 = [E_at_r(r) for r in rs2]
    r0 = rs2[es2.index(min(es2))]
    h = 0.005
    d2 = (E_at_r(r0+h) - 2*E_at_r(r0) + E_at_r(r0-h)) / (h*h)
    kb = d2 / 143.9325
    return r0, kb

def load(name):
    return open(f'scripts/val_set/{name}.sdf').read().split('$$$$')[0].rstrip() + '\nM  END\n'

thio = load('thioacetone')
r0, kb = bond_params(thio, 1, 2)
print(f'thioacetone C=S: r0={r0:.3f} kb={kb:.3f}  [webmm r0=1.610 kb=5.500]')

mpa = load('methylphosphonic_acid')
for label, (i,j) in [('C-P',(0,1)),('P=O',(1,2)),('P-O',(1,3))]:
    r0, kb = bond_params(mpa, i, j)
    print(f'methylphosphonic {label}: r0={r0:.3f} kb={kb:.3f}')

cs2 = load('carbon_disulfide')
r0, kb = bond_params(cs2, 0, 1)
print(f'carbon_disulfide C=S: r0={r0:.3f} kb={kb:.3f}')
