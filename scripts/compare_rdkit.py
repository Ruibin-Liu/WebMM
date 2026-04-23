#!/usr/bin/env python3
"""Compare WebMM ETKDG and MMFF94s against RDKit reference."""

import subprocess
import sys
import math

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

TEST_MOLECULES = {
    'water': '''Water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END''',
    'acetic': '''Acetic acid
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -0.9335   -0.0601   -0.2304 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4936    0.2789    0.0469 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0325    1.3566   -0.1361 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1814   -0.7645    0.5462 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4427    0.8203   -0.6327 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4305   -0.3576    0.6963 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9867   -0.8625   -0.9702 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0859   -0.4112    0.6800 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
M  END''',
    'benzene': '''Benzene
     RDKit          3D

 12 12  0  0  0  0  0  0  0  0999 V2000
    0.8035   -1.1401   -0.0082 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3889    0.1258   -0.0281 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5853    1.2659   -0.0199 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8035    1.1401    0.0083 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3889   -0.1258    0.0281 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5853   -1.2659    0.0199 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4295   -2.0284   -0.0147 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4709    0.2238   -0.0501 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0414    2.2522   -0.0354 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4295    2.0284    0.0147 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4709   -0.2238    0.0500 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0414   -2.2522    0.0354 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  2  8  1  0
  3  9  1  0
  4 10  1  0
  5 11  1  0
  6 12  1  0
M  END''',
    'ethanol': '''Ethanol
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.8883    0.1670   -0.0273 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4658   -0.5116   -0.0368 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4311    0.3229    0.5867 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8487    1.1175   -0.5695 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6471   -0.4704   -0.4896 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1964    0.3978    0.9977 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7920   -0.7224   -1.0597 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4246   -1.4559    0.5138 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4671    1.1550    0.0848 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END''',
    'aniline': '''Aniline
     RDKit          3D

 13 13  0  0  0  0  0  0  0  0999 V2000
   -1.2000    0.6930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -0.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3000   -1.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3000    0.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1500    1.3000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2500    2.2800    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.8000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4000   -1.1000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3500   -2.2700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7500   -1.5500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3000    1.0500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  7  8  1  0
  7  9  1  0
  2 10  1  0
  3 11  1  0
  4 12  1  0
  5 13  1  0
M  END''',
}

def run_rust_etkdg(sdf):
    proc = subprocess.run(
        ['./target/debug/compare_etkdg'],
        input=sdf,
        capture_output=True,
        text=True,
        cwd='/Users/rliu/Projects/WebMM'
    )
    if proc.returncode != 0:
        print(f"  Rust ETKDG failed: {proc.stderr}")
        return None
    lines = proc.stdout.strip().split('\n')
    n_atoms = int(lines[0].split()[0])
    coords = []
    for line in lines[1:]:
        parts = line.split()
        coords.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return coords

def run_rust_mmff(sdf):
    proc = subprocess.run(
        ['./target/debug/compare_mmff'],
        input=sdf,
        capture_output=True,
        text=True,
        cwd='/Users/rliu/Projects/WebMM'
    )
    if proc.returncode != 0:
        print(f"  Rust MMFF failed: {proc.stderr}")
        return None
    lines = proc.stdout.strip().split('\n')
    energy = float(lines[0].split()[1])
    n_atoms = int(lines[1].split()[1])
    coords = []
    grads = []
    for line in lines[2:]:
        parts = line.split()
        coords.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
        # format: "... | grad: x y z | norm: ..."
        grad_idx = parts.index('grad:')
        grads.append((float(parts[grad_idx+1]), float(parts[grad_idx+2]), float(parts[grad_idx+3])))
    return energy, coords, grads

def rdkit_etkdg(sdf, seed=42):
    mol = Chem.MolFromMolBlock(sdf, removeHs=False)
    if mol is None:
        return None
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.maxIterations = 100
    params.numThreads = 1
    params.clearConfs = True
    AllChem.EmbedMolecule(mol, params)
    conf = mol.GetConformer()
    coords = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append((mol.GetAtomWithIdx(i).GetSymbol(), pos.x, pos.y, pos.z))
    return coords

def rdkit_mmff(sdf, seed=42):
    mol = Chem.MolFromMolBlock(sdf, removeHs=False)
    if mol is None:
        return None, None, None
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.maxIterations = 100
    AllChem.EmbedMolecule(mol, params)
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    if ff is None:
        return None, None, None
    ff.Initialize()
    energy_before = ff.CalcEnergy()
    
    # Get gradient
    grad_rdkit = []
    for i in range(mol.GetNumAtoms()):
        pos = ff.Positions()
        # We can't easily get per-atom gradients from RDKit Python API
        # So we'll skip direct gradient comparison
        grad_rdkit.append((0.0, 0.0, 0.0))
    
    ff.Minimize(maxIts=500)
    energy_after = ff.CalcEnergy()
    conf = mol.GetConformer()
    coords = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append((mol.GetAtomWithIdx(i).GetSymbol(), pos.x, pos.y, pos.z))
    return coords, (energy_before, energy_after), grad_rdkit

def compute_distance_matrix(coords):
    n = len(coords)
    mat = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            d = math.dist(coords[i][1:], coords[j][1:])
            mat[i][j] = d
            mat[j][i] = d
    return mat

def rmsd(coords1, coords2):
    assert len(coords1) == len(coords2)
    n = len(coords1)
    total = 0.0
    for i in range(n):
        for j in range(1, 4):
            total += (coords1[i][j] - coords2[i][j])**2
    return math.sqrt(total / n)

def compare_gradients(grad1, grad2, label, threshold=1.0):
    max_diff = 0.0
    max_idx = -1
    for i, (g1, g2) in enumerate(zip(grad1, grad2)):
        for j in range(3):
            diff = abs(g1[j] - g2[j])
            if diff > max_diff:
                max_diff = diff
                max_idx = i
    print(f"  {label} max gradient diff: {max_diff:.4f} (atom {max_idx})")
    if max_diff > threshold:
        for i, (g1, g2) in enumerate(zip(grad1, grad2)):
            for j in range(3):
                diff = abs(g1[j] - g2[j])
                if diff > threshold:
                    print(f"    atom {i} dim {j}: ours={g1[j]:.6f} rdkit={g2[j]:.6f} diff={diff:.6f}")
    return max_diff

def compare_molecule(name, sdf):
    print(f"\n{'='*60}")
    print(f"Molecule: {name}")
    print(f"{'='*60}")
    
    # ETKDG comparison
    rust_coords = run_rust_etkdg(sdf)
    rdkit_coords = rdkit_etkdg(sdf)
    
    if rust_coords is None or rdkit_coords is None:
        print("  FAILED to generate coordinates")
        return
    
    rust_dm = compute_distance_matrix(rust_coords)
    rdkit_dm = compute_distance_matrix(rdkit_coords)
    
    max_bond_diff = 0.0
    for i in range(len(rust_coords)):
        for j in range(i+1, len(rust_coords)):
            d = rust_dm[i][j]
            if d < 2.2:
                diff = abs(rust_dm[i][j] - rdkit_dm[i][j])
                if diff > max_bond_diff:
                    max_bond_diff = diff
    
    print(f"  ETKDG max bond length diff: {max_bond_diff:.4f}")
    
    # MMFF comparison on our ETKDG coords
    rust_energy, rust_mmff_coords, rust_grad = run_rust_mmff(sdf)
    rdkit_mmff_coords, rdkit_energies, rdkit_grad = rdkit_mmff(sdf)
    
    if rust_energy is not None and rdkit_energies is not None:
        print(f"  MMFF94s energy (Rust ETKDG start): {rust_energy:.4f}")
        print(f"  MMFF94s energy (RDKit ETKDG start, before opt): {rdkit_energies[0]:.4f}")
        print(f"  MMFF94s energy (RDKit ETKDG start, after opt):  {rdkit_energies[1]:.4f}")
        
        # Compare gradient norms at our ETKDG starting positions
        rust_max_g = max(math.sqrt(g[0]**2 + g[1]**2 + g[2]**2) for g in rust_grad)
        print(f"  Max gradient norm (Rust):  {rust_max_g:.4f}")
    
    # Compare optimized structures: optimize our ETKDG with our MMFF, then compare to RDKit optimized
    # We can also run RDKit MMFF starting from our coords for a fair comparison
    mol = Chem.MolFromMolBlock(sdf, removeHs=False)
    if mol is not None and rust_coords is not None:
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, (rust_coords[i][1], rust_coords[i][2], rust_coords[i][3]))
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
        ff = AllChem.MMFFGetMoleculeForceField(mol, props)
        if ff is not None:
            ff.Initialize()
            e_ours_start = ff.CalcEnergy()
            ff.Minimize(maxIts=500)
            e_ours_end = ff.CalcEnergy()
            print(f"  RDKit MMFF on Rust ETKDG: {e_ours_start:.4f} -> {e_ours_end:.4f}")

def main():
    print("Comparing WebMM against RDKit")
    for name, sdf in TEST_MOLECULES.items():
        compare_molecule(name, sdf)

if __name__ == '__main__':
    main()
