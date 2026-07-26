"""Compare RDKit verbose torsion output to WebMM per-torsion for adamantane."""
import re, subprocess
# Parse RDKit verbose (already captured in /tmp/adammff.log)
rdkit = {}
in_tor = False
for line in open('/tmp/adammff.log'):
    if 'T O R S I O N A L' in line:
        in_tor = True; continue
    if in_tor and ('V A N' in line or 'TOTAL' in line):
        if 'TOTAL' in line: break
        in_tor = False; continue
    if not in_tor:
        continue
    nums = re.findall(r'#(\d+)', line)
    if len(nums) < 4:
        continue
    i,j,k,l = int(nums[0])-1, int(nums[1])-1, int(nums[2])-1, int(nums[3])-1
    # last 6 numeric fields: angle, energy, V1, V2, V3 (5 actually)
    fields = line.split()
    try:
        ang = float(fields[-5]); e = float(fields[-4])
    except (ValueError, IndexError):
        continue
    key = (min(i,l), min(j,k), max(j,k), max(i,l))
    rdkit[key] = (ang, e)

# Get WebMM per-torsion via a Rust example
rust_code = '''
use webmm::mmff::MMFFForceField;
use webmm::mmff::torsion::{get_torsion_params, torsion_energy};
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;
fn main(){
 let mol=parse_sdf(&std::fs::read_to_string("scripts/val_set/adamantane.sdf").unwrap()).unwrap();
 let c:Vec<[f64;3]>=mol.atoms.iter().map(|a|a.position).collect();
 let ff=MMFFForceField::new(&mol,MMFFVariant::MMFF94s);
 for (ti,t) in ff.torsions.iter().enumerate(){
  let p=get_torsion_params(ff.atom_types[t.atom1],ff.atom_types[t.atom2],ff.atom_types[t.atom3],ff.atom_types[t.atom4],ff.torsion_types[ti],ff.variant);
  if let Some(pp)=p{
   let e=torsion_energy(&c,t.atom1,t.atom2,t.atom3,t.atom4,&pp);
   let k=(t.atom1.min(t.atom4),t.atom2.min(t.atom3),t.atom2.max(t.atom3),t.atom1.max(t.atom4));
   println!("{} {} {} {} {:.4}", k.0,k.1,k.2,k.3, e);
  }
 }
}'''
open('/tmp/wt.rs','w').write(rust_code)
import shutil, os
shutil.copy('/tmp/wt.rs', 'examples/wt.rs')
out = subprocess.run(['cargo','run','--release','--example','wt'], capture_output=True, text=True).stdout
os.remove('examples/wt.rs')
webmm = {}
for line in out.split('\n'):
    p = line.split()
    if len(p)==5:
        try:
            key = (int(p[0]),int(p[1]),int(p[2]),int(p[3])); e=float(p[4])
            webmm[key]=e
        except: pass

# compare
print(f"RDKit torsions: {len(rdkit)}, WebMM torsions: {len(webmm)}")
rd_sum = sum(v[1] for v in rdkit.values())
wm_sum = sum(webmm.values())
print(f"RDKit total torsion: {rd_sum:.2f}, WebMM total: {wm_sum:.2f}")
# biggest mismatches
mism=[]
for key in rdkit:
    re_ang, re_e = rdkit[key]
    we = webmm.get(key)
    if we is None:
        continue
    mism.append((abs(we-re_e), key, re_e, we, re_ang))
mism.sort(reverse=True)
print("biggest |ΔE| mismatches (key, rdkit_E, webmm_E, rdkit_angle):")
for d,key,re,we,ang in mism[:12]:
    print(f"  {key}: rdkit={re:+.3f} webmm={we:+.3f} Δ={we-re:+.3f} (rdkit angle={ang:.1f})")
