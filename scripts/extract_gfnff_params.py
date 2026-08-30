#!/usr/bin/env python3
"""Extract GFN-FF parameters from the xtb source tree into data/gfnff_params.json.

Source: github.com/thfroitzheim/xtb branch gxtb (LGPL-3.0), files:
  src/gfnff/gfnff_param.f90    (element tables + generator defaults)
  src/gfnff/gfnff_ini.f90      (zeta() constants)
  src/param/covalentradd3.f90  (covalent radii, D3)
  src/param/sqrtzr4r2.f90      (sqrt(<r4>/<r2>) table)
  include/param_ref.fh         (D4 reference-system data)

Usage: python3 scripts/extract_gfnff_params.py [path-to-xtb-src] (default /tmp/gxtb-src/xtb-gxtb)
"""
import json, re, sys, os

SRC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/gxtb-src/xtb-gxtb"

def read(p): return open(os.path.join(SRC, p)).read()

def parse_array(text, name, count=None):
    """parse Fortran array init:  :: name(103) = [& ... &]  or  (/ ... /)"""
    m = re.search(r"::\s*%s\b[^=]*=" % name, text)
    if not m: raise KeyError(name)
    rest = text[m.end():]
    # find terminator: ']' or '/)' whichever first
    end = min([p for p in [rest.find("]"), rest.find("/)")] if p >= 0])
    body = rest[:end].replace("D", "E").replace("d", "e")
    vals = [float(x) for x in re.findall(r"-?\d+\.\d*(?:[eE][-+]?\d+)?", body)]
    return vals

param = read("src/gfnff/gfnff_param.f90")
ini   = read("src/gfnff/gfnff_ini.f90")

out = {}
for key in ["chi_angewChem2020","gam_angewChem2020","cnf_angewChem2020","alp_angewChem2020",
            "bond_angewChem2020","repa_angewChem2020","repan_angewChem2020",
            "angl_angewChem2020","angl2_angewChem2020","tors_angewChem2020","tors2_angewChem2020",
            "en","rad","repz"]:
    v = parse_array(param, key)
    assert len(v) == 103, (key, len(v))
    out[key.split("_")[0]] = v

for key, arr in [("metal",None),("group",None),("normcn",None)]:
    m = re.search(r"::\s*%s\b[^=]*=" % key, param)
    rest = param[m.end():]
    end = min(p for p in [rest.find("]"), rest.find("/)")] if p >= 0)
    v = [int(float(x)) for x in re.findall(r"-?\d+\.?\d*", rest[:end])]
    assert len(v) == 103, (key, len(v))
    out[key] = v

def grab_until_close(text, decl_re):
    m = re.search(decl_re, text)
    if not m: raise KeyError(decl_re)
    rest = text[m.end():]
    end = min(p for p in [rest.find("]"), rest.find("/)")] if p >= 0)
    return rest[:end].replace("D", "E").replace("d", "e")

body = grab_until_close(ini, r"zeff\(103\)[^=]*=")
out["zeta_zeff"] = [int(float(x)) for x in re.findall(r"-?\d+\.?\d*", body)]
body = grab_until_close(ini, r"c\(1:103\)[^=]*=")
out["zeta_c"] = [float(x) for x in re.findall(r"-?\d+\.\d*", body)]
assert len(out["zeta_zeff"]) == 103 and len(out["zeta_c"]) == 103

# ---- covalent radii (D3) and sqrt(Z r4r2) ----
cov = read("src/param/covalentradd3.f90")
out["rcov"] = parse_array(cov, "covalentRadD3")[:103]

r4r2 = read("src/param/sqrtzr4r2.f90")
# sqrtZr4r2(z) = sqrt(0.5 * r4Overr2(z) * sqrt(z))  (sqrtzr4r2.f90 line 68)
out["sqrt_zr4r2"] = [(0.5 * x * (z + 1) ** 0.5) ** 0.5
                     for z, x in enumerate(parse_array(r4r2, "r4Overr2")[:103])]

# ---- generator defaults (verified against newGFNFFGenerator, gfnff_param.f90) ----
out["generator"] = dict(
    cnmax=4.4, linthr=160.0, fcthr=1e-3, tdist_thr=12.0, rthr=1.25, rthr2=1.00,
    rqshrink=0.23, hqabthr=0.01, qabthr=0.10, srb1=0.3731, srb2=0.3171, srb3=0.2538,
    qrepscal=0.3480, nrepscal=-0.1270, hhfac=0.6290, hh13rep=1.4580, hh14rep=0.7080,
    bstren=[1.00,1.24,1.98,1.22,1.00,0.78,3.40,3.40,3.40],
    qfacBEN=-0.54, qfacTOR=12.0, fr3=0.3, fr4=1.0, fr5=1.5, fr6=5.7,
    torsf=[1.00,1.18,1.05,0.0,0.50,-0.90,0.70,-2.00],
    fbs1=0.50, batmscal=0.30, mchishift=-0.09, rabshift=-0.110, rabshifth=-0.050,
    hyper_shift=0.03, hshift3=-0.11, hshift4=-0.11, hshift5=-0.06,
    metal1_shift=0.2, metal2_shift=0.15, metal3_shift=0.05, eta_shift=0.040,
    qfacbm=[1.0,-0.2,-0.2,0.70,0.50], qfacbm0=0.047, rfgoed1=1.175,
    htriple=1.45, hueckelp2=1.00, hueckelp3=-0.24,
    hdiag={5:-0.5,6:0.0,7:0.14,8:-0.38,9:-0.29,16:-0.30,17:-0.30},
    hoffdiag={5:0.5,6:1.00,7:0.66,8:1.10,9:0.23,16:0.60,17:1.00},
    hiter=0.700, hueckelp=0.340, bzref=0.370, bzref2=0.315, pilpf=0.530, maxhiter=5,
    d3a1=0.58, d3a2=4.80, split0=0.670, fringbo=0.020, aheavy3=89.0, aheavy4=100.0,
)
g = out["generator"]
g["split1"] = 1.0 - g["split0"]
# bsmat (hyb_i x hyb_j) from bstren
s1,s2,s3 = g["bstren"][0], g["bstren"][1], g["bstren"][2]
bsmat = [[-999.0]*4 for _ in range(4)]
bsmat[0][0]=s1; bsmat[3][0]=s1; bsmat[3][3]=s1; bsmat[2][2]=s2; bsmat[1][1]=s3
bsmat[1][0]=g["split0"]*s1+g["split1"]*s3; bsmat[3][1]=bsmat[1][0]
bsmat[2][1]=g["split0"]*s2+g["split1"]*s3; bsmat[2][0]=g["split0"]*s1+g["split1"]*s2; bsmat[3][2]=bsmat[2][0]
g["bsmat"] = bsmat

# ---- gfnffrab tables (src/gfnff/gfnff_rab.f): en, r0, cnfak (103 each) ----
rab = read("src/gfnff/gfnff_rab.f")
def parse_data_block(text, name):
    m = re.search(r"data\s+%s\s*/" % name, text)
    if not m: raise KeyError(name)
    rest = text[m.end():]
    end = rest.find("/")
    body = rest[:end].replace("D","E").replace("d","e")
    vals = [float(x) for x in re.findall(r"-?\d+\.\d*", body)]
    return vals
out["rab_en"]    = parse_data_block(rab, "en")[:103]
out["rab_r0"]    = parse_data_block(rab, "r0")[:103]
out["rab_cnfak"] = parse_data_block(rab, "cnfak")[:103]

# ---- misc scalars from gfnff_set_param ----
out["misc"] = dict(
    cnmax=4.4, atcuta=0.595, atcutt=0.505, atcuta_nci=0.395, atcutt_nci=0.305,
    repscalb=1.7583, repscaln=0.4270, hbacut=49.0, hbscut=22.0, xbacut=70.0, xbscut=5.0,
    hbsf=1.0, hbst=15.0, xbsf=0.03, xbst=15.0, hbalp=6.0, hblongcut=85.0,
    hblongcut_xb=70.0, hbabmix=0.80, hbnbcut=11.20, tors_hb=0.94, bend_hb=0.20,
    vbond_scale=0.9, xhaci_globabh=0.268, xhaci_coh=0.350, xhaci_glob=1.50,
)
# HB basicity / acidity per element (xhbas, xhaci, xbaci)
xhbas = [0.0]*104; xhaci = [0.0]*104; xbaci = [0.0]*104
for z,v in {6:0.80,7:1.68,8:0.67,9:0.52,14:4.0,15:3.5,16:2.0,17:1.5,35:1.5,53:1.9,33:3.5,34:2.0,51:3.5,52:2.0}.items(): xhbas[z]=v
glob = out["misc"]["xhaci_glob"]
for z,v in {6:0.75,7:glob+0.1,8:glob,9:glob,15:glob,16:glob,17:glob+1.0,35:glob+1.0,53:glob+1.0}.items(): xhaci[z]=v
for z,v in {15:1.0,16:1.0,17:0.5,33:1.2,34:1.2,35:0.9,51:1.2,52:1.2,53:1.2}.items(): xbaci[z]=v
out["xhbas"]=xhbas; out["xhaci"]=xhaci; out["xbaci"]=xbaci

# ---- D4 reference data from include/param_ref.fh ----
fh = read("include/param_ref.fh")
def data_int(name, width):
    """collect `data name(idx) / v /` integer statements -> dict idx->v"""
    d = {}
    for m in re.finditer(r"data\s+%s\s*\(\s*(\d+)\s*\)\s*/\s*(\d+)\s*/" % name, fh):
        d[int(m.group(1))] = int(m.group(2))
    return d

def data_real_1d(name):
    m = re.search(r"data\s+%s\b(.*?)(?=data\s+\w+\s*\(|!)", fh, re.S)
    vals = {}
    # entries look like: data name(i[,j]) / v [, v] /
    for mm in re.finditer(r"data\s+%s\s*\(([^)]*)\)\s*/\s*([^/]*)/" % name, fh):
        idx = [int(x) for x in re.findall(r"\d+", mm.group(1))]
        vs = [float(x) for x in re.findall(r"-?\d+\.\d+(?:[eEdD][-+]?\d+)?", mm.group(2))]
        vals[tuple(idx)] = vs
    return vals

refn   = data_int("refn", None)
sscale_d = data_real_1d("sscale")
sscale = {k[0]: v[0] for k,v in sscale_d.items() if len(k)==1}
# refsys: `data refsys (iref, iat) / int /`
refsys = {}
for mm in re.finditer(r"data\s+refsys\s*\((\d+)\s*,\s*(\d+)\)\s*/\s*(-?\d+)\s*/", fh):
    refsys[(int(mm.group(1)), int(mm.group(2)))] = int(mm.group(3))
secaiw = data_real_1d("secaiw")          # (23, isys)
refcn  = data_real_1d("refcn")           # (iref, iat)
ascale = data_real_1d("ascale")          # (iref, iat)
hcount = data_real_1d("hcount")          # (iref, iat)
alphaiw= data_real_1d("alphaiw")         # (23, iref, iat) — indices (k, iref, iat)?

d4 = dict(refn={str(k):v for k,v in refn.items()},
          refsys={",".join(map(str,k)):v for k,v in refsys.items()},
          sscale={str(k):v for k,v in sscale.items()},
          secaiw={",".join(map(str,k)):v for k,v in secaiw.items()},
          refcn={",".join(map(str,k)):v for k,v in refcn.items()},
          ascale={",".join(map(str,k)):v for k,v in ascale.items()},
          hcount={",".join(map(str,k)):v for k,v in hcount.items()},
          alphaiw={",".join(map(str,k)):v for k,v in alphaiw.items()})
out["d4"] = d4

json.dump(out, open("data/gfnff_params.json","w"))
print("elements tables OK (103)")
print("rcov[:5] =", [round(x,4) for x in out["rcov"][:5]], "(Bohr expected ~0.585 for H)")
print("d4: refn[1,6,8] =", refn.get(1), refn.get(6), refn.get(8),
      "| alphaiw entries:", len(alphaiw), "| secaiw:", len(secaiw))
print("wrote data/gfnff_params.json")
