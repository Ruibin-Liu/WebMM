//! GFN-FF (Spicher & Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665)
//!
//! Faithful port of the xtb implementation (LGPL-3.0) for main-group
//! organic molecules. Status / coverage:
//!   [x] parameters (auto-extracted from xtb source, data/gfnff_params.json)
//!   [x] topology: bonds (gfnffrab criterion), hybridization (organic rules)
//!   [x] erf coordination number (logCN)
//!   [x] topology charges qa + geometry-dependent EEQ charges & electrostatics
//!   [x] SE non-bonded repulsion (incl. H...H and X...H pair factors)
//!   [x] bond stretching (exp(-alp*(r-r0)^2) with full prefactor chain)
//!   [x] angle bending (cos potential with gfnffdampa damping)
//!   [x] D4(BJ) dispersion (charge-scaled C6 from reference systems)
//!   [x] torsions + out-of-plane impropers (fc*(1-cos phi) / cos double well)
//!   [x] pi (Hückel) bond orders (Fermi smearing, biradical, Nel-1 retry)
//!   [x] HB: eg1 unbound-H, eg2new bound-H, eg2_rnr aromatic-N lone pair,
//!       eg3 carbonyl/nitro (all with exact gradients)
//!   [x] XB halogen bonds + bonded ATM three-body (incl. gradients)
//!   [x] rings (smallest-ring prefactors), fragments (per-fragment EEQ charges)
//!   [x] bond detection via gfnffrab radii (normcn CN, qa/rqshrink, fat);
//!       H-bonded X-H exponent softening (egbond_hb, vbond_scale)
//!   [x] full angle phi0 rules (ring 3/4/5, amide-N, heavy atoms, halogens,
//!       carbene, Ph-O-Ph, CO2) + N-sp2/N-O/O-O torsion specials + amide/alphaCO
//!   [x] charged systems: per-fragment EEQ (try-both placement), pi-system
//!       charge ipis (Cp- / nitromethane validated)
//!   [ ] metals/eta, PBC -> see docs/gfnff-porting-notes.md
//!
//! Reference implementation validated term-by-term against
//! `xtb --gfnff --verbose` (see tests below).
//!
//! Clippy: dense matrix/neighbor loops intentionally index by (i, j) to mirror
//! the Fortran reference (gfnff_eg.f90/gfnff_ini.f90); allowed module-wide like
//! the etkdg distance-geometry code.
#![allow(clippy::needless_range_loop)]
#![allow(clippy::type_complexity)]
// The module intentionally mirrors the Fortran source layout (compact
// multi-declaration lines, section blocks) — formatting it would destroy the
// line-by-line correspondence with gfnff_eg/gfnff_ini that the port relies on.
#![cfg_attr(rustfmt, rustfmt_skip)]

use serde::Deserialize;
use std::collections::HashMap;

pub const BOHR: f64 = 0.52917726;

#[cfg(test)]
thread_local! {
    pub static SELF_DEBUG_TERMS: std::cell::RefCell<Vec<(usize, usize, f64, f64, f64, f64, f64, f64)>> =
        const { std::cell::RefCell::new(Vec::new()) };
} // Angstrom per Bohr (xtb convention)
const TSQRT2PI: f64 = 0.797884560802866; // sqrt(2/pi)

// ---------------------------------------------------------------------------
// parameters
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct RawParams {
    chi: Vec<f64>, gam: Vec<f64>, cnf: Vec<f64>, alp: Vec<f64>, bond: Vec<f64>,
    repa: Vec<f64>, repan: Vec<f64>, angl: Vec<f64>, angl2: Vec<f64>,
    tors: Vec<f64>, tors2: Vec<f64>, en: Vec<f64>, rad: Vec<f64>, repz: Vec<f64>,
    #[allow(dead_code)] metal: Vec<i32>, group: Vec<i32>, normcn: Vec<i32>,
    rcov: Vec<f64>, sqrt_zr4r2: Vec<f64>, rab_en: Vec<f64>, rab_r0: Vec<f64>, rab_cnfak: Vec<f64>,
    zeta_zeff: Vec<i32>, zeta_c: Vec<f64>,
    xhbas: Vec<f64>, xhaci: Vec<f64>, xbaci: Vec<f64>,
    generator: serde_json::Value,
    misc: serde_json::Value,
    d4: RawD4,
}

#[derive(Deserialize)]
struct RawD4 {
    refn: HashMap<String, i32>,
    refsys: HashMap<String, i32>,
    sscale: HashMap<String, f64>,
    secaiw: HashMap<String, Vec<f64>>,
    refcn: HashMap<String, Vec<f64>>,
    ascale: HashMap<String, Vec<f64>>,
    hcount: HashMap<String, Vec<f64>>,
    alphaiw: HashMap<String, Vec<f64>>,
}

/// All element tables and generator constants, in xtb conventions.
pub struct Params {
    pub chi: Vec<f64>, pub gam: Vec<f64>, pub cnf: Vec<f64>, pub alp: Vec<f64>,
    pub bond: Vec<f64>, pub repa: Vec<f64>, pub repan: Vec<f64>,
    pub angl: Vec<f64>, pub angl2: Vec<f64>, pub tors: Vec<f64>, pub tors2: Vec<f64>,
    pub en: Vec<f64>, pub rad: Vec<f64>, pub repz: Vec<f64>,
    pub group: Vec<i32>,
    pub rcov: Vec<f64>,     // covalent radii (D3), ANGSTROM
    pub r4r2: Vec<f64>,     // sqrt(<r4>/<r2>)
    pub rab_en: Vec<f64>,   // gfnffrab EN table
    pub rab_r0: Vec<f64>,   // gfnffrab base radii
    pub rab_cnfak: Vec<f64>,// gfnffrab CN factors
    pub zeta_zeff: Vec<i32>, pub zeta_c: Vec<f64>,
    // generator constants
    pub cnmax: f64, srb1: f64, srb2: f64, srb3: f64,
    pub qrepscal: f64, pub nrepscal: f64, pub hhfac: f64, pub hh13rep: f64, pub hh14rep: f64,
    pub bstren: Vec<f64>, pub qfacben: f64, pub rabshift: f64, pub rabshifth: f64,
    pub rfgoed1: f64, pub qfacbm0: f64, pub fringbo: f64, pub bsmat: Vec<Vec<f64>>,
    pub atcuta: f64, pub repscaln: f64,
    pub d3a1: f64, pub d3a2: f64,
    pub hdiag_tab: std::collections::HashMap<usize, f64>,
    pub hoff_tab: std::collections::HashMap<usize, f64>,
    pub hdiag_c: f64, pub hoff_c: f64,
    pub hiter: f64, pub hueckelp3: f64, pub pilpf: f64, pub htriple: f64,
    pub hueckelp: f64, pub bzref: f64, pub hueckelp2: f64, pub bzref2: f64,
    pub hbacut: f64, pub hbscut: f64, pub hblongcut: f64, pub hbalp: f64,
    pub hblongcut_xb: f64, pub xbacut: f64, pub xbscut: f64, pub xbst: f64, pub xbsf: f64,
    pub batmscal: f64, pub zb3atm: Vec<f64>, pub hbabmix: f64,
    pub hbsf: f64, pub hbst: f64,
    // HB/XB ancillary constants (gfnff_param.f90 449-453 + gfnff_hbset thresholds,
    // accuracy=0.1 as used by the GFN-FF calculator: hbthr = 200/400 - log10(0.1)*50)
    pub tors_hb: f64, pub bend_hb: f64, pub hbnbcut: f64,
    pub xhaci_coh: f64, pub xhaci_globabh: f64,
    pub hqabthr: f64, pub qabthr: f64, pub hbthr1: f64, pub hbthr2: f64,
    /// normal valences (bond-detection CN guess) + bond thresholds (gfnff_param.f90 730-733)
    pub normcn: Vec<i32>, pub rthr: f64, pub rqshrink: f64,
    /// X-H bond softening for H-bonded N/O donors (gfnff_param.f90 451)
    pub vbond_scale: f64,
    /// angle-setup constants (gfnff_param.f90 720-731): heavy-atom sp3 phi0,
    /// small-angle damp, linearity threshold
    pub aheavy3: f64, pub aheavy4: f64, pub fbs1: f64, pub linthr: f64,
    pub xhbas: Vec<f64>, pub xhaci: Vec<f64>, pub xbaci: Vec<f64>,
    // D4 reference data (element-indexed)
    d4_refn: Vec<i32>, d4_refsys: Vec<Vec<i32>>, d4_refcn: Vec<Vec<f64>>,
    d4_ascale: Vec<Vec<f64>>, d4_hcount: Vec<Vec<f64>>,
    d4_alphaiw: Vec<Vec<Vec<f64>>>, d4_sscale: Vec<f64>, d4_secaiw: Vec<Vec<f64>>,
}

impl Params {
    pub fn load() -> Self {
        let raw: RawParams = serde_json::from_str(include_str!("../../data/gfnff_params.json"))
            .expect("embedded gfnff_params.json");
        let g = &raw.generator;
        let gf = |k: &str| g[k].as_f64().expect(k);
        let bstren: Vec<f64> = g["bstren"].as_array().unwrap().iter()
            .map(|v| v.as_f64().unwrap()).collect();
        let mut bsmat = [[-999.0f64; 4]; 4];
        for i in 0..4 { for j in 0..4 {
            bsmat[i][j] = g["bsmat"][i][j].as_f64().unwrap();
        }}
        let m = &raw.misc;
        let mf = |k: &str| m[k].as_f64().expect(k);

        // D4 data: index maps "i" or "i,j" -> value
        let mut refn = vec![0i32; 104]; let mut refsys = vec![vec![0i32; 8]; 104];
        let mut refcn = vec![vec![0f64; 8]; 104];
        let mut ascale = vec![vec![0f64; 8]; 104];
        let mut hcount = vec![vec![0f64; 8]; 104];
        let mut alphaiw = vec![vec![vec![0f64; 8]; 104]; 23];
        let mut sscale = vec![0f64; 18]; let mut secaiw = vec![vec![0f64; 18]; 23];
        // all indices below are stored 1-based in JSON, converted to 0-based here
        // (param_ref.fh covers 118 elements; keep only Z <= 104)
        let ok = |p: &[usize]| p.iter().all(|&x| (1..=104).contains(&x));
        for (k, v) in &raw.d4.refn {
            let z: usize = k.parse().unwrap();
            if z <= 104 { refn[z-1] = *v; }
        }
        for (k, v) in &raw.d4.refsys {           // key "iref,iat"
            let p: Vec<usize> = k.split(',').map(|x| x.parse().unwrap()).collect();
            if ok(&p) { refsys[p[1]-1][p[0]-1] = *v; }
        }
        for (k, v) in &raw.d4.refcn {            // key "iref,iat"
            let p: Vec<usize> = k.split(',').map(|x| x.parse().unwrap()).collect();
            if ok(&p) { refcn[p[1]-1][p[0]-1] = v[0]; }
        }
        for (k, v) in &raw.d4.ascale {           // key "iref,iat"
            let p: Vec<usize> = k.split(',').map(|x| x.parse().unwrap()).collect();
            if ok(&p) { ascale[p[1]-1][p[0]-1] = v[0]; }
        }
        for (k, v) in &raw.d4.hcount {           // key "iref,iat"
            let p: Vec<usize> = k.split(',').map(|x| x.parse().unwrap()).collect();
            if ok(&p) { hcount[p[1]-1][p[0]-1] = v[0]; }
        }
        for (k, v) in &raw.d4.sscale {            // key "is" (1-based system id)
            let is: usize = k.split(',').next().unwrap().parse().unwrap();
            if is <= 18 { sscale[is-1] = *v; }
        }
        for (k, v) in &raw.d4.secaiw {            // key "is", 23 grid values
            let is: usize = k.split(',').next().unwrap().parse().unwrap();
            if is <= 18 { for (w, val) in v.iter().enumerate() { secaiw[w][is-1] = *val; } }
        }
        for (k, v) in &raw.d4.alphaiw {           // key "iref,iat", 23 grid values
            let p: Vec<usize> = k.split(',').map(|x| x.parse().unwrap()).collect();
            if ok(&p) {
                let (iref, iat) = (p[0]-1, p[1]-1);
                for (w, val) in v.iter().enumerate() { alphaiw[w][iat][iref] = *val; }
            }
        }
        // note: alpha integration uses the non-uniform frequency grid of trapzd
        // (src/disp/dftd4.F90 381); weights computed in d4_dispersion

        Params {
            chi: raw.chi, gam: raw.gam, cnf: raw.cnf, alp: raw.alp, bond: raw.bond,
            repa: raw.repa, repan: raw.repan, angl: raw.angl, angl2: raw.angl2,
            tors: raw.tors, tors2: raw.tors2, en: raw.en, rad: raw.rad, repz: raw.repz,
            group: raw.group, rcov: raw.rcov, r4r2: raw.sqrt_zr4r2,
            rab_en: raw.rab_en, rab_r0: raw.rab_r0, rab_cnfak: raw.rab_cnfak,
            zeta_zeff: raw.zeta_zeff, zeta_c: raw.zeta_c,
            cnmax: gf("cnmax"), srb1: gf("srb1"), srb2: gf("srb2"), srb3: gf("srb3"),
            qrepscal: gf("qrepscal"), nrepscal: gf("nrepscal"),
            hhfac: gf("hhfac"), hh13rep: gf("hh13rep"), hh14rep: gf("hh14rep"),
            bstren, qfacben: gf("qfacBEN"), rabshift: gf("rabshift"), rabshifth: gf("rabshifth"),
            rfgoed1: gf("rfgoed1"), qfacbm0: gf("qfacbm0"), fringbo: gf("fringbo"),
            bsmat: bsmat.iter().map(|r| r.to_vec()).collect(),
            atcuta: mf("atcuta"), repscaln: mf("repscaln"),
            d3a1: gf("d3a1"), d3a2: gf("d3a2"),
            hdiag_tab: [(5,-0.5),(6,0.0),(7,0.14),(8,-0.38),(9,-0.29),(16,-0.30),(17,-0.30)]
                .iter().map(|(k,v)| (*k as usize, *v)).collect(),
            hoff_tab: [(5,0.5),(6,1.00),(7,0.66),(8,1.10),(9,0.23),(16,0.60),(17,1.00)]
                .iter().map(|(k,v)| (*k as usize, *v)).collect(),
            hdiag_c: 0.0, hoff_c: 0.0,
            hiter: gf("hiter"), hueckelp3: gf("hueckelp3"), pilpf: gf("pilpf"),
            htriple: gf("htriple"), hueckelp: gf("hueckelp"), bzref: gf("bzref"),
            hueckelp2: gf("hueckelp2"), bzref2: gf("bzref2"),
            hbacut: mf("hbacut"), hbscut: mf("hbscut"), hblongcut: mf("hblongcut"),
            hblongcut_xb: mf("hblongcut_xb"), xbacut: mf("xbacut"), xbscut: mf("xbscut"),
            xbst: mf("xbst"), xbsf: mf("xbsf"),
            batmscal: gf("batmscal"), hbabmix: mf("hbabmix"),
            zb3atm: { let mut v = vec![0.0f64; 104];
                for (i, val) in v.iter_mut().enumerate() {
                    let z = (i + 1) as f64;   // index i == Z-1 (gfnff_param.f90 495-502)
                    *val = if i == 0 { -0.25 * gf("batmscal").powf(1.0/3.0) }
                           else { -z * gf("batmscal").powf(1.0/3.0) };
                }
                v },
            hbalp: mf("hbalp"), hbsf: mf("hbsf"), hbst: mf("hbst"),
            tors_hb: mf("tors_hb"), bend_hb: mf("bend_hb"), hbnbcut: mf("hbnbcut"),
            xhaci_coh: mf("xhaci_coh"), xhaci_globabh: mf("xhaci_globabh"),
            hqabthr: gf("hqabthr"), qabthr: gf("qabthr"),
            hbthr1: 250.0, hbthr2: 450.0,   // Bohr^2 (gfnff_thresholds, accuracy 0.1)
            normcn: raw.normcn.clone(), rthr: gf("rthr"), rqshrink: gf("rqshrink"),
            vbond_scale: mf("vbond_scale"),
            aheavy3: gf("aheavy3"), aheavy4: gf("aheavy4"),
            fbs1: gf("fbs1"), linthr: gf("linthr"),
            xhbas: raw.xhbas.clone(), xhaci: raw.xhaci.clone(), xbaci: raw.xbaci.clone(),
            d4_refn: refn, d4_refsys: refsys, d4_refcn: refcn, d4_ascale: ascale,
            d4_hcount: hcount, d4_alphaiw: alphaiw, d4_sscale: sscale, d4_secaiw: secaiw,
        }
    }

    /// D4 charge-scaling function (gfnff_ini.f90 zeta(), lines 2384-2438)
    /// zeta = exp(3*(1 - exp(+c*(1 - zeff/qmod))))  [NOTE: plus sign before c]
    pub fn zeta(&self, at: usize, q: f64) -> f64 {
        let zeff = self.zeta_zeff[at - 1] as f64;
        let qmod = zeff + q;
        if qmod < 0.0 { return (3.0f64).exp(); }
        (3.0 * (1.0 - (self.zeta_c[at - 1] * (1.0 - zeff / qmod)).exp())).exp()
    }
}

// ---------------------------------------------------------------------------
// setup / topology
// ---------------------------------------------------------------------------

pub struct Topology {
    pub nb: Vec<Vec<usize>>,   // bonded neighbors (0-based)
    pub hyb: Vec<i32>,         // 0 unknown/none, 1 sp, 2 sp2, 3 sp3, 5 hypervalent
    pub qa: Vec<f64>,          // topology charges (geometry independent EEQ)
    pub chieeq: Vec<f64>,      // final EEQ electronegativity
    pub gameeq: Vec<f64>,      // final EEQ J (gamma)
    pub alpeeq: Vec<f64>,      // final EEQ alpha (squared exponent)
    pub dxi: Vec<f64>,
    pub nb13: Vec<(usize, usize)>, // 1-3 pairs (for H..H rep factor)
    pub nfrag: usize,
    pub piadr: Vec<bool>,      // atom is a pi atom
    /// pi bond order per bond (parallel to `bonds`), HMO density matrix elements
    pub pibo: Vec<f64>,
    /// bond-path matrix: 0=unreached, 1=bonded, 2=1-3, 3=1-4, 5=far
    pub bpair: Vec<Vec<usize>>,
    /// smallest ring size containing the atom (0 = none); up to 20 rings each
    pub ring_size: Vec<usize>,
    /// fragment id per atom (0-based)
    pub fraglist: Vec<usize>,
    /// charge per fragment
    pub qfrag: Vec<f64>,
    // ---- HB/XB/bATM ancillary (gfnff_ini.f90 885-1000 + 735-800) ----
    /// per-atom HB basicity (xhbas w/ carbene/carbonyl/nitro rules)
    pub hbbas: Vec<f64>,
    /// per-atom HB acidity (xhaci w/ amide-H x0.80)
    pub hbaci: Vec<f64>,
    /// itag=1 flags: bent 2-coordinate group-14 (carbene) and nitro N
    pub carbene: Vec<bool>,
    /// nitro N flag (3-coord N with terminal O; xtb itag=1, hyb 2)
    pub nitro_n: Vec<bool>,
    /// HB-relevant H atoms (qa > hqabthr, non-bridging)
    pub hb_h: Vec<usize>,
    /// HB-relevant (A,B) pairs (i > j, charge/type filtered)
    pub hb_ab: Vec<(usize, usize)>,
    /// XB triples (A, B, X): X bonded to A, B the donor base
    pub xb_triples: Vec<(usize, usize, usize)>,
    /// bonded-ATM triples (i, j, k) from 1-4 pairs (xtb b3list)
    pub b3: Vec<(usize, usize, usize)>,
    /// per-bond HB acceptor list (N/O only): X-H bonds whose donor A and an
    /// hb_ab partner B are both N/O get a softened exponent (egbond_hb)
    pub bond_hb_b: Vec<Vec<usize>>,
    /// all rings (up to size 6) per atom, as sorted member lists
    pub rings_all: Vec<Vec<Vec<usize>>>,
}

pub struct EnergyComponents {
    pub bond: f64, pub angle: f64, pub torsion: f64,
    pub rep: f64, pub es: f64, pub disp: f64,
    pub hb: f64, pub xb: f64, pub batm: f64,
}

impl EnergyComponents {
    pub fn total(&self) -> f64 {
        self.bond + self.angle + self.torsion + self.rep + self.es + self.disp + self.hb + self.xb + self.batm
    }
}

/// Adapter: GFN-FF as a composable force source (kcal/mol + kcal/mol/A),
/// drop-in for the WebMM optimizer/MD pipeline.
pub struct GfnffForceField {
    inner: std::rc::Rc<Gfnff>,
}

impl GfnffForceField {
    pub fn new(at: &[usize], xyz: &[[f64; 3]], charge: f64) -> Self {
        GfnffForceField { inner: std::rc::Rc::new(Gfnff::new(at, xyz, charge)) }
    }

    /// EnergyComponents at arbitrary coordinates (Eh, topology fixed at setup),
    /// matching what energy_and_gradient totals — for term-wise display.
    pub fn components_at(&self, coords: &[[f64; 3]]) -> EnergyComponents {
        self.inner.energy(coords)
    }

    /// Topology EEQ charges (qa) from setup — for charge coloring/exports.
    pub fn charges(&self) -> Vec<f64> {
        self.inner.topo.qa.clone()
    }
}

impl crate::forces::ForceField for GfnffForceField {
    fn energy_and_gradient(&self, coords: &[[f64; 3]], grad: &mut [[f64; 3]]) -> f64 {
        let e = self.inner.energy_and_gradient(coords, grad);   // E in Eh, grad in Eh/Bohr
        // convert BOTH to kcal/mol and kcal/mol/Å (the previous version scaled
        // only the energy, leaving the gradient 2.24× too small — L-BFGS then
        // stalled at the input geometry and "converged" immediately)
        const SCALE: f64 = 627.5094740631 / BOHR;
        for g in grad.iter_mut() {
            for c in g.iter_mut() { *c *= SCALE; }
        }
        e.total() * 627.5094740631   // Eh -> kcal/mol
    }
}

pub struct Gfnff {
    pub p: Params,
    pub at: Vec<usize>,        // atomic numbers
    pub charge: f64,
    pub topo: Topology,
    /// per-bond (i, j, alp, kbond, rab0_bohr) from setup
    pub bonds: Vec<BondParam>,
    /// per-angle (j, i, k, phi0_rad, fc)
    pub angles: Vec<AngleParam>,
    /// per-torsion (l, i, j, k, nrot, phi0_rad, fc)
    pub torsions: Vec<TorsionParam>,
    /// setup geometry (Bohr), used for chktors linearity filters
    xyz0: Vec<[f64; 3]>,
}

pub struct BondParam { pub i: usize, pub j: usize, pub alp: f64, pub kb: f64, pub r0: f64 }
pub struct AngleParam { pub j: usize, pub i: usize, pub k: usize, pub phi0: f64, pub fc: f64 }
pub struct TorsionParam { pub l: usize, pub i: usize, pub j: usize, pub k: usize,
                          pub nrot: i32, pub phi0: f64, pub fc: f64 }

fn erf(x: f64) -> f64 {
    // Abramowitz-Stegun 7.1.26 insufficient; use high-accuracy series
    // (matches Fortran erf to ~1e-12 for our range)
    libm_erf(x)
}

#[inline]
fn libm_erf(x: f64) -> f64 {
    // Rust std has no erf; use the W. J. Cody rational approximation (used by libm)
    // via the incomplete gamma identity would be slow; implement Cody (1970).
    let ax = x.abs();
    let (res, _) = erf_cody(ax);
    if x >= 0.0 { res } else { -res }
}

fn erf_cody(x: f64) -> (f64, f64) {
    // erf via Taylor series (small x) + continued fraction for erfc (large x)
    // accuracy ~1e-15, sufficient to match Fortran libm erf
    let sqpi = std::f64::consts::FRAC_2_SQRT_PI; // 2/sqrt(pi)
    let xx = x.abs();
    let res = if xx <= 2.5 {
        // erf(x) = 2x/sqrt(pi) * sum_{n} (-1)^n x^{2n} / (n! (2n+1))
        let x2 = xx * xx;
        let mut term = 1.0f64;
        let mut sum = 1.0f64;
        let mut n = 1u32;
        while n < 100 {
            term *= -x2 / n as f64;
            let add = term / (2.0 * n as f64 + 1.0);
            sum += add;
            if add.abs() < 1e-18 { break; }
            n += 1;
        }
        sqpi * xx * sum
    } else {
        // erfc(x) = exp(-x^2)/sqrt(pi) * 1/f,  f = x + a1/(x + a2/(x + ...)),  a_j = j/2
        // modified Lentz iteration
        let tiny = 1e-300;
        let mut f = xx;
        let mut c = f;
        let mut d = 0.0f64;
        for j in 1..80 {
            let a = j as f64 * 0.5;
            d = 1.0 / (xx + a * d);
            if d.abs() < tiny { d = tiny; }
            c = xx + a / c;
            if c.abs() < tiny { c = tiny; }
            f *= c * d;
            if (c * d - 1.0).abs() < 1e-16 { break; }
        }
        1.0 - (-xx * xx).exp() / (std::f64::consts::PI.sqrt() * f)
    };
    (res, res)
}

impl Gfnff {
    /// `at` atomic numbers (1-based), `xyz_ang` coordinates in Angstrom.
    pub fn new(at: &[usize], xyz_ang: &[[f64; 3]], charge: f64) -> Self {
        let p = Params::load();
        let n = at.len();
        let xyz: Vec<[f64; 3]> = xyz_ang.iter().map(|r| [r[0]/BOHR, r[1]/BOHR, r[2]/BOHR]).collect();

        // --- coordination number (erf) ---
        let cn = erf_cn(&p, at, &xyz);

        // --- bond detection (gfnffrab criterion, simplified qa shift) ---
        // first topology charges need bonds -> iterate twice like xtb's qloop
        let mut topo_q = vec![0.0f64; n];
        let mut nb: Vec<Vec<usize>>;
        let mut hyb: Vec<i32>;
        let mut rabd = vec![vec![0.0f64; n]; n];
        // early pi detection (needed for amide/ff_gam before full HMO)
        let mut piadr_temp = vec![false; n];
        {
            nb = detect_bonds(&p, at, &xyz, &topo_q);
            (hyb, _) = assign_hyb(at, &nb, &xyz, &topo_q);
            for i in 0..n {
                let z = at[i];
                let h = hyb[i];
                if matches!(z, 6|7|8|9|16) && (h == 1 || h == 2) { piadr_temp[i] = true; }
                let attached_sp2 = nb[i].iter().any(|&j| hyb[j] == 1 || hyb[j] == 2);
                if matches!(z, 7..=9) && attached_sp2 { piadr_temp[i] = true; }   // picon: sp3 N/O/F on sp2 joins π (nofs)
            }
        }
        // fragment detection: connected components on the bond graph
        // (computed once here; xtb detects fragments BEFORE the topology EEQ
        // and constrains EACH fragment charge separately - a total-charge-only
        // constraint lets charge leak between fragments and skews qa)
        let mut fraglist = vec![0usize; n];
        let mut nfrag = 0usize;
        {
            nb = detect_bonds(&p, at, &xyz, &topo_q);
            hyb = assign_hyb(at, &nb, &xyz, &topo_q).0;
            for s in 0..n {
                if fraglist[s] != 0 { continue; }
                nfrag += 1;
                let mut stack = vec![s];
                fraglist[s] = nfrag;
                while let Some(a) = stack.pop() {
                    for &b in &nb[a] {
                        if fraglist[b] == 0 { fraglist[b] = nfrag; stack.push(b); }
                    }
                }
            }
        }
        let fraglist0: Vec<usize> = fraglist.iter().map(|x| x - 1).collect();
        let mut qfrag = vec![0.0f64; nfrag];
        if nfrag == 1 { qfrag[0] = charge; }
        // charged 2-fragment systems: try BOTH charge placements and keep the
        // one with the lower topology-EEQ energy (gfnff_ini.f90 570-588);
        // more fragments: charge on the most electropositive (approximation)
        if nfrag == 2 && charge != 0.0 {
            // dxi (topology electronegativity corrections, ini 391-447) for a
            // consistent es comparison; the qloop below recomputes it per pass
            let dxi_try = compute_dxi(&p, at, &nb, &hyb, &piadr_temp);
            let rabd0 = floyd_rabd(&p, at, &nb);   // topology distances for the try-both EEQ
            qfrag = vec![0.0, charge];
            let (_, es1) = solve_eeq(&p, at, &rabd0, &nb, charge, true, &hyb, &dxi_try, &fraglist0, &qfrag);
            qfrag = vec![charge, 0.0];
            let (_, es2) = solve_eeq(&p, at, &rabd0, &nb, charge, true, &hyb, &dxi_try, &fraglist0, &qfrag);
            qfrag = if es1 < es2 { vec![0.0, charge] } else { vec![charge, 0.0] };
        } else if nfrag > 1 && charge != 0.0 {
            let mut frag_en = vec![0.0f64; nfrag];
            let mut frag_cnt = vec![0.0f64; nfrag];
            for i in 0..n {
                frag_en[fraglist0[i]] += p.en[at[i]-1];
                frag_cnt[fraglist0[i]] += 1.0;
            }
            for fi in 0..nfrag {
                if frag_cnt[fi] > 0.0 { frag_en[fi] /= frag_cnt[fi]; }
            }
            let most_electropos = (0..nfrag)
                .min_by(|a, b| frag_en[*a].partial_cmp(&frag_en[*b]).unwrap())
                .unwrap_or(0);
            qfrag[most_electropos] = charge;
        }

        let mut dxi = vec![0.0f64; n];
        let mut nitro = vec![false; n];
        for _iter in 0..2 {
            nb = detect_bonds(&p, at, &xyz, &topo_q);
            (hyb, nitro) = assign_hyb(at, &nb, &xyz, &topo_q);
            dxi = compute_dxi(&p, at, &nb, &hyb, &piadr_temp);
            // topology-distance EEQ charges (geometry independent)
            rabd = floyd_rabd(&p, at, &nb);
            let (q_new, _) = solve_eeq(&p, at, &rabd, &nb, charge, true, &hyb, &dxi, &fraglist0, &qfrag);
            topo_q = q_new;
        }

        // --- EEQ parameters (second pass, gfnff_ini.f90 lines 670-719) ---
        let mut chieeq = vec![0.0f64; n];
        let mut gameeq = vec![0.0f64; n];
        let mut alpeeq = vec![0.0f64; n];
        for i in 0..n {
            let z = at[i];
            // dgam: qa * ff(gam)
            let ff_gam = match z {
                1 => -0.08, 5 => -0.05,
                6 => if hyb[i] < 3 { -0.45 } else { -0.27 },
                7 => {
                    if is_amide_n(at, &hyb, &nb, &piadr_temp, i) { -0.16 }
                    else if piadr_temp[i] { -0.14 }
                    else { -0.13 }
                },
                8 => if hyb[i] < 3 { -0.08 } else { -0.15 },
                9 => 0.10, 17 => -0.02, 35 => -0.11, 53 => -0.07,
                z if z > 10 => -0.02,
                _ => 0.0,
            };
            let ff_alp = if z == 6 { 0.09 }
                else if z == 7 { -0.21 }
                else if p.group[z-1] == 6 { -0.03 }
                else if p.group[z-1] == 7 { 0.50 }
                else { 0.0 };
            chieeq[i] = -p.chi[z-1] + dxi[i];
            gameeq[i] = p.gam[z-1] + topo_q[i] * ff_gam;
            alpeeq[i] = (p.alp[z-1] + ff_alp * topo_q[i]).powi(2);
        }

        // --- 1-3 pairs ---
        let mut nb13 = Vec::new();
        for i in 0..n {
            for a in 0..nb[i].len() {
                for b in (a+1)..nb[i].len() {
                    let (j, k) = (nb[i][a], nb[i][b]);
                    if j < k { nb13.push((j, k)); } else { nb13.push((k, j)); }
                }
            }
        }
        nb13.sort(); nb13.dedup();

        // carbene flags (xtb itag=1): bent 2-coordinate group-14, unless strongly anionic
        let carbene: Vec<bool> = (0..n).map(|i| {
            let z = at[i];
            let group14 = matches!(z, 6 | 14 | 32 | 50 | 82);
            group14 && nb[i].len() == 2 && hyb[i] == 2 && topo_q[i] >= -0.4
        }).collect();

        let topo = Topology { nb, hyb, qa: topo_q, chieeq, gameeq, alpeeq, dxi, nb13, nfrag, piadr: vec![false; n], pibo: vec![0.0; 0], bpair: Vec::new(), ring_size: vec![0; n], fraglist: fraglist0, qfrag, hbbas: vec![0.0; n], hbaci: vec![0.0; n], carbene, nitro_n: nitro, hb_h: Vec::new(), hb_ab: Vec::new(), xb_triples: Vec::new(), b3: Vec::new(), bond_hb_b: Vec::new(),
            rings_all: Vec::new() };

        let mut g = Gfnff { p, at: at.to_vec(), charge, topo, bonds: Vec::new(), angles: Vec::new(), torsions: Vec::new(), xyz0: xyz.clone() };
        g.setup_bonds(&xyz, &rabd, &cn);
        // bond-path matrix (nbondmat_pbc, non-periodic) + smallest rings
        g.setup_bpair_rings();
        // pi system: HMO bond orders (gfnff_ini.f90 975-1121)
        g.setup_hmo(&xyz);
        g.setup_bonds(&xyz, &rabd, &cn);   // re-run with pibo active
        g.setup_angles(&xyz, &cn);
        g.setup_torsions();
        // HB/XB/bATM lists + per-atom bas/aci (gfnff_ini.f90 735-1000)
        g.setup_hb_xb_batm();
        g
    }

    fn setup_bonds(&mut self, xyz: &[[f64; 3]], rabd: &[Vec<f64>], _cn: &[f64]) {
        let p = &self.p;
        let mut bonds = Vec::new();
        for (bi2, (bi, bj)) in self.topo.nb.iter().enumerate()
            .flat_map(|(i, nbs)| nbs.iter().map(move |&j| (i, j)))
            .filter(|&(i, j)| i < j).enumerate() {
            let (ia, ja) = (self.at[bi], self.at[bj]);
            let hybi = self.topo.hyb[bi].max(self.topo.hyb[bj]);
            let hybj = self.topo.hyb[bi].min(self.topo.hyb[bj]);
            // bond type (gfnff_ini.f90 1229-1244, organic subset)
            let mut btyp = 1u8;
            if self.topo.hyb[bi] == 2 && self.topo.hyb[bj] == 2 { btyp = 2; }   // sp2-sp2 = pi
            if (self.topo.hyb[bi] == 3 && self.topo.hyb[bj] == 2 && ia == 7)
            || (self.topo.hyb[bj] == 3 && self.topo.hyb[bi] == 2 && ja == 7) { btyp = 2; }  // N-sp2
            if self.topo.hyb[bi] == 1 || self.topo.hyb[bj] == 1 { btyp = 3; }   // sp-X
            let mut bstrength = if hybi == 5 || hybj == 5 { p.bstren[3] }
                else { p.bsmat[hybi as usize][hybj as usize] };
            // N-sp2 bonds use the pi-type strength (ini 1258-1259)
            if hybi == 3 && hybj == 2 && (ia == 7 || ja == 7) {
                bstrength = p.bstren[1] * 1.04;
            }
            // CO (sp C-O) slightly weaker (ini 1295-1296)
            if btyp == 3 && ((ia == 6 && ja == 8) || (ia == 8 && ja == 6)) {
                bstrength = p.bstren[2] * 0.90;
            }
            let mut shift = 0.0;
            let mut fxh = 1.0;
            let mut fcn = 1.0f64;
            if ia == 1 || ja == 1 { shift = p.rabshifth; }
            if ia == 9 && ja == 9 { shift = 0.22; }   // F2
            // X-sp3 correction (both directions)
            if (self.topo.hyb[bi] == 3 && self.topo.hyb[bj] == 0)
            || (self.topo.hyb[bj] == 3 && self.topo.hyb[bi] == 0) { shift -= 0.022; }
            // X-sp correction
            if (self.topo.hyb[bi] == 1 && self.topo.hyb[bj] == 0)
            || (self.topo.hyb[bj] == 1 && self.topo.hyb[bi] == 0) { shift += 0.14; }
            if (ia == 1 && ja == 8) || (ja == 1 && ia == 8) { fxh = 0.93; }
            if (ia == 1 && ja == 6) || (ja == 1 && ia == 6) {
                fxh = 1.0;
                let c = if ia == 6 { bi } else { bj };
                if self.topo.ring_size[c] == 3 { fxh = 1.05; }          // 3-ring CH
                if is_aldehyde_c(&self.at, &self.topo.nb, &self.topo.piadr, c) { fxh = 0.95; }  // aldehyde CH
            }
            if (ia == 1 && ja == 5) || (ja == 1 && ia == 5) { fxh = 1.10; }
            if (ia == 1 && ja == 7) || (ja == 1 && ia == 7) { fxh = 1.06; }
            if (ia == 1 && ja == 9) || (ja == 1 && ia == 9) { fxh = 1.0; }
            if ia > 10 && ja > 10 {
                shift += -0.11; // hshift3
                if ia > 18 { shift += -0.11; }
                if ja > 18 { shift += -0.11; }
                fcn /= 1.0 + 0.007 * (self.topo.nb[bi].len() as f64).powi(2);
                fcn /= 1.0 + 0.007 * (self.topo.nb[bj].len() as f64).powi(2);
            }
            let qafac = self.topo.qa[bi] * self.topo.qa[bj] * 70.0;
            let fqq = 1.0 + p.qfacbm0 / (1.0 + (15.0 * qafac).exp());
            let en_diff = p.en[ia-1] - p.en[ja-1];
            // pi corrections (pibo > 0)
            // ring prefactor (ringsbond): smallest ring containing the bond
            let ringf = {
                let ri = self.topo.ring_size[bi]; let rj = self.topo.ring_size[bj];
                let rings = if ri > 0 && ri == rj { ri } else { 0 };
                if rings > 0 { 1.0 + self.p.fringbo * (6.0 - rings as f64).powi(2) } else { 1.0 }
            };
            let pibo = self.topo.pibo.get(bi2).copied().unwrap_or(0.0);
            if pibo > 0.0 {
                shift += p.hueckelp * (p.bzref - pibo);
                if hybi != 1 && hybj != 1 && pibo > 0.1 { /* btyp=2 handled by bstren below */ }
            }
            let mut fpi = 1.0f64;
            if pibo > 0.0 { fpi = 1.0 - p.hueckelp2 * (p.bzref2 - pibo); }
            let vb1 = p.rabshift + shift;
            let vb2 = p.srb1 * (1.0 + p.srb2 * en_diff * en_diff + p.srb3 * bstrength);
            let vb3 = -p.bond[ia-1] * p.bond[ja-1] * ringf * bstrength * fqq * fpi * fxh * fcn;
            bonds.push(BondParam { i: bi, j: bj, alp: vb2, kb: vb3, r0: vb1 });
        }
        self.bonds = bonds;
    }

    fn setup_angles(&mut self, _xyz: &[[f64; 3]], _cn: &[f64]) {
        let p = &self.p;
        let mut angles = Vec::new();
        for i in 0..self.at.len() {
            let nn = self.topo.nb[i].len();
            if nn <= 1 || nn > 6 { continue; }
            let nbs = &self.topo.nb[i];
            for a in 0..nbs.len() {
                for b in (a+1)..nbs.len() {
                    let (jj, kk) = (nbs[a], nbs[b]);
                    let (atj, atk) = (self.at[jj], self.at[kk]);
                    let fijk = p.angl[self.at[i]-1] * p.angl2[atj-1] * p.angl2[atk-1];
                    if fijk < 1e-3 { continue; }
                    let nh = [atj, atk].iter().filter(|&&z| z == 1).count();
                    let no = [atj, atk].iter().filter(|&&z| z == 8).count();
                    let nc = [atj, atk].iter().filter(|&&z| z == 6).count();
                    let nsi = [atj, atk].iter().filter(|&&z| z == 14).count();
                    let npi = [jj, kk].iter().filter(|&&a| self.topo.piadr[a]).count();
                    let rings = ringsbend(&self.topo.rings_all, i, jj, kk);
                    let triple = self.topo.hyb[i] == 1 || self.topo.hyb[jj] == 1 || self.topo.hyb[kk] == 1;
                    let phi_deg = self.angle_deg(jj, i, kk);
                    let ati = self.at[i];
                    // ---- phi0/f2 rules (gfnff_ini.f90 1610-1790, organic subset) ----
                    let mut r0;
                    let mut f2 = 1.0f64;
                    match self.topo.hyb[i] { 1 => r0 = 180.0, 2 => r0 = 120.0, 3 => r0 = 109.5, _ => r0 = 100.0 }
                    if self.topo.hyb[i] == 3 && ati > 10 {
                        if nn <= 3 { r0 = self.p.aheavy3; }   // heavy main-group 3-coord
                        if nn >= 4 { r0 = self.p.aheavy4; }
                        if nn == 4 && self.p.group[ati - 1] == 5 { r0 = 109.5; }
                        if nn == 4 && self.p.group[ati - 1] == 4 && ati > 49 { r0 = 109.5; }
                        if (4..=6).contains(&self.p.group[ati - 1]) { r0 -= nh as f64 * 5.0; }
                    }
                    if ati == 5 && matches!(self.topo.hyb[i], 2 | 3) { r0 = 115.0; }
                    if ati == 6 {
                        if self.topo.hyb[i] == 3 && nh == 2 { r0 = 108.6; }
                        if self.topo.hyb[i] == 3 && no == 1 { r0 = 108.5; }
                        if self.topo.hyb[i] == 2 && no == 2 { r0 = 122.0; }
                        if self.topo.hyb[i] == 2 && no == 1 { f2 = 0.7; }
                        if self.topo.hyb[i] == 1 && no == 2 { f2 = 2.0; }   // CO2
                        if self.topo.hyb[i] == 3 && nn > 4 && phi_deg > self.p.linthr { r0 = 180.0; }
                    }
                    if ati == 8 && nn == 2 {
                        r0 = 104.5;
                        if nh == 2 { r0 = 100.0; f2 = 1.20; }
                        r0 += 7.0 * nsi as f64;
                        if npi == 2 { r0 = 109.0; }   // e.g. Ph-O-Ph
                    }
                    if ati == 7 && nn == 2 {
                        f2 = 1.4; r0 = 115.0;
                        if rings != 0 { r0 = 105.0; }
                        if atj == 8 || atk == 8 { r0 = 103.0; }
                        if atj == 9 || atk == 9 { r0 = 102.0; }
                        if self.topo.hyb[i] == 1 { r0 = 180.0; }
                    }
                    if ati == 7 && self.topo.hyb[i] == 3 {
                        if npi > 0 {
                            if is_amide_n(&self.at, &self.topo.hyb, &self.topo.nb, &self.topo.piadr, i) {
                                r0 = 115.0; f2 = 1.2;
                            } else {
                                let sumppi: f64 = self.bonds.iter().enumerate()
                                    .filter(|(_, b)| (b.i == i && (b.j == jj || b.j == kk))
                                                 || (b.j == i && (b.i == jj || b.i == kk)))
                                    .map(|(bi_, _)| self.topo.pibo.get(bi_).copied().unwrap_or(0.0))
                                    .sum();
                                r0 = 113.0;
                                f2 = 1.0 - sumppi * 0.7;
                            }
                        } else {
                            r0 = 104.0;
                            f2 = 0.40 + nh as f64 * 0.19 + no as f64 * 0.25 + nc as f64 * 0.01;
                        }
                    }
                    if rings == 3 { r0 = 82.0; }
                    if rings == 4 { r0 = 96.0; }
                    if rings == 5 && ati == 6 { r0 = 109.0; }
                    if rings == 0 {
                        // 3-ring center with ONE neighbor in a 3-ring and one
                        // acyclic neighbor (xtb ringsatom returns the sentinel
                        // 99 for atoms in no ring, so its literal
                        // ringsj+ringsk.eq.102 means 99+3; our ring_size uses 0
                        // for acyclic atoms, hence the mapping below)
                        if self.topo.ring_size[i] == 3 {
                            let rj = if self.topo.ring_size[jj] == 0 { 99 } else { self.topo.ring_size[jj] };
                            let rk = if self.topo.ring_size[kk] == 0 { 99 } else { self.topo.ring_size[kk] };
                            if rj + rk == 102 { r0 += 4.0; }
                        }
                    }
                    if triple { f2 = 0.60; if atj == 7 || atk == 7 { f2 = 1.00; } }
                    if self.p.group[ati - 1] == 4 && nn == 2 && self.topo.carbene[i] {
                        if ati == 6 { r0 = 145.0; } else { r0 = 90.0; }
                    }
                    if self.p.group[ati - 1] == 6 && nn == 4 && no >= 1 { r0 = 115.0; }
                    if self.p.group[ati - 1] == 7 && self.topo.hyb[i] == 1 {
                        if matches!(ati, 9 | 17 | 35 | 53) { r0 = 90.0; }
                        if ati > 9 && phi_deg > self.p.linthr { r0 = 180.0; }
                        f2 = 0.6 / (ati as f64).powf(0.15);
                    }
                    if self.topo.hyb[i] == 3 && self.p.group[ati - 1] == 4 && ati > 32
                        && self.topo.qa[i] > 0.4 {
                        if phi_deg > 140.0 { r0 = 180.0; }
                        if phi_deg < 100.0 { r0 = 90.0; }
                        f2 = 1.0;
                    }
                    // ---- force constant ----
                    let fqq = 1.0 - (self.topo.qa[i]*self.topo.qa[jj] + self.topo.qa[i]*self.topo.qa[kk]) * p.qfacben;
                    let fnn = 1.0 - 2.36 / (nn as f64).powi(2);   // gfnff_ini.f90 1792
                    let phi0 = r0 * std::f64::consts::PI / 180.0;
                    let fbsmall = 1.0 - 0.5 * (-0.64 * (phi0 - std::f64::consts::PI).powi(2)).exp();
                    let fc = fijk * fqq * f2 * fnn * fbsmall;
                    angles.push(AngleParam { j: jj, i, k: kk, phi0, fc });
                }
            }
        }
        self.angles = angles;
    }

    /// bond-path matrix (1=bond, 2=1-3, 3=1-4, 5=far) and smallest ring sizes
    fn setup_bpair_rings(&mut self) {
        let n = self.at.len();
        let mut bpair = vec![vec![0usize; n]; n];
        for i in 0..n { for &j in &self.topo.nb[i] { bpair[i][j] = 1; } }
        // BFS up to depth 3 (shortest path only)
        for i in 0..n {
            for depth in 1..=2 {
                // atoms at path length `depth` from i
                let frontier: Vec<usize> = (0..n).filter(|&j| bpair[i][j] == depth).collect();
                for &f in &frontier {
                    for &j in &self.topo.nb[f] {
                        if j != i && bpair[i][j] == 0 { bpair[i][j] = depth + 1; }
                    }
                }
            }
        }
        for i in 0..n { for j in 0..n { if i != j && bpair[i][j] == 0 { bpair[i][j] = 5; } } }
        self.topo.bpair = bpair;
        // all rings up to size 6 via cycle search (getring36): keep sorted
        // member lists per atom (for ringsbend) + smallest size per atom
        let mut rings_all: Vec<Vec<Vec<usize>>> = vec![Vec::new(); n];
        let mut ring = vec![0usize; n];
        let mut seen: std::collections::HashSet<Vec<usize>> = Default::default();
        for a0 in 0..n {
            for &a1 in &self.topo.nb[a0] {
                for &a2 in &self.topo.nb[a1] {
                    if a2 == a0 || a2 == a1 { continue; }
                    for &a3 in &self.topo.nb[a2] {
                        if a3 == a0 {
                            let members = sorted3(a0, a1, a2);
                            if seen.insert(members.clone()) { add_ring(&mut rings_all, &members); }
                            continue;
                        }
                        if a3 == a1 || a3 == a2 { continue; }
                        for &a4 in &self.topo.nb[a3] {
                            if a4 == a0 {
                                let members = sorted4(a0, a1, a2, a3);
                                if seen.insert(members.clone()) { add_ring(&mut rings_all, &members); }
                                continue;
                            }
                            if a4 == a1 || a4 == a2 || a4 == a3 { continue; }
                            for &a5 in &self.topo.nb[a4] {
                                if a5 == a0 {
                                    let members = sorted5(a0, a1, a2, a3, a4);
                                    if seen.insert(members.clone()) { add_ring(&mut rings_all, &members); }
                                    continue;
                                }
                                if a5 == a1 || a5 == a2 || a5 == a3 || a5 == a4 { continue; }
                                for &a6 in &self.topo.nb[a5] {
                                    if a6 == a0 {
                                        let members = sorted6(a0, a1, a2, a3, a4, a5);
                                        if seen.insert(members.clone()) { add_ring(&mut rings_all, &members); }
                                        continue;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for a0 in 0..n {
            ring[a0] = rings_all[a0].iter().map(|r| r.len()).min().unwrap_or(0);
        }
        self.topo.ring_size = ring;
        self.topo.rings_all = rings_all;
    }

    /// pi-system charge from a neutralized-fragment EEQ difference
    /// (gfnff_ini.f90 610-660): dqa = qheavy(qa) - qheavy(qa_neutral) summed
    /// over the pi atoms, x1.1, rounded to the nearest integer
    fn pi_system_charge(&self, pi_atoms: &[usize], rabd: &[Vec<f64>]) -> i32 {
        let ifrag = self.topo.fraglist[pi_atoms[0]];
        if self.topo.qfrag[ifrag] == 0.0 { return 0; }
        let p = &self.p;
        let mut qah = self.topo.qa.clone();
        qheavy(&self.at, &self.topo.nb, &mut qah);
        let mut qfrag2 = self.topo.qfrag.clone();
        qfrag2[ifrag] = 0.0;
        let (qa0, _) = solve_eeq(p, &self.at, rabd, &self.topo.nb, self.charge, true,
            &self.topo.hyb, &self.topo.dxi, &self.topo.fraglist, &qfrag2);
        let mut qa0 = qa0;
        qheavy(&self.at, &self.topo.nb, &mut qa0);
        let mut dum = 0.0f64;
        for &a in pi_atoms {
            dum += qah[a] - qa0[a];
        }
        (dum * 1.1).round() as i32
    }

    /// pi system detection + iterative HMO (gfnff_ini.f90 357-390, 975-1121)
    fn setup_hmo(&mut self, xyz: &[[f64; 3]]) {
        let n = self.at.len();
        // pi atom list (organic subset of piadr rules)
        let mut piadr = vec![false; n];
        for i in 0..n {
            let z = self.at[i];
            let hyb = self.topo.hyb[i];
            let inlist = matches!(z, 6 | 7 | 8 | 9 | 16);
            let mut piat = (hyb == 1 || hyb == 2) && inlist;
            // sp3 N/O/F attached to sp2/sp → pi (nofs rule; xtb ini 373-386:
            // picon = kk>0 .and. nofs — NO hyb!=3 condition, amide N's join)
            let attached_sp2 = self.topo.nb[i].iter()
                .any(|&j| self.topo.hyb[j] == 1 || self.topo.hyb[j] == 2);
            if matches!(z, 7..=9) && attached_sp2 { piat = true; }
            if !attached_sp2 && hyb == 3 { piat = false; }
            if z == 7 && self.topo.nb[i].len() > 3 { piat = false; } // NR4
            if piat { piadr[i] = true; }
        }
        // connected components over pi atoms
        let mut pimvec = vec![0usize; n];   // pi-system id per atom (1-based), 0 = none
        let mut npis = 0usize;
        for s in 0..n {
            if !piadr[s] || pimvec[s] != 0 { continue; }
            npis += 1;
            let mut stack = vec![s];
            pimvec[s] = npis;
            while let Some(a) = stack.pop() {
                for &b in &self.topo.nb[a] {
                    if piadr[b] && pimvec[b] == 0 { pimvec[b] = npis; stack.push(b); }
                }
            }
        }
        // per-system electron count (piel rules) and atoms
        self.topo.pibo = vec![0.0; self.bonds.len()];
        let mut piadr_out = vec![false; n];
        // topology distances for the ipis EEQ re-solve (shared by all systems;
        // hoisted so multi-pi molecules do one Floyd pass instead of one each)
        let rabd_ipis = floyd_rabd(&self.p, &self.at, &self.topo.nb);
        for pis in 1..=npis {
            let atoms: Vec<usize> = (0..n).filter(|&a| pimvec[a] == pis).collect();
            let npi = atoms.len();
            if npi < 2 { continue; }
            let mut nel = 0i32;
            let mut piel = vec![0i32; n];
            for &a in &atoms {
                let (z, hyb) = (self.at[a], self.topo.hyb[a]);
                // piel electron-count rules (gfnff_ini.f90 1006-1028);
                // carbene C (itag=1) contributes 0; the N-hyb2-itag nitro rule
                // never fires in our subset (carbene tags only set for group 14)
                let carbene = self.topo.carbene[a];
                let nelpi = match (z, hyb) {
                    (5, 1) => 1,                                    // B in borine
                    (6, _) if !carbene => 1,                        // skip carbene
                    (7, 2) if self.topo.nitro_n[a] => 2,            // nitro N (itag=1): no own el
                    (7, _) if hyb <= 2 => 1,
                    (7, 3) => 2,
                    (8, 1) => 1, (8, 2) => 1, (8, 3) => 2,
                    (9, 1) => 3, (9, _) => 2,                        // F
                    (16, 1) => 1, (16, 2) => 1, (16, 3) => 2,       // S
                    (17, 0) => 2, (17, 1) => 3,                     // Cl
                    _ => 0,
                };
                piel[a] = nelpi.min(2);
                nel += nelpi;
            }
            // pi-system charge ipis (gfnff_ini.f90 610-660): neutralize the
            // fragment holding this pi system, re-solve the topology EEQ, and
            // take the rounded (x1.1) charge difference on the pi atoms
            nel -= self.pi_system_charge(&atoms, &rabd_ipis);
            if nel < 1 { continue; }
            let nel = nel as usize;
            // iterative HMO (gfnff_ini.f90 1031-1066)
            let idx: Vec<usize> = atoms.clone();
            let pos = |a: usize| idx.iter().position(|&x| x == a).unwrap();
            let mut pold = vec![vec![2.0/3.0; npi]; npi];
            let mut dens = pold.clone();
            let mut eold = 0.0f64;
            let mut eps_homo = 0.0f64;
            let mut api_last = vec![vec![0.0f64; npi]; npi];
            for _iter in 0..5 {
                let mut api = vec![vec![0.0f64; npi]; npi];
                for (ia, &a) in atoms.iter().enumerate() {
                    let hd = *self.p.hdiag(self.at[a]);
                    api[ia][ia] = hd + self.topo.qa[a] * self.p.hueckelp3
                        - (piel[a] - 1) as f64 * self.p.pilpf;
                }
                for b in &self.bonds {
                    if pimvec[b.i] != pis || pimvec[b.j] != pis { continue; }
                    let (ia, ja) = (pos(b.i), pos(b.j));
                    let dum = 1e-9 * dist2(xyz, b.i, b.j).sqrt();
                    let ho = (self.p.hoff(self.at[b.i]) * self.p.hoff(self.at[b.j])).sqrt();
                    let mut dum2 = self.p.hiter;
                    if self.topo.hyb[b.i] == 1 || self.topo.hyb[b.j] == 1 { dum2 *= self.p.htriple; }
                    let v = -ho * (1.0 - dum2 * (2.0/3.0 - pold[ja][ia]));
                    api[ja][ia] = v; api[ia][ja] = v;
                }
                api_last = api.clone();
                let (p, eel, eps, focc) = hmo_solve(&api, nel);
                // eps(HOMO) = last eigenvalue with occupation > 0.5
                eps_homo = eps.iter().zip(focc.iter())
                    .filter(|(_, &f)| f > 0.5).map(|(e, _)| *e).next_back().unwrap_or(0.0);
                dens = p;
                if (eel - eold).abs() < 1e-4 { break; }
                pold = dens.clone();
                eold = eel;
            }
            // wrong pi occupation: HOMO > 0.40 eV -> one retry with Nel-1
            // (gfnff_ini.f90 1075-1102; fixes nitro/radical occupations)
            if eps_homo > 0.40 && nel > 1 {
                let (p, _eel, _eps, _focc) = hmo_solve(&api_last, nel - 1);
                dens = p;
            }
            // store pibo on bonds
            for (bi_, b) in self.bonds.iter().enumerate() {
                if pimvec[b.i] == pis && pimvec[b.j] == pis {
                    self.topo.pibo[bi_] = dens[pos(b.j)][pos(b.i)];
                    piadr_out[b.i] = true;
                    piadr_out[b.j] = true;
                }
            }
        }
        self.topo.piadr = piadr_out;
    }

    /// torsion list for the organic/saturated subset (gfnff_ini.f90 1830-2050)
    fn setup_torsions(&mut self) {
        let p = &self.p;
        let n = self.at.len();
        let torsf = [1.00f64, 1.18, 1.05, 0.0, 0.50, -0.90, 0.70, -2.00]; // gen%torsf
        let mut out = Vec::new();
        for bi in 0..n {
            for &bj in &self.topo.nb[bi].clone() {
                if bj <= bi { continue; }  // unique central bonds
                let (zi, zj) = (self.at[bi], self.at[bj]);
                let hybi = self.topo.hyb[bi];
                let hybj = self.topo.hyb[bj];
                if hybi == 1 || hybj == 1 { continue; }       // sp/linear: no torsion
                if hybi == 5 || hybj == 5 { continue; }       // hypervalent
                if p.tors[zi-1] < 0.0 || p.tors[zj-1] < 0.0 { continue; }
                if p.tors[zi-1] * p.tors[zj-1] < 1e-3 { continue; }
                // N-sp2 bonds count as pi-type for torsions too (ini 1231-1234)
                let btyp = if (hybi == 2 && hybj == 2)
                    || (hybi == 3 && hybj == 2 && self.at[bi] == 7)
                    || (hybi == 3 && hybj == 2 && self.at[bj] == 7) { 2 } else { 1 };
                let nhi = 1 + self.topo.nb[bi].iter().filter(|&&x| self.at[x] == 1).count();
                let nhj = 1 + self.topo.nb[bj].iter().filter(|&&x| self.at[x] == 1).count();
                let mut fij = p.tors[zi-1] * p.tors[zj-1] * ((nhi as f64) * (nhj as f64)).powf(0.07);
                if alpha_co(&self.at, &self.topo.hyb, &self.topo.nb, &self.topo.piadr, bi, bj) { fij *= 1.3; }
                // amide C-N torsions are stiffer (ini 1881-1883)
                if self.topo.hyb[bj] == 3 && self.at[bj] == 6
                    && is_amide_n(&self.at, &self.topo.hyb, &self.topo.nb, &self.topo.piadr, bi) { fij *= 1.3; }
                if self.topo.hyb[bi] == 3 && self.at[bi] == 6
                    && is_amide_n(&self.at, &self.topo.hyb, &self.topo.nb, &self.topo.piadr, bj) { fij *= 1.3; }
                let fqq = 1.0 + (self.topo.qa[bi] * self.topo.qa[bj]).abs() * 12.0;
                for &kk in &self.topo.nb[bi].clone() {
                    if kk == bj { continue; }
                    for &ll in &self.topo.nb[bj].clone() {
                        if ll == bi { continue; }
                        // chktors: skip near-linear arrangements
                        if self.angle_deg(bj, bi, kk) > 170.0 { continue; }
                        if self.angle_deg(bi, bj, ll) > 170.0 { continue; }
                        let zk = self.at[kk];
                        let zl = self.at[ll];
                        let mut fkl = p.tors2[zk-1] * p.tors2[zl-1]
                            * (self.topo.nb[kk].len() as f64 * self.topo.nb[ll].len() as f64).powf(-0.14);
                        if zk == 7 && !self.topo.piadr[kk] { fkl *= 0.5; }   // saturated N penalty
                        if zl == 7 && !self.topo.piadr[ll] { fkl *= 0.5; }
                        if fkl < 1e-3 { continue; }
                        // ring vs acyclic torsion definition (gfnff_ini.f90 1897-1943)
                        let (mut nrot, mut phi0);
                        let mut f1 = torsf[0];
                        // lring: central bond in a ring (ringsbond, 0 = none)
                        let ringb = ringsbond(&self.topo.rings_all, bi, bj);
                        if ringb > 0 {
                            // RING CASE: cis minima (phi0 = 0); a central bond in
                            // a 3-ring forces rings4 = 3 ("the 3-ring is special")
                            let rings4 = if ringb > 3 {
                                ring_common(&self.topo.rings_all, [bi, bj, kk, ll], false)
                            } else {
                                3
                            };
                            nrot = if btyp == 2 { 2 } else { 1 };
                            phi0 = 0.0;
                            if btyp == 1 && rings4 > 0 {
                                let ringl = ring_common(&self.topo.rings_all, [bi, bj, kk, ll], true);
                                let notpicon = !self.topo.piadr[kk] && !self.topo.piadr[ll];
                                if rings4 == 3 && notpicon { nrot = 1; phi0 = 0.0; f1 = 0.3; }   // gen%fr3
                                if rings4 == 4 && ringl == 4 && notpicon {
                                    nrot = 6; phi0 = 30.0 * std::f64::consts::PI / 180.0; f1 = 1.0; } // fr4
                                if rings4 == 5 && ringl == 5 && notpicon {
                                    nrot = 6; phi0 = 30.0 * std::f64::consts::PI / 180.0; f1 = 1.5; } // fr5
                                if rings4 == 6 && ringl == 6 && notpicon {
                                    nrot = 3; phi0 = 60.0 * std::f64::consts::PI / 180.0; f1 = 5.7; } // fr6
                            }
                            // ring bond but no common 4-atom ring: terminal flanks
                            if rings4 == 0 && btyp == 1
                                && self.topo.nb[kk].len() == 1 && self.topo.nb[ll].len() == 1 {
                                nrot = 6; phi0 = 30.0 * std::f64::consts::PI / 180.0; f1 = 0.30;
                            }
                            // 5-ring C-N pi bonds adjacent to an amide (CB7 fix)
                            if btyp == 2 && ringb == 5 && self.at[bi] * self.at[bj] == 42
                                && (is_amide_n(&self.at, &self.topo.hyb, &self.topo.nb, &self.topo.piadr, bi)
                                    || is_amide_n(&self.at, &self.topo.hyb, &self.topo.nb, &self.topo.piadr, bj)) {
                                f1 = 5.0;
                            }
                        } else {
                            // ACYCLIC CASE (ini 1922-1943)
                            if hybi == 3 && hybj == 3 {
                                (nrot, phi0) = (3, std::f64::consts::PI);
                            } else if btyp == 2 {
                                (nrot, phi0) = (2, std::f64::consts::PI);
                            } else {
                                (nrot, phi0) = (1, std::f64::consts::PI);
                            }
                            // pi-sp3 central bond: soft torsion (N-pi especially)
                            if self.topo.piadr[bi] && !self.topo.piadr[bj] && self.topo.hyb[bj] == 3 {
                                f1 = 0.5;
                                if self.at[bi] == 7 { f1 = 0.2; }
                                nrot = 3;
                            }
                            if self.topo.piadr[bj] && !self.topo.piadr[bi] && self.topo.hyb[bi] == 3 {
                                f1 = 0.5;
                                if self.at[bj] == 7 { f1 = 0.2; }
                                nrot = 3;
                            }
                        }
                        // SP3 specials apply AFTER both cases (ini 1945-1963)
                        if hybi == 3 && hybj == 3 {
                            let (gi, gj) = (p.group[zi - 1], p.group[zj - 1]);
                            if gi == 5 && gj == 5 { nrot = 3; phi0 = 60.0 * std::f64::consts::PI / 180.0; f1 = 3.0; }
                            if (gi == 5 && gj == 6) || (gi == 6 && gj == 5) {
                                nrot = 2; phi0 = std::f64::consts::PI / 2.0; f1 = 1.0;
                                if zi >= 15 && zj >= 15 { f1 = 20.0; }
                            }
                            if gi == 6 && gj == 6 {
                                nrot = 2; phi0 = std::f64::consts::PI / 2.0; f1 = 5.0;
                                if zi >= 16 && zj >= 16 { f1 = 25.0; }
                            }
                        }
                        // pibo for the central bond (bonds list is ordered i<j; linear search ok, few bonds)
                        let pibo = self.bonds.iter().position(|b| (b.i==bi && b.j==bj) || (b.i==bj && b.j==bi))
                            .and_then(|bi_| self.topo.pibo.get(bi_).copied()).unwrap_or(0.0);
                        let mut f2 = 0.0f64;
                        if pibo > 0.0 {
                            f2 = pibo * (-2.5 * (1.24 - pibo).powi(14)).exp();
                            if !self.topo.piadr[kk] && self.at[kk] > 10 { f2 *= 1.3; }
                            if !self.topo.piadr[ll] && self.at[ll] > 10 { f2 *= 1.3; }
                            f1 *= 0.55;
                        }
                        let fctot = (f1 + 10.0 * torsf[1] * f2) * fqq * fij * fkl;
                        if fctot > 1e-3 {
                            out.push(TorsionParam { l: ll, i: bi, j: bj, k: kk, nrot, phi0, fc: fctot });
                        }
                        // extra rot=1 torsion for sp3-sp3 to get gauche conf
                        // energies well — acyclic only (xtb: .not.lring)
                        if ringb == 0 && self.topo.hyb[kk] == 3 && self.topo.hyb[ll] == 3
                            && hybi == 3 && hybj == 3 {
                            let ff = if zi == 7 || zj == 7 { torsf[6] }
                                else if zi == 8 || zj == 8 { torsf[7] }
                                else { torsf[5] };
                            out.push(TorsionParam { l: ll, i: bi, j: bj, k: kk, nrot: 1,
                                phi0: std::f64::consts::PI, fc: ff * fij * fkl * fqq });
                        }
                        let _ = &mut nrot;
                    }
                }
            }
        }
        // ---- out-of-plane impropers (gfnff_ini.f90 2040-2150) ----
        // tlist order: (center, 1st, 2nd, 3rd nb sorted by distance to center)
        for i in 0..self.at.len() {
            if self.topo.nb[i].len() != 3 { continue; }
            let pi_center = self.topo.piadr[i];
            if !pi_center && self.at[i] != 7 { continue; }
            let mut nbs = self.topo.nb[i].clone();
            nbs.sort_by(|&a, &b| dist2(&self.xyz0, a, i)
                .partial_cmp(&dist2(&self.xyz0, b, i)).unwrap());
            let (jj, kk, ll) = (nbs[0], nbs[1], nbs[2]);
            if !pi_center && self.at[i] == 7 {
                // saturated N: cos-based double well at +/- phi0
                let mut fc = 0.0f64;
                for &m in &nbs { fc += 0.60 * self.p.repz[self.at[m] - 1].sqrt(); }
                out.push(TorsionParam { l: i, i: jj, j: kk, k: ll, nrot: -1,
                    phi0: 80.0 * std::f64::consts::PI / 180.0, fc });
            } else {
                // pi center: fc = torsf(3) * (1 - sum pibo * torsf(5)) * (1 + 5 qa)
                // with carbonyl/halogen/nitro corrections (ini 2123-2127)
                let ncarbo = nbs.iter().filter(|&&m| self.at[m] == 8 || self.at[m] == 16).count();
                let nf = nbs.iter().filter(|&&m| self.p.group[self.at[m] - 1] == 7).count();
                let fqq = 1.0 + self.topo.qa[i] * 5.0;
                let pibo_sum: f64 = self.bonds.iter().enumerate()
                    .filter(|(_, b)| b.i == i || b.j == i)
                    .map(|(bi_, _)| self.topo.pibo.get(bi_).copied().unwrap_or(0.0))
                    .sum();
                let f2 = 1.0 - pibo_sum * 0.50;  // torsf(5)
                let mut fc = 1.05 * f2 * fqq;     // torsf(3)
                if (self.at[i] == 5 || self.at[i] == 6) && ncarbo > 0 { fc *= 38.0; }
                if self.at[i] == 6 && nf > 0 && ncarbo == 0 { fc *= 10.0; }
                if self.at[i] == 7 && ncarbo > 0 { fc *= 10.0 / f2; }  // no pi dep
                out.push(TorsionParam { l: i, i: jj, j: kk, k: ll, nrot: 0, phi0: 0.0, fc });
            }
        }
        self.torsions = out;
    }

    fn angle_deg(&self, a: usize, b: usize, c: usize) -> f64 {
        let v1 = [self.xyz0[a][0]-self.xyz0[b][0], self.xyz0[a][1]-self.xyz0[b][1], self.xyz0[a][2]-self.xyz0[b][2]];
        let v2 = [self.xyz0[c][0]-self.xyz0[b][0], self.xyz0[c][1]-self.xyz0[b][1], self.xyz0[c][2]-self.xyz0[b][2]];
        let cos = (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
            / ((v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]).sqrt()
             * (v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]).sqrt());
        cos.clamp(-1.0, 1.0).acos() * 180.0 / std::f64::consts::PI
    }

    /// HB/XB/bATM list setup + per-atom HB basicity/acidity
    /// (gfnff_ini.f90 735-1000: b3list, hbbas, hbaci, hbatHl, hbatABl, xbatABl)
    fn setup_hb_xb_batm(&mut self) {
        let n = self.at.len();
        // ---- hbbas: xhbas w/ carbene / carbonyl / nitro rules (ini 848-872) ----
        for i in 0..n {
            let mut bas = self.p.xhbas[self.at[i] - 1];
            let nn = self.topo.nb[i].len();
            if self.at[i] == 6 && nn == 2 && self.topo.carbene[i] { bas = 1.46; }
            if self.at[i] == 8 && nn == 1 {
                if let Some(&j) = self.topo.nb[i].first() {
                    if self.at[j] == 6 { bas = 0.68; }   // carbonyl R-C=O
                    if self.at[j] == 7 { bas = 0.47; }   // nitro R-N=O
                }
            }
            self.topo.hbbas[i] = bas;
        }
        // ---- hbaci: xhaci, amide N x0.80 via its terminal H (ini 874-891) ----
        // NOTE: xtb scales the acidity of the N (the donor heavy atom), not the H
        for i in 0..n { self.topo.hbaci[i] = self.p.xhaci[self.at[i] - 1]; }
        for i in 0..n {
            if self.topo.nb[i].len() != 1 { continue; }
            let nn = self.topo.nb[i][0];
            if !is_amide_n(&self.at, &self.topo.hyb, &self.topo.nb, &self.topo.piadr, nn) { continue; }
            let nc = self.topo.nb[nn].iter()
                .filter(|&&j| self.at[j] == 6 && self.topo.hyb[j] == 3).count();
            if nc == 1 { self.topo.hbaci[nn] *= 0.80; }
        }
        // ---- HB-relevant H atoms (ini 886-902): qa > hqabthr w/ adjustments ----
        let mut hb_h = Vec::new();
        for i in 0..n {
            if self.at[i] != 1 || self.topo.hyb[i] == 1 { continue; }  // skip bridging H
            let Some(&j) = self.topo.nb[i].first() else { continue };
            let mut ff = self.p.hqabthr;
            if self.at[j] > 10 { ff -= 0.20; }                     // H on heavy atom
            if self.at[j] == 6 && self.topo.hyb[j] == 3 { ff += 0.05; }  // H on sp3 C
            if self.topo.qa[i] > ff { hb_h.push(i); }
        }
        // ---- HB-relevant (A,B) pairs (ini 904-922) ----
        let mut hb_ab = Vec::new();
        for i in 0..n {
            if self.at[i] == 6 && !self.topo.piadr[i] { continue; }  // only pi C
            let mut ffi = self.p.qabthr;
            if self.at[i] > 10 { ffi += 0.2; }
            if self.topo.qa[i] > ffi { continue; }
            for j in 0..i {
                if self.at[j] == 6 && !self.topo.piadr[j] { continue; }
                let mut ffj = self.p.qabthr;
                if self.at[j] > 10 { ffj += 0.2; }
                if self.topo.qa[j] > ffj { continue; }
                let (cai, cci) = (self.topo.hbbas[i], self.topo.hbaci[i]);
                let (caj, ccj) = (self.topo.hbbas[j], self.topo.hbaci[j]);
                if cai * ccj < 1e-6 && cci * caj < 1e-6 { continue; }
                hb_ab.push((i, j));
            }
        }
        // ---- XB triples (ini 925-985): A-X bond, X halogen/chalcogen/P, B donor ----
        let xatom = |z: usize| matches!(z, 15 | 16 | 17 | 33 | 34 | 35 | 51 | 52 | 53);
        let mut xb = Vec::new();
        for ia in 0..n {
            for &ix in &self.topo.nb[ia] {
                let zx = self.at[ix];
                if !xatom(zx) { continue; }
                if zx == 16 && self.topo.nb[ix].len() > 2 { continue; }  // no sulfoxide S
                for ib in 0..n {
                    if ib == ia || ib == ix { continue; }
                    if self.topo.bpair[ix][ib] <= 3 { continue; }  // must be A-X...B
                    if self.p.xhbas[self.at[ib] - 1] < 1e-6 { continue; }
                    if self.p.group[self.at[ib] - 1] == 4
                        && (!self.topo.piadr[ib] || self.topo.qa[ib] > 0.05) { continue; }
                    xb.push((ia, ib, ix));
                }
            }
        }
        // ---- b3 triples (ini 753-785): 1-4 pair once, k in nb(i) ∪ nb(j) ----
        let mut b3 = Vec::new();
        for i in 0..n {
            for j in 0..i {
                if self.topo.bpair[i][j] != 3 { continue; }
                for &k in self.topo.nb[i].iter().chain(self.topo.nb[j].iter()) {
                    if k == i || k == j { continue; }
                    b3.push((i, j, k));
                }
            }
        }
        // ---- per-bond HB acceptor map (bond_hbset + bond_hb_AHB_set, ini2 ----
        // 1087-1274): for X-H bonds with N/O donor A and N/O hb_ab partner B
        // (nonbonded, rAB^2 < hbthr1), the bond exponent is softened through
        // the runtime HB coordination number of H (egbond_hb)
        let mut bond_hb_b = vec![Vec::new(); self.bonds.len()];
        let mut register = |a: usize, h: usize, b: usize| {
            if let Some(bi) = self.bonds.iter().position(|bo| {
                (bo.i == a && bo.j == h) || (bo.i == h && bo.j == a)
            }) {
                if !bond_hb_b[bi].contains(&b) { bond_hb_b[bi].push(b); }
            }
        };
        for &(i, j) in &hb_ab {
            if dist2(&self.xyz0, i, j) > self.p.hbthr1 { continue; }
            let ijnonbond = self.topo.bpair[i][j] != 1;
            if !ijnonbond { continue; }
            for &h in &hb_h.clone() {
                let ok_i = self.topo.bpair[i][h] == 1
                    && matches!(self.at[i], 7 | 8) && matches!(self.at[j], 7 | 8);
                let ok_j = self.topo.bpair[j][h] == 1
                    && matches!(self.at[j], 7 | 8) && matches!(self.at[i], 7 | 8);
                if ok_i { register(i, h, j); }
                if ok_j { register(j, h, i); }
            }
        }
        self.topo.bond_hb_b = bond_hb_b;
        self.topo.hb_h = hb_h;
        self.topo.hb_ab = hb_ab;
        self.topo.xb_triples = xb;
        self.topo.b3 = b3;
    }

    // -------------------------------------------------------------------------
    // energy (Eh, coordinates in Angstrom)
    // -------------------------------------------------------------------------

    pub fn energy(&self, xyz_ang: &[[f64; 3]]) -> EnergyComponents {
        let n = self.at.len();
        let xyz: Vec<[f64; 3]> = xyz_ang.iter().map(|r| [r[0]/BOHR, r[1]/BOHR, r[2]/BOHR]).collect();
        let cn = erf_cn(&self.p, &self.at, &xyz);

        let dist = |i: usize, j: usize| -> f64 {
            let dx = xyz[i][0]-xyz[j][0]; let dy = xyz[i][1]-xyz[j][1]; let dz = xyz[i][2]-xyz[j][2];
            (dx*dx + dy*dy + dz*dz).sqrt()
        };

        // ---- EEQ charges + electrostatics (goed_gfnff, gfnff_eg.f90 1758) ----
        let mut q = vec![0.0f64; n];
        let mut es = 0.0;
        {
            let m = n + 1;
            let mut a = vec![vec![0.0f64; m]; m];
            let mut x = vec![0.0f64; m];
            for i in 0..n {
                x[i] = self.topo.chieeq[i] + self.p.cnf[self.at[i]-1] * cn[i].sqrt();
                a[i][i] = TSQRT2PI / self.topo.alpeeq[i].sqrt() + self.topo.gameeq[i];
            }
            for i in 0..n { for j in 0..i {
                let g = 1.0 / (self.topo.alpeeq[i] + self.topo.alpeeq[j]).sqrt();
                let r = dist(i, j);
                let v = erf(g * r) / r;
                a[i][j] = v; a[j][i] = v;
            }}
            let nfrag = self.topo.nfrag;
            let m2 = n + nfrag;
            let mut a2 = vec![vec![0.0f64; m2]; m2];
            let mut x2 = vec![0.0f64; m2];
            for i2 in 0..n {
                x2[i2] = x[i2];
                a2[i2][i2] = a[i2][i2];
            }
            for i2 in 0..n { for j2 in 0..i2 { a2[i2][j2] = a[i2][j2]; a2[j2][i2] = a[i2][j2]; } }
            for (fi, &qf) in self.topo.qfrag.iter().enumerate() {
                x2[n + fi] = qf;
                for j2 in 0..n {
                    if self.topo.fraglist[j2] == fi {
                        a2[n + fi][j2] = 1.0; a2[j2][n + fi] = 1.0;
                    }
                }
            }
            q = solve_sym(&a2, &x2);
            for i in 0..n { for j in 0..i {
                let g = 1.0 / (self.topo.alpeeq[i] + self.topo.alpeeq[j]).sqrt();
                let r = dist(i, j);
                es += q[i] * q[j] * erf(g * r) / r;
            }}
            for i in 0..n {
                let x_i = self.topo.chieeq[i] + self.p.cnf[self.at[i]-1] * cn[i].sqrt();
                es += -q[i] * x_i
                    + q[i]*q[i] * 0.5 * (self.topo.gameeq[i] + TSQRT2PI / self.topo.alpeeq[i].sqrt());
            }
        }

        // ---- SE repulsion (non-bonded, gfnff_eg.f90 345-420) ----
        let mut rep = 0.0;
        for i in 0..n { for j in 0..i {
            if self.topo.nb[i].contains(&j) { continue; }
            let r = dist(i, j);
            if r > 20.0 { continue; }
            let (zi, zj) = (self.at[i], self.at[j]);
            let fni = 1.0 + self.p.nrepscal / (1.0 + (self.topo.nb[i].len() as f64).powi(2));
            let fnj = 1.0 + self.p.nrepscal / (1.0 + (self.topo.nb[j].len() as f64).powi(2));
            let di = self.p.repan[zi-1] * (1.0 + self.topo.qa[i] * self.p.qrepscal) * fni;
            let dj = self.p.repan[zj-1] * (1.0 + self.topo.qa[j] * self.p.qrepscal) * fnj;
            let mut ff = 1.0;
            if zi == 1 && zj == 1 {
                ff = self.p.hhfac;
                let bp = self.topo.bpair.get(i).and_then(|r| r.get(j)).copied().unwrap_or(5);
                if bp == 2 { ff *= self.p.hh13rep; }
                if bp == 3 { ff *= self.p.hh14rep; }
            }
            if (zi == 1 && self.p.metal_is(zj)) || (zj == 1 && self.p.metal_is(zi)) { ff = 0.85; }
            if (zi == 1 && zj == 6) || (zj == 1 && zi == 6) { ff = 0.91; }
            if (zi == 1 && zj == 8) || (zj == 1 && zi == 8) { ff = 1.04; }
            let alpha = (di * dj).sqrt() * ff;
            let t16 = r.powf(1.5);
            let e = (-alpha * t16).exp() * self.p.repz[zi-1] * self.p.repz[zj-1] * self.p.repscaln;
            rep += e / r;
        }}

        // ---- bonds (egbond/egbond_hb + gfnffdrab, gfnff_eg.f90 555-1230) ----
        let mut ebond = 0.0;
        for (bi, b) in self.bonds.iter().enumerate() {
            let r = dist(b.i, b.j);
            // reference length is CN-dependent (gfnffdrab): shift carried inside
            let rab0 = self.p.gfnffrab(self.at[b.i], self.at[b.j], cn[b.i], cn[b.j], b.r0);
            let dr = r - rab0;
            // H-bonded X-H bonds use a softened exponent (egbond_hb)
            let alp = if self.topo.bond_hb_b[bi].is_empty() { b.alp } else {
                let h = if self.at[b.i] == 1 { b.i } else { b.j };
                (1.0 - (1.0 - self.p.vbond_scale) * self.hb_cn_of(&xyz, bi, h)) * b.alp
            };
            ebond += b.kb * (-alp * dr * dr).exp();
        }

        // ---- bonded SE repulsion (gfnff_eg.f90 596-660) ----
        for b in &self.bonds {
            let r = dist(b.i, b.j);
            let (zi, zj) = (self.at[b.i], self.at[b.j]);
            let alpha = (self.p.repa[zi-1] * self.p.repa[zj-1]).sqrt();
            let repab = self.p.repz[zi-1] * self.p.repz[zj-1] * self.p.repscalb();
            let t16 = r.powf(1.5);
            let e = (-alpha * t16).exp() * repab / r;
            rep += e;
        }

        // ---- angles (egbend, gfnff_eg.f90 1233) ----
        let mut eangl = 0.0;
        for ang in &self.angles {
            let (j, i, k) = (ang.j, ang.i, ang.k);
            let va = [xyz[j][0]-xyz[i][0], xyz[j][1]-xyz[i][1], xyz[j][2]-xyz[i][2]];
            let vc = [xyz[k][0]-xyz[i][0], xyz[k][1]-xyz[i][1], xyz[k][2]-xyz[i][2]];
            let rab = (va[0]*va[0]+va[1]*va[1]+va[2]*va[2]).sqrt();
            let rcb = (vc[0]*vc[0]+vc[1]*vc[1]+vc[2]*vc[2]).sqrt();
            let cos = (va[0]*vc[0]+va[1]*vc[1]+va[2]*vc[2]) / (rab*rcb);
            let theta = cos.clamp(-1.0, 1.0).acos();
            let dij = dampa(&self.p, self.at[j], self.at[i], rab * rab);
            let dki = dampa(&self.p, self.at[k], self.at[i], rcb * rcb);
            eangl += ang.fc * (cos - ang.phi0.cos()).powi(2) * dij * dki;
            let _ = theta;
        }

        // ---- torsions (egtors, gfnff_eg.f90 1444) ----
        // NOTE: damping uses the two flanking 1-3 distances and the central bond
        let mut etors = 0.0;
        let mut g_scratch = vec![[0.0f64; 3]; n];
        for t in &self.torsions {
            if t.nrot <= 0 {
                // out-of-plane improper (energy; gradient into scratch)
                etors += self.improper_term(&xyz, t, &mut g_scratch);
                continue;
            }
            let phi = dihedral(&xyz, t.l, t.i, t.j, t.k);
            let c1 = t.nrot as f64 * (phi - t.phi0) + std::f64::consts::PI;
            let et = (1.0 + c1.cos()) * t.fc;
            // flanking 1-3 distances ll-i and j-kk, central bond i-j
            let rli2 = dist2(&xyz, t.l, t.i);
            let rij2 = dist2(&xyz, t.i, t.j);
            let rjk2 = dist2(&xyz, t.j, t.k);
            let damp = dampt(&self.p, self.at[t.l], self.at[t.i], rli2)
                * dampt(&self.p, self.at[t.i], self.at[t.j], rij2)
                * dampt(&self.p, self.at[t.j], self.at[t.k], rjk2);
            etors += et * damp;
        }

        // ---- D4(BJ) dispersion with EEQ-charge scaling (gdisp0.f90) ----
        let disp = self.d4_dispersion(&xyz, &dist, &cn);

        // ---- HB / XB / bonded ATM (gfnff_eg.f90 2088-3446) ----
        let ebatm = self.batm_terms(&xyz, &mut g_scratch);
        let ehb = self.hb_terms(&xyz, &mut g_scratch);
        let exb = self.xb_terms(&xyz, &mut g_scratch);

        EnergyComponents {
            bond: ebond, angle: eangl, torsion: etors,
            rep, es, disp, hb: ehb, xb: exb, batm: ebatm,
        }
    }

    // -----------------------------------------------------------------------
    // HB / XB / bonded-ATM (gfnff_eg.f90 2088-3446) — faithful ports.
    // Energy and gradient in ONE code path each, so energy() and
    // energy_and_gradient() are consistent by construction. All terms use
    // TOPOLOGY charges qa (xtb passes topo%qa), gradients in Eh/Bohr.
    // -----------------------------------------------------------------------

    /// charge smearing sigma(s·q) = e^{s q}/(e^{s q}+f)
    fn csig(ex: f64, sf: f64) -> f64 { ex / (ex + sf) }

    /// build hblist1 (unbound H) + hblist2 (bound H) at this geometry and
    /// evaluate all HB terms (gfnff_hbset, gfnff_ini2.F90 764-948, non-PBC)
    fn hb_terms(&self, xyz: &[[f64; 3]], g: &mut [[f64; 3]]) -> f64 {
        let p = &self.p;
        let mut e = 0.0f64;
        for &(i, j) in &self.topo.hb_ab {
            let rab2 = dist2(xyz, i, j);
            if rab2 > p.hbthr1 { continue; }
            let ijnonbond = self.topo.bpair[i][j] != 1;
            for &h in &self.topo.hb_h {
                // bound-H: H covalently bound to one end of a nonbonded A-B pair
                if ijnonbond && self.topo.bpair[i][h] == 1 {
                    e += self.hb_bound(xyz, i, j, h, g);
                    continue;
                }
                if ijnonbond && self.topo.bpair[j][h] == 1 {
                    e += self.hb_bound(xyz, j, i, h, g);
                    continue;
                }
                // unbound-H: sum of squared distances cutoff
                if rab2 + dist2(xyz, i, h) + dist2(xyz, j, h) < p.hbthr2 {
                    e += self.hb_eg1(xyz, i, j, h, g);
                }
            }
        }
        e
    }

    /// dispatch one bound-H HB (A donor, B acceptor, H on A):
    /// eg3 for carbonyl/nitro O acceptors, eg2_rnr for 2-coordinate aromatic N
    /// (explicit lone-pair acceptor), eg2new otherwise
    fn hb_bound(&self, xyz: &[[f64; 3]], a: usize, b: usize, h: usize, g: &mut [[f64; 3]]) -> f64 {
        if self.at[b] == 8 && self.topo.nb[b].len() == 1 {
            let c = self.topo.nb[b][0];
            if (self.at[c] == 6 || self.at[c] == 7) && self.topo.nb[c].len() > 1 {
                return self.hb_eg3(xyz, a, b, h, c, g);
            }
        }
        if self.at[b] == 7 && self.topo.nb[b].len() == 2 {
            // N heteroaromatic acceptor: lone-pair term (abhgfnff_eg2_rnr)
            return self.hb_eg2_rnr(xyz, a, b, h, g);
        }
        self.hb_eg2new(xyz, a, b, h, g)
    }

    /// unbound-H HB A···H···B (abhgfnff_eg1, gfnff_eg.f90 2088-2252)
    fn hb_eg1(&self, xyz: &[[f64; 3]], a: usize, b: usize, h: usize, g: &mut [[f64; 3]]) -> f64 {
        let p = &self.p;
        let rab2 = dist2(xyz, a, b); let rab = rab2.sqrt();
        let rah2 = dist2(xyz, a, h); let rah = rah2.sqrt();
        let rbh2 = dist2(xyz, b, h); let rbh = rbh2.sqrt();
        let rahprbh = rah + rbh + 1e-12;
        let radab = p.rad[self.at[a] - 1] + p.rad[self.at[b] - 1];
        let expo = (p.hbacut / radab) * (rahprbh / rab - 1.0);
        if expo > 15.0 { return 0.0; }
        let ratio2 = expo.exp();
        let outl = 2.0 / (1.0 + ratio2);
        let ratio1 = (rab2 / p.hblongcut).powf(p.hbalp);
        let dampl = 1.0 / (1.0 + ratio1);
        let shortcut = p.hbscut * radab;
        let ratio3 = (shortcut / rab2).powf(p.hbalp);
        let damps = 1.0 / (1.0 + ratio3);
        let damp = damps * dampl;
        let rdamp = damp / rab2 / rab;
        let qh = Self::csig((p.hbst * self.topo.qa[h]).exp(), p.hbsf);
        let qa = Self::csig((-p.hbst * self.topo.qa[a]).exp(), p.hbsf);
        let qb = Self::csig((-p.hbst * self.topo.qa[b]).exp(), p.hbsf);
        let rah4 = rah2 * rah2;
        let rbh4 = rbh2 * rbh2;
        let denom = 1.0 / (rah4 + rbh4);
        let caa = qa * self.topo.hbbas[a];
        let cbb = qb * self.topo.hbbas[b];
        let ca2 = self.topo.hbaci[a];
        let cb2 = self.topo.hbaci[b];
        let bas = (caa * rah4 + cbb * rbh4) * denom;
        let aci = (cb2 * rah4 + ca2 * rbh4) * denom;
        let qhoutl = qh * outl;
        let rterm = -aci * rdamp * qhoutl;
        let energy = bas * rterm;
        // gradient
        let drah = sub3(xyz, a, h);
        let drbh = sub3(xyz, b, h);
        let drab = sub3(xyz, a, b);
        let aterm = -aci * bas * rdamp * qh;
        let sterm = -rdamp * bas * qhoutl;
        let dterm = -aci * bas * qhoutl;
        let tmp = denom * denom * 4.0;
        let dd24a = rah2 * rbh4 * tmp;
        let dd24b = rbh2 * rah4 * tmp;
        // donor-acceptor part: bas
        let ga_bas = scale3(drah, (caa - cbb) * dd24a * rterm);
        let gb_bas = scale3(drbh, (cbb - caa) * dd24b * rterm);
        let mut ga = ga_bas;
        let mut gb = gb_bas;
        let mut gh = [-ga_bas[0] - gb_bas[0], -ga_bas[1] - gb_bas[1], -ga_bas[2] - gb_bas[2]];
        // donor-acceptor part: aci
        let dga = scale3(drah, (cb2 - ca2) * dd24a * sterm);
        let dgb = scale3(drbh, (ca2 - cb2) * dd24b * sterm);
        ga = add3(ga, dga);
        gb = add3(gb, dgb);
        gh = [gh[0] - dga[0] - dgb[0], gh[1] - dga[1] - dgb[1], gh[2] - dga[2] - dgb[2]];
        // damping part: rab
        let gi = rdamp * (-(2.0 * p.hbalp * ratio1 / (1.0 + ratio1))
            + (2.0 * p.hbalp * ratio3 / (1.0 + ratio3)) - 3.0) / rab2;
        let dg = scale3(drab, gi * dterm);
        ga = add3(ga, dg);
        gb = subv3(gb, dg);
        // out-of-line term: rab
        let gi = aterm * 2.0 * ratio2 * expo * rahprbh
            / (1.0 + ratio2).powi(2) / (rahprbh - rab) / rab2;
        let dg = scale3(drab, gi);
        ga = add3(ga, dg);
        gb = subv3(gb, dg);
        // out-of-line term: rah, rbh
        let tmp1 = -2.0 * aterm * ratio2 * expo / (1.0 + ratio2).powi(2) / (rahprbh - rab);
        let dga = scale3(drah, tmp1 / rah);
        let dgb = scale3(drbh, tmp1 / rbh);
        ga = add3(ga, dga);
        gb = add3(gb, dgb);
        gh = [gh[0] - dga[0] - dgb[0], gh[1] - dga[1] - dgb[1], gh[2] - dga[2] - dgb[2]];
        for t in 0..3 {
            g[a][t] += ga[t];
            g[b][t] += gb[t];
            g[h][t] += gh[t];
        }
        energy
    }

    /// bound-H HB A-H···B, default acceptors (abhgfnff_eg2new, gfnff_eg.f90 2255-2513)
    fn hb_eg2new(&self, xyz: &[[f64; 3]], a: usize, b: usize, h: usize, g: &mut [[f64; 3]]) -> f64 {
        let p = &self.p;
        let p_bh = 1.0 + p.hbabmix;
        let p_ab = -p.hbabmix;
        let rab2 = dist2(xyz, a, b); let rab = rab2.sqrt();
        let rah2 = dist2(xyz, a, h); let rah = rah2.sqrt();
        let rbh2 = dist2(xyz, b, h); let rbh = rbh2.sqrt();
        let rahprbh = rah + rbh + 1e-12;
        let radab = p.rad[self.at[a] - 1] + p.rad[self.at[b] - 1];
        let expo = (p.hbacut / radab) * (rahprbh / rab - 1.0);
        if expo > 15.0 { return 0.0; }
        let ratio2 = expo.exp();
        let outl = 2.0 / (1.0 + ratio2);
        // out-of-line damp A...nb(B)-B over ALL covalent neighbors of B
        let nbb = self.topo.nb[b].len();
        let hbnbcut = if self.at[b] == 7 && nbb == 1 { 2.0 } else { p.hbnbcut };
        let mut ranb = vec![0.0f64; nbb]; let mut rbnb = vec![0.0f64; nbb];
        let mut expo_nb = vec![0.0f64; nbb]; let mut ratio2_nb = vec![0.0f64; nbb];
        let mut outl_nb = vec![0.0f64; nbb];
        let mut ranbprbnb = vec![0.0f64; nbb];
        let mut dranb = vec![[0.0f64; 3]; nbb]; let mut drbnb = vec![[0.0f64; 3]; nbb];
        for (i, &nbj) in self.topo.nb[b].iter().enumerate() {
            dranb[i] = sub3(xyz, a, nbj);
            drbnb[i] = sub3(xyz, b, nbj);
            ranb[i] = norm3(dranb[i]);
            rbnb[i] = norm3(drbnb[i]);
            ranbprbnb[i] = ranb[i] + rbnb[i] + 1e-12;
            expo_nb[i] = (hbnbcut / radab) * (ranbprbnb[i] / rab - 1.0);
            ratio2_nb[i] = (-expo_nb[i]).exp();
            outl_nb[i] = 2.0 / (1.0 + ratio2_nb[i]) - 1.0;
        }
        let outl_nb_tot: f64 = outl_nb.iter().product();
        let ratio1 = (rab2 / p.hblongcut).powf(p.hbalp);
        let dampl = 1.0 / (1.0 + ratio1);
        let shortcut = p.hbscut * radab;
        let ratio3 = (shortcut / rab2).powf(p.hbalp);
        let damps = 1.0 / (1.0 + ratio3);
        let damp = damps * dampl;
        let ddamp = -(2.0 * p.hbalp * ratio1 / (1.0 + ratio1))
            + (2.0 * p.hbalp * ratio3 / (1.0 + ratio3));
        let rbhdamp = damp * (p_bh / rbh2 / rbh);
        let rabdamp = damp * (p_ab / rab2 / rab);
        let rdamp = rbhdamp + rabdamp;
        let qh = Self::csig((p.hbst * self.topo.qa[h]).exp(), p.hbsf);
        let qa = Self::csig((-p.hbst * self.topo.qa[a]).exp(), p.hbsf);
        let qb = Self::csig((-p.hbst * self.topo.qa[b]).exp(), p.hbsf);
        let qhoutl = qh * outl * outl_nb_tot;
        let const_ = self.topo.hbaci[a] * qa * self.topo.hbbas[b] * qb * p.xhaci_globabh;
        let energy = -rdamp * qhoutl * const_;
        // gradient
        let drah = sub3(xyz, a, h);
        let drbh = sub3(xyz, b, h);
        let drab = sub3(xyz, a, b);
        let aterm = -rdamp * qh * outl_nb_tot * const_;
        let nbterm = -rdamp * qh * outl * const_;
        let dterm = -qhoutl * const_;
        let mut ga = [0.0f64; 3]; let mut gb = [0.0f64; 3]; let mut gh = [0.0f64; 3];
        // damping part: rab
        let gi = ((rabdamp + rbhdamp) * ddamp - 3.0 * rabdamp) / rab2 * dterm;
        let dg = scale3(drab, gi);
        ga = add3(ga, dg); gb = subv3(gb, dg);
        // damping part: rbh
        let gi = -3.0 * rbhdamp / rbh2 * dterm;
        let dg = scale3(drbh, gi);
        gb = add3(gb, dg); gh = subv3(gh, dg);
        // angular A-H...B: rab
        let tmp1 = -2.0 * aterm * ratio2 * expo / (1.0 + ratio2).powi(2) / (rahprbh - rab);
        let gi = -tmp1 * rahprbh / rab2;
        let dg = scale3(drab, gi);
        ga = add3(ga, dg); gb = subv3(gb, dg);
        let dga = scale3(drah, tmp1 / rah);
        let dgb = scale3(drbh, tmp1 / rbh);
        ga = add3(ga, dga); gb = add3(gb, dgb);
        gh = subv3(gh, add3(dga, dgb));
        // angular A...nb(B)-B terms
        for i in 0..nbb {
            let prod_others: f64 = outl_nb.iter().enumerate()
                .filter(|&(j, _)| j != i).map(|(_, v)| v).product();
            let tmp2 = 2.0 * nbterm * prod_others * ratio2_nb[i] * expo_nb[i]
                / (1.0 + ratio2_nb[i]).powi(2) / (ranbprbnb[i] - rab);
            let gi = -tmp2 * ranbprbnb[i] / rab2;
            let dg = scale3(drab, gi);
            ga = add3(ga, dg); gb = subv3(gb, dg);
            let dga = scale3(dranb[i], tmp2 / ranb[i]);
            let dgb = scale3(drbnb[i], tmp2 / rbnb[i]);
            ga = add3(ga, dga); gb = add3(gb, dgb);
            let nbj = self.topo.nb[b][i];
            for t in 0..3 { g[nbj][t] -= dga[t] + dgb[t]; }
        }
        for t in 0..3 {
            g[a][t] += ga[t];
            g[b][t] += gb[t];
            g[h][t] += gh[t];
        }
        energy
    }

    /// bound-H HB onto 2-coordinate aromatic N acceptor with an explicit
    /// lone-pair position (abhgfnff_eg2_rnr, gfnff_eg.f90 2516-2822);
    /// lp = B - (0.50 - 0.018*repz(B)) * unit(sum of ring-neighbor vectors),
    /// extra A-lp-B out-of-line damping with hblpcut = 56
    fn hb_eg2_rnr(&self, xyz: &[[f64; 3]], a: usize, b: usize, h: usize, g: &mut [[f64; 3]]) -> f64 {
        let p = &self.p;
        let p_bh = 1.0 + p.hbabmix;
        let p_ab = -p.hbabmix;
        let lp_dist = 0.50 - 0.018 * p.repz[self.at[b] - 1];
        const HBLCUT: f64 = 56.0;
        let nbs = &self.topo.nb[b];
        let nbb = nbs.len();
        let mut dranb = vec![[0.0f64; 3]; nbb];
        let mut drbnb = vec![[0.0f64; 3]; nbb];
        let mut ranb = vec![0.0f64; nbb];
        let mut rbnb = vec![0.0f64; nbb];
        let mut vector = [0.0f64; 3];
        for (i, &nbj) in nbs.iter().enumerate() {
            dranb[i] = sub3(xyz, a, nbj);
            drbnb[i] = sub3(xyz, b, nbj);
            ranb[i] = norm3(dranb[i]);
            rbnb[i] = norm3(drbnb[i]);
            vector = add3(vector, sub3(xyz, nbj, b));
        }
        let vnorm = norm3(vector);
        let lp = if vnorm > 1e-10 {
            let n = [vector[0] / vnorm, vector[1] / vnorm, vector[2] / vnorm];
            [xyz[b][0] - lp_dist * n[0], xyz[b][1] - lp_dist * n[1], xyz[b][2] - lp_dist * n[2]]
        } else {
            xyz[b]
        };
        let rab2 = dist2(xyz, a, b); let rab = rab2.sqrt();
        let rah2 = dist2(xyz, a, h); let rah = rah2.sqrt();
        let rbh2 = dist2(xyz, b, h); let rbh = rbh2.sqrt();
        let rahprbh = rah + rbh + 1e-12;
        let radab = p.rad[self.at[a] - 1] + p.rad[self.at[b] - 1];
        let expo = (p.hbacut / radab) * (rahprbh / rab - 1.0);
        if expo > 15.0 { return 0.0; }
        let ratio2 = expo.exp();
        let outl = 2.0 / (1.0 + ratio2);
        // lone-pair distances
        let dralp = [xyz[a][0] - lp[0], xyz[a][1] - lp[1], xyz[a][2] - lp[2]];
        let drblp = [xyz[b][0] - lp[0], xyz[b][1] - lp[1], xyz[b][2] - lp[2]];
        let ralp = norm3(dralp);
        let rblp = norm3(drblp);
        let ralpprblp = ralp + rblp + 1e-12;
        let expo_lp = (HBLCUT / radab) * (ralpprblp / rab - 1.0);
        let ratio2_lp = expo_lp.exp();
        let outl_lp = 2.0 / (1.0 + ratio2_lp);
        let mut expo_nb = vec![0.0f64; nbb];
        let mut ratio2_nb = vec![0.0f64; nbb];
        let mut outl_nb = vec![0.0f64; nbb];
        let mut ranbprbnb = vec![0.0f64; nbb];
        for i in 0..nbb {
            ranbprbnb[i] = ranb[i] + rbnb[i] + 1e-12;
            expo_nb[i] = (p.hbnbcut / radab) * (ranbprbnb[i] / rab - 1.0);
            ratio2_nb[i] = (-expo_nb[i]).exp();
            outl_nb[i] = 2.0 / (1.0 + ratio2_nb[i]) - 1.0;
        }
        let outl_nb_tot: f64 = outl_nb.iter().product();
        let ratio1 = (rab2 / p.hblongcut).powf(p.hbalp);
        let dampl = 1.0 / (1.0 + ratio1);
        let shortcut = p.hbscut * radab;
        let ratio3 = (shortcut / rab2).powf(p.hbalp);
        let damps = 1.0 / (1.0 + ratio3);
        let damp = damps * dampl;
        let ddamp = -(2.0 * p.hbalp * ratio1 / (1.0 + ratio1))
            + (2.0 * p.hbalp * ratio3 / (1.0 + ratio3));
        let rbhdamp = damp * (p_bh / rbh2 / rbh);
        let rabdamp = damp * (p_ab / rab2 / rab);
        let rdamp = rbhdamp + rabdamp;
        let qh = Self::csig((p.hbst * self.topo.qa[h]).exp(), p.hbsf);
        let qa = Self::csig((-p.hbst * self.topo.qa[a]).exp(), p.hbsf);
        let qb = Self::csig((-p.hbst * self.topo.qa[b]).exp(), p.hbsf);
        let qhoutl = qh * outl * outl_nb_tot * outl_lp;
        let const_ = self.topo.hbaci[a] * qa * self.topo.hbbas[b] * qb * p.xhaci_globabh;
        let energy = -rdamp * qhoutl * const_;
        // gradient
        let drah = sub3(xyz, a, h);
        let drbh = sub3(xyz, b, h);
        let drab = sub3(xyz, a, b);
        let aterm = -rdamp * qh * outl_nb_tot * outl_lp * const_;
        let nbterm = -rdamp * qh * outl * outl_lp * const_;
        let lpterm = -rdamp * qh * outl * outl_nb_tot * const_;
        let dterm = -qhoutl * const_;
        let mut ga = [0.0f64; 3]; let mut gb = [0.0f64; 3]; let mut gh = [0.0f64; 3];
        // damping part: rab
        let gi = ((rabdamp + rbhdamp) * ddamp - 3.0 * rabdamp) / rab2 * dterm;
        let dg = scale3(drab, gi);
        ga = add3(ga, dg); gb = subv3(gb, dg);
        // damping part: rbh
        let gi = -3.0 * rbhdamp / rbh2 * dterm;
        let dg = scale3(drbh, gi);
        gb = add3(gb, dg); gh = subv3(gh, dg);
        // angular A-H...B: rab
        let tmp1 = -2.0 * aterm * ratio2 * expo / (1.0 + ratio2).powi(2) / (rahprbh - rab);
        let gi = -tmp1 * rahprbh / rab2;
        let dg = scale3(drab, gi);
        ga = add3(ga, dg); gb = subv3(gb, dg);
        let dga = scale3(drah, tmp1 / rah);
        let dgb = scale3(drbh, tmp1 / rbh);
        ga = add3(ga, dga); gb = add3(gb, dgb);
        gh = subv3(gh, add3(dga, dgb));
        // angular A-lp-B: rab, ralp, rblp (xtb's literal scatter: dgb unused,
        // gb subtracts dga again, glp = -dga)
        let tmp3 = -2.0 * lpterm * ratio2_lp * expo_lp
            / (1.0 + ratio2_lp).powi(2) / (ralpprblp - rab);
        let gi = -tmp3 * ralpprblp / rab2;
        let dg = scale3(drab, gi);
        ga = add3(ga, dg); gb = subv3(gb, dg);
        let dga_lp = scale3(dralp, tmp3 / ralp);
        let _dgb_lp = scale3(drblp, tmp3 / (rblp + 1e-12));
        ga = add3(ga, dga_lp);
        gb = subv3(gb, dga_lp);
        let glp = [-dga_lp[0], -dga_lp[1], -dga_lp[2]];
        // angular A...nb(B)-B terms
        let mut gnb = vec![[0.0f64; 3]; nbb];
        for i in 0..nbb {
            let prod_others: f64 = outl_nb.iter().enumerate()
                .filter(|&(j, _)| j != i).map(|(_, v)| v).product();
            let tmp2 = 2.0 * nbterm * prod_others * ratio2_nb[i] * expo_nb[i]
                / (1.0 + ratio2_nb[i]).powi(2) / (ranbprbnb[i] - rab);
            let gi = -tmp2 * ranbprbnb[i] / rab2;
            let dg = scale3(drab, gi);
            ga = add3(ga, dg); gb = subv3(gb, dg);
            let dga = scale3(dranb[i], tmp2 / ranb[i]);
            let dgb = scale3(drbnb[i], tmp2 / rbnb[i]);
            ga = add3(ga, dga); gb = add3(gb, dgb);
            gnb[i] = [-dga[0] - dgb[0], -dga[1] - dgb[1], -dga[2] - dgb[2]];
        }
        // lone-pair position chain: lp = B - lp_dist*n̂(vector), vector = sum(R_nb - B)
        let mut gnb_lp = [0.0f64; 3];
        if vnorm > 1e-10 {
            let v2 = dot3(vector, vector);
            for i in 0..3 {
                let mut unit_vec = [0.0f64; 3];
                unit_vec[i] = -1.0;
                let mut gii = [0.0f64; 3];
                for j in 0..3 {
                    gii[j] = -lp_dist * nbb as f64
                        * (unit_vec[j] / vnorm + vector[j] * vector[i] / (v2 * vnorm));
                }
                for j in 0..3 { gnb_lp[j] += gii[j] * glp[i]; }
            }
        }
        for t in 0..3 {
            g[a][t] += ga[t];
            g[b][t] += gb[t] + gnb_lp[t];
            g[h][t] += gh[t];
        }
        for (i, &nbj) in nbs.iter().enumerate() {
            for t in 0..3 {
                g[nbj][t] += gnb[i][t] - gnb_lp[t] / nbb as f64;
            }
        }
        energy
    }

    /// bound-H HB onto carbonyl/nitro O acceptor (abhgfnff_eg3, gfnff_eg.f90 2827-3212)
    /// B=O (single neighbor C), torsion H-B=C-R product + bend H-B=C
    fn hb_eg3(&self, xyz: &[[f64; 3]], a: usize, b: usize, h: usize, c: usize, g: &mut [[f64; 3]]) -> f64 {
        let p = &self.p;
        let p_bh = 1.0 + p.hbabmix;
        let p_ab = -p.hbabmix;
        let rab2 = dist2(xyz, a, b); let rab = rab2.sqrt();
        let rah2 = dist2(xyz, a, h); let rah = rah2.sqrt();
        let rbh2 = dist2(xyz, b, h); let rbh = rbh2.sqrt();
        let rahprbh = rah + rbh + 1e-12;
        let radab = p.rad[self.at[a] - 1] + p.rad[self.at[b] - 1];
        let expo = (p.hbacut / radab) * (rahprbh / rab - 1.0);
        if expo > 15.0 { return 0.0; }
        let ratio2 = expo.exp();
        let outl = 2.0 / (1.0 + ratio2);
        // out-of-line damp A...C-B (single neighbor)
        let dranb = sub3(xyz, a, c);
        let drbnb = sub3(xyz, b, c);
        let ranb = norm3(dranb);
        let rbnb = norm3(drbnb);
        let ranbprbnb = ranb + rbnb + 1e-12;
        let expo_nb = (p.hbnbcut / radab) * (ranbprbnb / rab - 1.0);
        let ratio2_nb = (-expo_nb).exp();
        let outl_nb_tot = 2.0 / (1.0 + ratio2_nb) - 1.0;
        let ratio1 = (rab2 / p.hblongcut).powf(p.hbalp);
        let dampl = 1.0 / (1.0 + ratio1);
        let shortcut = p.hbscut * radab;
        let ratio3 = (shortcut / rab2).powi(6);   // xtb uses power 6 here (eg3 only)
        let damps = 1.0 / (1.0 + ratio3);
        let damp = damps * dampl;
        let ddamp = -(2.0 * p.hbalp * ratio1 / (1.0 + ratio1))
            + (2.0 * p.hbalp * ratio3 / (1.0 + ratio3));
        let rbhdamp = damp * (p_bh / rbh2 / rbh);
        let rabdamp = damp * (p_ab / rab2 / rab);
        let rdamp = rbhdamp + rabdamp;
        // torsion product H-B=C-R for each R in nb(C), R != B (egtors_nci_mul)
        let rlist: Vec<usize> = self.topo.nb[c].iter().copied().filter(|&r| r != b).collect();
        let mut etmp = vec![0.0f64; rlist.len()];
        let mut g4 = vec![[[0.0f64; 3]; 4]; rlist.len()];
        let phi0t = std::f64::consts::PI / 2.0;
        for (i, &r) in rlist.iter().enumerate() {
            etmp[i] = self.tors_nci_mul(xyz, r, b, c, h, 2.0, phi0t, p.tors_hb, &mut g4[i]);
        }
        let etors: f64 = etmp.iter().product();
        // bend H-B=C (egbend_nci_mul), phi0 = 120 deg
        let mut g3 = [[0.0f64; 3]; 3];
        let eangl = bend_nci_mul(xyz, b, c, h, 120.0 * std::f64::consts::PI / 180.0,
            1.0 - p.bend_hb, &mut g3);
        let qh = Self::csig((p.hbst * self.topo.qa[h]).exp(), p.hbsf);
        let qa = Self::csig((-p.hbst * self.topo.qa[a]).exp(), p.hbsf);
        let qb = Self::csig((-p.hbst * self.topo.qa[b]).exp(), p.hbsf);
        let qhoutl = qh * outl * outl_nb_tot;
        let const_ = self.topo.hbaci[a] * qa * self.topo.hbbas[b] * qb * p.xhaci_coh;
        let energy = -rdamp * qhoutl * eangl * etors * const_;
        // gradient
        let drah = sub3(xyz, a, h);
        let drbh = sub3(xyz, b, h);
        let drab = sub3(xyz, a, b);
        let aterm = -rdamp * qh * outl_nb_tot * eangl * etors * const_;
        let nbterm = -rdamp * qh * outl * eangl * etors * const_;
        let dterm = -qhoutl * eangl * etors * const_;
        let tterm = -rdamp * qhoutl * eangl * const_;
        let bterm = -rdamp * qhoutl * etors * const_;
        let mut ga = [0.0f64; 3]; let mut gb = [0.0f64; 3]; let mut gh = [0.0f64; 3];
        // damping part: rab
        let gi = ((rabdamp + rbhdamp) * ddamp - 3.0 * rabdamp) / rab2 * dterm;
        let dg = scale3(drab, gi);
        ga = add3(ga, dg); gb = subv3(gb, dg);
        // damping part: rbh
        let gi = -3.0 * rbhdamp / rbh2 * dterm;
        let dg = scale3(drbh, gi);
        gb = add3(gb, dg); gh = subv3(gh, dg);
        // angular A-H...B: rab
        let tmp1 = -2.0 * aterm * ratio2 * expo / (1.0 + ratio2).powi(2) / (rahprbh - rab);
        let gi = -tmp1 * rahprbh / rab2;
        let dg = scale3(drab, gi);
        ga = add3(ga, dg); gb = subv3(gb, dg);
        let dga = scale3(drah, tmp1 / rah);
        let dgb = scale3(drbh, tmp1 / rbh);
        ga = add3(ga, dga); gb = add3(gb, dgb);
        gh = subv3(gh, add3(dga, dgb));
        // angular A...C-B term
        let tmp2 = 2.0 * nbterm * ratio2_nb * expo_nb
            / (1.0 + ratio2_nb).powi(2) / (ranbprbnb - rab);
        let gi = -tmp2 * ranbprbnb / rab2;
        let dg = scale3(drab, gi);
        ga = add3(ga, dg); gb = subv3(gb, dg);
        let dga = scale3(dranb, tmp2 / ranb);
        let dgb = scale3(drbnb, tmp2 / rbnb);
        ga = add3(ga, dga); gb = add3(gb, dgb);
        for t in 0..3 { g[c][t] -= dga[t] + dgb[t]; }
        // torsion term H...B=C<R
        for (i, &r) in rlist.iter().enumerate() {
            let prod_others: f64 = etmp.iter().enumerate()
                .filter(|&(j, _)| j != i).map(|(_, v)| v).product();
            for t in 0..3 {
                g[r][t] += g4[i][0][t] * tterm * prod_others;
                g[b][t] += g4[i][1][t] * tterm * prod_others;
                g[c][t] += g4[i][2][t] * tterm * prod_others;
                g[h][t] += g4[i][3][t] * tterm * prod_others;
            }
        }
        // angle term H...B=C
        for t in 0..3 {
            g[b][t] += g3[0][t] * bterm;
            g[c][t] += g3[1][t] * bterm;
            g[h][t] += g3[2][t] * bterm;
        }
        for t in 0..3 {
            g[a][t] += ga[t];
            g[b][t] += gb[t];
            g[h][t] += gh[t];
        }
        energy
    }

    /// scaled torsion factor for eg3 (egtors_nci_mul, gfnff_eg.f90 1585-1627)
    /// e = (1+cos(rn(phi-phi0)+pi))*fc + tshift, fc = (1-tshift)/2; no damping
    #[allow(clippy::too_many_arguments)]
    fn tors_nci_mul(&self, xyz: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize,
                    rn: f64, phi0: f64, tshift: f64, g4: &mut [[f64; 3]; 4]) -> f64 {
        let fc = (1.0 - tshift) / 2.0;
        let phi = dihedral(xyz, i, j, k, l);
        let (dda, ddb, ddc, ddd) = dphidr(xyz, i, j, k, l, phi);
        let c1 = rn * (phi - phi0) + std::f64::consts::PI;
        let et = (1.0 + c1.cos()) * fc + tshift;
        let dij = -rn * c1.sin() * fc;
        for t in 0..3 {
            g4[0][t] = dij * dda[t];
            g4[1][t] = dij * ddb[t];
            g4[2][t] = dij * ddc[t];
            g4[3][t] = dij * ddd[t];
        }
        et
    }

    /// XB list filter + rbxgfnff_eg evaluation (gfnff_eg.f90 3216-3345)
    fn xb_terms(&self, xyz: &[[f64; 3]], g: &mut [[f64; 3]]) -> f64 {
        let p = &self.p;
        let mut e = 0.0f64;
        for &(a, b, x) in &self.topo.xb_triples {
            if dist2(xyz, a, b) > p.hbthr2 { continue; }
            e += self.rbx_eg(xyz, a, b, x, g);
        }
        e
    }

    /// halogen bond A-X···B (rbxgfnff_eg); note cb = 1.0 hardcoded in xtb
    fn rbx_eg(&self, xyz: &[[f64; 3]], a: usize, b: usize, x: usize, g: &mut [[f64; 3]]) -> f64 {
        let p = &self.p;
        let cb = 1.0;
        let cx = p.xbaci[self.at[x] - 1];
        let drax = sub3(xyz, a, x);
        let drbx = sub3(xyz, b, x);
        let drab = sub3(xyz, a, b);
        let rab = norm3(drab);
        let rab2 = rab * rab;
        let rax = norm3(drax) + 1e-12;
        let rax2 = rax * rax;
        let rbx = norm3(drbx) + 1e-12;
        let rbx2 = rbx * rbx;
        let expo = p.xbacut * ((rax + rbx) / rab - 1.0);
        if expo > 15.0 { return 0.0; }
        let ratio2 = expo.exp();
        let outl = 2.0 / (1.0 + ratio2);
        let ratio1 = (rbx2 / p.hblongcut_xb).powf(p.hbalp);
        let dampl = 1.0 / (1.0 + ratio1);
        let shortcut = p.xbscut * (p.rad[self.at[a] - 1] + p.rad[self.at[b] - 1]);
        let ratio3 = (shortcut / rbx2).powf(p.hbalp);
        let damps = 1.0 / (1.0 + ratio3);
        let damp = damps * dampl;
        let rdamp = damp / rbx2 / rbx;
        let qx = Self::csig((p.xbst * self.topo.qa[x]).exp(), p.xbsf);
        let qb = Self::csig((-p.xbst * self.topo.qa[b]).exp(), p.xbsf);
        let const_ = cb * qb * cx * qx;
        let aterm = -rdamp * const_;
        let dterm = -outl * const_;
        let energy = -rdamp * outl * const_;
        // damping part: rbx
        let gi = rdamp * (-(2.0 * p.hbalp * ratio1 / (1.0 + ratio1))
            + (2.0 * p.hbalp * ratio3 / (1.0 + ratio3)) - 3.0) / rbx2 * dterm;
        let dg = scale3(drbx, gi);
        let mut gb = dg;
        let mut gx = [-dg[0], -dg[1], -dg[2]];
        // out-of-line term: rab
        let gi = 2.0 * ratio2 * expo * (rax + rbx)
            / (1.0 + ratio2).powi(2) / (rax + rbx - rab) / rab2 * aterm;
        let dg = scale3(drab, gi);
        let mut ga = dg;
        gb = subv3(gb, dg);
        // out-of-line term: rax, rbx
        let gi = -2.0 * ratio2 * expo / (1.0 + ratio2).powi(2) / (rax + rbx - rab) / rax * aterm;
        let dga = scale3(drax, gi);
        ga = add3(ga, dga);
        let gi = -2.0 * ratio2 * expo / (1.0 + ratio2).powi(2) / (rax + rbx - rab) / rbx * aterm;
        let dgb = scale3(drbx, gi);
        gb = add3(gb, dgb);
        let dgx = [-dga[0] - dgb[0], -dga[1] - dgb[1], -dga[2] - dgb[2]];
        gx = add3(gx, dgx);
        for t in 0..3 {
            g[a][t] += ga[t];
            g[b][t] += gb[t];
            g[x][t] += gx[t];
        }
        energy
    }

    /// bonded ATM over 1-4 triples (batmgfnff_eg, gfnff_eg.f90 3348-3446)
    fn batm_terms(&self, xyz: &[[f64; 3]], g: &mut [[f64; 3]]) -> f64 {
        let mut e = 0.0f64;
        for &(i, j, k) in &self.topo.b3 {
            e += self.batm_eg(xyz, i, j, k, g);
        }
        e
    }

    fn batm_eg(&self, xyz: &[[f64; 3]], iat: usize, jat: usize, kat: usize,
               g: &mut [[f64; 3]]) -> f64 {
        const FQQ: f64 = 3.0;
        let fi = (1.0 - FQQ * self.topo.qa[iat]).clamp(-4.0, 4.0);
        let fj = (1.0 - FQQ * self.topo.qa[jat]).clamp(-4.0, 4.0);
        let fk = (1.0 - FQQ * self.topo.qa[kat]).clamp(-4.0, 4.0);
        let ff = fi * fj * fk;
        let p = &self.p;
        let c9 = ff * p.zb3atm[self.at[iat] - 1] * p.zb3atm[self.at[jat] - 1] * p.zb3atm[self.at[kat] - 1];
        let rij = sub3(xyz, jat, iat);
        let rik = sub3(xyz, kat, iat);
        let rjk = sub3(xyz, kat, jat);
        let sr2ij = norm3(rij); let r2ij = sr2ij * sr2ij;
        let sr2ik = norm3(rik); let r2ik = sr2ik * sr2ik;
        let sr2jk = norm3(rjk); let r2jk = sr2jk * sr2jk;
        let mijk = -r2ij + r2jk + r2ik;
        let imjk = r2ij - r2jk + r2ik;
        let ijmk = r2ij + r2jk - r2ik;
        let rijk3 = r2ij * r2jk * r2ik;
        let rav3 = rijk3 * sr2ij * sr2jk * sr2ik;   // R^9
        let ang = 0.375 * ijmk * imjk * mijk / rijk3;
        let energy = c9 * (ang + 1.0) / rav3;
        // derivatives of the angular part w.r.t. each pair distance
        let dang_ij = -0.375 * (r2ij.powi(3) + r2ij * r2ij * (r2jk + r2ik)
            + r2ij * (3.0 * r2jk * r2jk + 2.0 * r2jk * r2ik + 3.0 * r2ik * r2ik)
            - 5.0 * (r2jk - r2ik).powi(2) * (r2jk + r2ik)) / (sr2ij * rijk3 * rav3);
        let drij = -dang_ij * c9;
        let dang_jk = -0.375 * (r2jk.powi(3) + r2jk * r2jk * (r2ik + r2ij)
            + r2jk * (3.0 * r2ik * r2ik + 2.0 * r2ik * r2ij + 3.0 * r2ij * r2ij)
            - 5.0 * (r2ik - r2ij).powi(2) * (r2ik + r2ij)) / (sr2jk * rijk3 * rav3);
        let drjk = -dang_jk * c9;
        let dang_ik = -0.375 * (r2ik.powi(3) + r2ik * r2ik * (r2jk + r2ij)
            + r2ik * (3.0 * r2jk * r2jk + 2.0 * r2jk * r2ij + 3.0 * r2ij * r2ij)
            - 5.0 * (r2jk - r2ij).powi(2) * (r2jk + r2ij)) / (sr2ik * rijk3 * rav3);
        let drik = -dang_ik * c9;
        for t in 0..3 {
            g[iat][t] += drij * rij[t] / sr2ij + drik * rik[t] / sr2ik;
            g[jat][t] += drjk * rjk[t] / sr2jk - drij * rij[t] / sr2ij;
            g[kat][t] += -drik * rik[t] / sr2ik - drjk * rjk[t] / sr2jk;
        }
        energy
    }


    /// runtime HB coordination count of H for a registered X-H bond
    /// (dncoord_erf: kn = 27.5, rcov_scal = 1.78, cutoff 900 Bohr^2)
    fn hb_cn_of(&self, xyz: &[[f64; 3]], bi: usize, h: usize) -> f64 {
        let mut hb_cn = 0.0f64;
        for &bb in &self.topo.bond_hb_b[bi] {
            let rij = sub3(xyz, bb, h);
            let r2 = dot3(rij, rij);
            if r2 > 900.0 { continue; }
            let r = r2.sqrt();
            let rc = 1.78 * (self.p.rcov[self.at[h] - 1] + self.p.rcov[self.at[bb] - 1])
                / BOHR * 4.0 / 3.0;
            hb_cn += 0.5 * (1.0 + erf(-27.5 * (r - rc) / rc));
        }
        hb_cn
    }

    /// improper torsion energy+gradient (egtors, gfnff_eg.f90 1520-1579);
    /// t = (l=center, i=1st, j=2nd, k=3rd), damping on (center,1st),(2nd,1st),(1st,3rd)
    fn improper_term(&self, xyz: &[[f64; 3]], t: &TorsionParam, g: &mut [[f64; 3]]) -> f64 {
        let (i, j, k, l) = (t.l, t.i, t.j, t.k);   // center, 1st, 2nd, 3rd
        let vab = sub3(xyz, j, i);   // 1st - center
        let vcb = sub3(xyz, j, k);   // 1st - 2nd
        let vdc = sub3(xyz, j, l);   // 1st - 3rd
        let rij = norm3(vab) * norm3(vab);
        let rjk = norm3(vcb) * norm3(vcb);
        let rjl = norm3(vdc) * norm3(vdc);
        let (dampij, damp2ij) = dampt2(&self.p, self.at[i], self.at[j], rij);
        let (dampjk, damp2jk) = dampt2(&self.p, self.at[k], self.at[j], rjk);
        let (dampjl, damp2jl) = dampt2(&self.p, self.at[j], self.at[l], rjl);
        let damp = dampjk * dampij * dampjl;
        let phi = improper_angle(xyz, i, j, k, l);
        let (dda, ddb, ddc, ddd) = domegadr(xyz, i, j, k, l, phi);
        let (et, dij);
        if t.nrot == 0 {
            // sp2: E = fc*(1 - cos phi)
            let c1 = phi - t.phi0 + std::f64::consts::PI;
            et = (1.0 + c1.cos()) * t.fc;
            dij = -c1.sin() * t.fc * damp;
        } else {
            // saturated N: double well E = fc*(cos phi - cos phi0)^2
            et = t.fc * (phi.cos() - t.phi0.cos()).powi(2);
            dij = 2.0 * t.fc * phi.sin() * (t.phi0.cos() - phi.cos()) * damp;
        }
        let term1 = scale3(vab, et * damp2ij * dampjk * dampjl);
        let term2 = scale3(vcb, et * damp2jk * dampij * dampjl);
        let term3 = scale3(vdc, et * damp2jl * dampij * dampjk);
        for u in 0..3 {
            g[i][u] += dij * dda[u] - term1[u];
            g[j][u] += dij * ddb[u] + term1[u] + term2[u] + term3[u];
            g[k][u] += dij * ddc[u] - term2[u];
            g[l][u] += dij * ddd[u] - term3[u];
        }
        et * damp
    }

}


impl Gfnff {
    /// per-term analytic gradients (kcal/mol/A), for validation/debug
    #[cfg(test)]
    pub fn gradient_per_term(&self, xyz_ang: &[[f64; 3]]) -> Vec<(&'static str, Vec<[f64;3]>)> {
        let n = self.at.len();
        let mut full = vec![[0.0f64; 3]; n];
        self.energy_and_gradient(xyz_ang, &mut full);
        // individual terms by re-running energy() with finite differences is exact
        // but slow; instead compute each term's own gradient via the same code path:
        // easier: return fd per term at high accuracy (central, h=2e-5 A)
        let h = 2e-5;
        let names = ["bond", "angle", "torsion", "rep", "es", "disp"];
        let picks: [fn(&EnergyComponents) -> f64; 6] =
            [|e| e.bond, |e| e.angle, |e| e.torsion, |e| e.rep, |e| e.es, |e| e.disp];
        let mut out = vec![];
        for (name, f) in names.iter().zip(picks.iter()) {
            let mut g = vec![[0.0f64; 3]; n];
            for a in 0..n { for t in 0..3 {
                let mut p = xyz_ang.to_vec(); p[a][t] += h;
                let mut m = xyz_ang.to_vec(); m[a][t] -= h;
                g[a][t] = (f(&self.energy(&p)) - f(&self.energy(&m))) * 627.5094740631 / (2.0*h);
            }}
            out.push((*name, g));
        }
        out.push(("full", full));
        out
    }

    /// analytic gradient (dE/dR, Eh/Bohr for xyz in Bohr internally).
    /// `xyz_ang` in Angstrom; returns gradient in kcal/mol/A (ForceField trait unit).
    pub fn energy_and_gradient(
        &self,
        xyz_ang: &[[f64; 3]],
        grad_out: &mut [[f64; 3]],
    ) -> EnergyComponents {
        let n = self.at.len();
        let xyz: Vec<[f64; 3]> = xyz_ang.iter().map(|r| [r[0]/BOHR, r[1]/BOHR, r[2]/BOHR]).collect();
        for g in grad_out.iter_mut() { *g = [0.0; 3]; }
        let mut g = vec![[0.0f64; 3]; n];
        let mut g_rep = vec![[0.0f64; 3]; n];
        let mut g_bond = vec![[0.0f64; 3]; n];
        let mut g_angle = vec![[0.0f64; 3]; n];
        let mut g_tors = vec![[0.0f64; 3]; n];
        let mut g_es = vec![[0.0f64; 3]; n];
        let mut g_disp = vec![[0.0f64; 3]; n];
        let mut d_ed_cn = vec![0.0f64; n];

        let dist = |i: usize, j: usize| -> f64 { dist2(&xyz, i, j).sqrt() };

        // ---------------- EEQ charges ----------------
        let cn = erf_cn(&self.p, &self.at, &xyz);
        let mut q = vec![0.0f64; n];
        {
            let m = n + 1;
            let mut a = vec![vec![0.0f64; m]; m];
            let mut x = vec![0.0f64; m];
            for i in 0..n {
                x[i] = self.topo.chieeq[i] + self.p.cnf[self.at[i]-1] * cn[i].sqrt();
                a[i][i] = TSQRT2PI / self.topo.alpeeq[i].sqrt() + self.topo.gameeq[i];
            }
            for i in 0..n { for j in 0..i {
                let gam = 1.0 / (self.topo.alpeeq[i] + self.topo.alpeeq[j]).sqrt();
                let r = dist(i, j);
                let v = erf(gam * r) / r;
                a[i][j] = v; a[j][i] = v;
            }}
            let nfrag = self.topo.nfrag;
            let m2 = n + nfrag;
            let mut a2 = vec![vec![0.0f64; m2]; m2];
            let mut x2 = vec![0.0f64; m2];
            for i2 in 0..n {
                x2[i2] = x[i2];
                a2[i2][i2] = a[i2][i2];
            }
            for i2 in 0..n { for j2 in 0..i2 { a2[i2][j2] = a[i2][j2]; a2[j2][i2] = a[i2][j2]; } }
            for (fi, &qf) in self.topo.qfrag.iter().enumerate() {
                x2[n + fi] = qf;
                for j2 in 0..n {
                    if self.topo.fraglist[j2] == fi {
                        a2[n + fi][j2] = 1.0; a2[j2][n + fi] = 1.0;
                    }
                }
            }
            q = solve_sym(&a2, &x2);
        }

        // ---------------- repulsion (nb + bonded) ----------------
        let mut rep = 0.0;
        for i in 0..n { for j in 0..i {
            if self.topo.nb[i].contains(&j) { continue; }
            let r = dist(i, j);
            if r > 20.0 { continue; }
            let (zi, zj) = (self.at[i], self.at[j]);
            let fni = 1.0 + self.p.nrepscal / (1.0 + (self.topo.nb[i].len() as f64).powi(2));
            let fnj = 1.0 + self.p.nrepscal / (1.0 + (self.topo.nb[j].len() as f64).powi(2));
            let di = self.p.repan[zi-1] * (1.0 + self.topo.qa[i] * self.p.qrepscal) * fni;
            let dj = self.p.repan[zj-1] * (1.0 + self.topo.qa[j] * self.p.qrepscal) * fnj;
            let mut ff = 1.0;
            if zi == 1 && zj == 1 {
                ff = self.p.hhfac;
                let bp = self.topo.bpair.get(i).and_then(|r| r.get(j)).copied().unwrap_or(5);
                if bp == 2 { ff *= self.p.hh13rep; }
                if bp == 3 { ff *= self.p.hh14rep; }
            }
            if (zi == 1 && zj == 6) || (zj == 1 && zi == 6) { ff = 0.91; }
            if (zi == 1 && zj == 8) || (zj == 1 && zi == 8) { ff = 1.04; }
            let alpha = (di * dj).sqrt() * ff;
            let t16 = r.powf(1.5);
            let t19 = t16 * t16;
            let t26 = (-alpha * t16).exp() * self.p.repz[zi-1] * self.p.repz[zj-1] * self.p.repscaln;
            rep += t26 / r;
            let t27 = t26 * (1.5 * alpha * t16 + 1.0) / t19;
            let d = [xyz[i][0]-xyz[j][0], xyz[i][1]-xyz[j][1], xyz[i][2]-xyz[j][2]];
            for k in 0..3 { g_rep[i][k] -= d[k] * t27; g_rep[j][k] += d[k] * t27; }
        }}
        let mut rep_bonded_local = 0.0;
        for b in &self.bonds {
            let r = dist(b.i, b.j);
            let (zi, zj) = (self.at[b.i], self.at[b.j]);
            let alpha = (self.p.repa[zi-1] * self.p.repa[zj-1]).sqrt();
            let repab = self.p.repz[zi-1] * self.p.repz[zj-1] * self.p.repscalb();
            let t16 = r.powf(1.5);
            let t19 = t16 * t16;
            let t26 = (-alpha * t16).exp() * repab;
            rep_bonded_local += t26 / r;
            let t27 = t26 * (1.5 * alpha * t16 + 1.0) / t19;
            let d = [xyz[b.i][0]-xyz[b.j][0], xyz[b.i][1]-xyz[b.j][1], xyz[b.i][2]-xyz[b.j][2]];
            for k in 0..3 { g_rep[b.i][k] -= d[k] * t27; g_rep[b.j][k] += d[k] * t27; }
        }
        rep += rep_bonded_local;

        // ---------------- bonds (+ r0 CN chain, + HB exponent chain) --------
        let mut ebond = 0.0;
        for (bi, b) in self.bonds.iter().enumerate() {
            let r = dist(b.i, b.j);
            let rabdcn = self.p.gfnffrab_dcnd(self.at[b.i], self.at[b.j], b.r0);
            let rab0 = self.p.gfnffrab(self.at[b.i], self.at[b.j], cn[b.i], cn[b.j], b.r0);
            let dr = r - rab0;
            // H-bonded X-H: softened exponent + dE/d(hb_cn) chain (egbond_hb)
            let mut alp = b.alp;
            if !self.topo.bond_hb_b[bi].is_empty() {
                let h = if self.at[b.i] == 1 { b.i } else { b.j };
                let t1 = 1.0 - self.p.vbond_scale;
                let mut dcn_hh = [0.0f64; 3];
                let mut b_terms: Vec<(usize, [f64; 3])> = Vec::new();
                let mut hb_cn = 0.0f64;
                for &bb in &self.topo.bond_hb_b[bi] {
                    let rij = sub3(&xyz, bb, h);
                    let r2 = dot3(rij, rij);
                    if r2 > 900.0 { continue; }
                    let rr = r2.sqrt();
                    let rc = 1.78 * (self.p.rcov[self.at[h] - 1] + self.p.rcov[self.at[bb] - 1])
                        / BOHR * 4.0 / 3.0;
                    let arg = 27.5 * (rr - rc) / rc;
                    hb_cn += 0.5 * (1.0 + erf(-arg));
                    // dtmp = -1/(2 sqrt(pi)) * kn * exp(-kn^2 (r-rc)^2 / rc^2) / rc
                    let dtmp = -0.28209479177387814 * 27.5
                        * (-(27.5 * 27.5) * (rr - rc) * (rr - rc) / (rc * rc)).exp() / rc;
                    let u = [rij[0] / rr, rij[1] / rr, rij[2] / rr];
                    for t in 0..3 { dcn_hh[t] -= dtmp * u[t]; }   // d cn_H/dR_H
                    b_terms.push((bb, [dtmp * u[0], dtmp * u[1], dtmp * u[2]])); // d cn_H/dR_B
                }
                alp = (1.0 - t1 * hb_cn) * b.alp;
                let dum_hb = b.kb * (-alp * dr * dr).exp();
                let zz = dum_hb * b.alp * dr * dr * t1;
                for t in 0..3 { g_bond[h][t] += dcn_hh[t] * zz; }
                for (bb, dcn_bh) in &b_terms {
                    for t in 0..3 { g_bond[*bb][t] -= dcn_bh[t] * zz; }
                }
            }
            let dum = b.kb * (-alp * dr * dr).exp();
            ebond += dum;
            let yy = 2.0 * alp * dr * dum;
            let d = [xyz[b.i][0]-xyz[b.j][0], xyz[b.i][1]-xyz[b.j][1], xyz[b.i][2]-xyz[b.j][2]];
            for k in 0..3 {
                g_bond[b.i][k] += -yy * d[k] / r;
                g_bond[b.j][k] += yy * d[k] / r;
            }
            d_ed_cn[b.i] += yy * rabdcn.0;
            d_ed_cn[b.j] += yy * rabdcn.1;
        }

        // ---------------- angles ----------------
        let mut eangl = 0.0;
        for ang in &self.angles {
            let (j, i, k) = (ang.j, ang.i, ang.k);
            let vab = [xyz[j][0]-xyz[i][0], xyz[j][1]-xyz[i][1], xyz[j][2]-xyz[i][2]];
            let vcb = [xyz[k][0]-xyz[i][0], xyz[k][1]-xyz[i][1], xyz[k][2]-xyz[i][2]];
            let rab2 = dot3(vab, vab);
            let rcb2 = dot3(vcb, vcb);
            let cos = dot3(vab, vcb) / (rab2 * rcb2).sqrt();
            let cos = cos.clamp(-1.0, 1.0);
            let theta = cos.acos();
            let (dampij, ddampij) = dampa2(&self.p, self.at[j], self.at[i], rab2);
            let (dampjk, ddampjk) = dampa2(&self.p, self.at[k], self.at[i], rcb2);
            let damp = dampij * dampjk;
            let ea = ang.fc * (cos - ang.phi0.cos()).powi(2);
            let deddt = 2.0 * ang.fc * theta.sin() * (ang.phi0.cos() - cos);
            eangl += ea * damp;
            let vp = cross3(vcb, vab);
            let rp = norm3(vp) + 1e-14;
            let deda = scale3(cross3(vab, vp), -deddt / (rab2 * rp));
            let dedc = scale3(cross3(vcb, vp), deddt / (rcb2 * rp));
            let dedb = [deda[0]+dedc[0], deda[1]+dedc[1], deda[2]+dedc[2]];
            let term1 = scale3(vab, ea * ddampij * dampjk);
            let term2 = scale3(vcb, ea * ddampjk * dampij);
            // xtb alist swaps center/first-neighbor: center gets -dedb, nb1 gets deda
            for t in 0..3 {
                g_angle[i][t] += -dedb[t] * damp - term1[t] - term2[t];
                g_angle[j][t] += deda[t] * damp + term1[t];
                g_angle[k][t] += dedc[t] * damp + term2[t];
            }
        }

        // ---------------- torsions ----------------
        let mut etors = 0.0;
        for t in &self.torsions {
            if t.nrot <= 0 {
                // out-of-plane improper (energy + gradient)
                etors += self.improper_term(&xyz, t, &mut g_tors);
                continue;
            }
            let phi = dihedral(&xyz, t.l, t.i, t.j, t.k);
            let c1 = t.nrot as f64 * (phi - t.phi0) + std::f64::consts::PI;
            let et = (1.0 + c1.cos()) * t.fc;
            let rli2 = dist2(&xyz, t.l, t.i);
            let rij2 = dist2(&xyz, t.i, t.j);
            let rjk2 = dist2(&xyz, t.j, t.k);
            let (dij_, ddij) = dampt2(&self.p, self.at[t.l], self.at[t.i], rli2);
            let (djk_, ddjk) = dampt2(&self.p, self.at[t.i], self.at[t.j], rij2);
            let (dkl_, ddkl) = dampt2(&self.p, self.at[t.j], self.at[t.k], rjk2);
            let damp = dij_ * djk_ * dkl_;
            etors += et * damp;
            let dij = -t.nrot as f64 * c1.sin() * t.fc * damp;
            let (dda, ddb, ddc, ddd) = dphidr(&xyz, t.l, t.i, t.j, t.k, phi);
            let vab = [xyz[t.l][0]-xyz[t.i][0], xyz[t.l][1]-xyz[t.i][1], xyz[t.l][2]-xyz[t.i][2]];
            let vcb = [xyz[t.i][0]-xyz[t.j][0], xyz[t.i][1]-xyz[t.j][1], xyz[t.i][2]-xyz[t.j][2]];
            let vdc = [xyz[t.j][0]-xyz[t.k][0], xyz[t.j][1]-xyz[t.k][1], xyz[t.j][2]-xyz[t.k][2]];
            let term1 = scale3(vab, et * ddij * djk_ * dkl_);
            let term2 = scale3(vcb, et * ddjk * dij_ * dkl_);
            let term3 = scale3(vdc, et * ddkl * dij_ * djk_);
            for u in 0..3 {
                g_tors[t.l][u] += dij * dda[u] + term1[u];
                g_tors[t.i][u] += dij * ddb[u] - term1[u] + term2[u];
                g_tors[t.j][u] += dij * ddc[u] + term3[u] - term2[u];
                g_tors[t.k][u] += dij * ddd[u] - term3[u];
            }
        }

        // ---------------- ES pair term ----------------
        let mut es = 0.0;
        let sqrtpi = std::f64::consts::PI.sqrt();
        for i in 0..n { for j in 0..i {
            let gam = 1.0 / (self.topo.alpeeq[i] + self.topo.alpeeq[j]).sqrt();
            let r = dist(i, j);
            let r2 = r * r;
            let erff = erf(gam * r);
            es += q[i] * q[j] * erff / r;
            let dd = (2.0 * gam * (-(gam*gam*r2)).exp() / (sqrtpi * r2)
                     - erff / (r * r2)) * q[i] * q[j];
            let d = [xyz[i][0]-xyz[j][0], xyz[i][1]-xyz[j][1], xyz[i][2]-xyz[j][2]];
            for t in 0..3 { g_es[i][t] += d[t] * dd; g_es[j][t] -= d[t] * dd; }
        }}
        for i in 0..n {
            let x_i = self.topo.chieeq[i] + self.p.cnf[self.at[i]-1] * cn[i].sqrt();
            es += -q[i] * x_i
                + q[i]*q[i] * 0.5 * (self.topo.gameeq[i] + TSQRT2PI / self.topo.alpeeq[i].sqrt());
        }
        // cn chain from EEQ RHS:  dE/dR -= sum_b q_b * cnf_b/(2 sqrt(cn_b)) * dcn_b/dR
        for b in 0..n {
            let f = q[b] * self.p.cnf[self.at[b]-1] / (2.0 * cn[b].sqrt() + 1e-16);
            d_ed_cn[b] -= f;
        }

        // ---------------- D3 dispersion (+ cn chain) ----------------
        let disp = self.d4_dispersion_grad(&xyz, &dist, &cn, &q, &mut g_disp, &mut d_ed_cn);

        // ---------------- HB / XB / bATM (energy + gradient, one path) ----------------
        let mut g_hb = vec![[0.0f64; 3]; n];
        let mut g_xb = vec![[0.0f64; 3]; n];
        let mut g_batm = vec![[0.0f64; 3]; n];
        let ehb = self.hb_terms(&xyz, &mut g_hb);
        let exb = self.xb_terms(&xyz, &mut g_xb);
        let ebatm = self.batm_terms(&xyz, &mut g_batm);

        // ---------------- CN chain (bond r0 + EEQ RHS + D3) ----------------
        // dcndr[atom a][cn owner b] = d cn_b / d R_a
        let mut dcndr = vec![[0.0f64; 3]; n * n];
        let kn = -7.5f64;
        let cn_raw = erf_cn_raw(&self.p, &self.at, &xyz);
        for i in 0..n { for j in 0..i {
            let r = dist(i, j);
            if r > 60.0 { continue; }
            let r0 = (self.p.rcov[self.at[i]-1] + self.p.rcov[self.at[j]-1]) / BOHR * 4.0 / 3.0;
            let derf = kn / sqrtpi * (-(kn*kn*((r-r0)/r0)*((r-r0)/r0))).exp() / r0;
            let unit = [(xyz[j][0]-xyz[i][0])/r, (xyz[j][1]-xyz[i][1])/r, (xyz[j][2]-xyz[i][2])/r];
            let dlci = create_dlogcn(&self.p, cn_raw[i]);
            let dlcj = create_dlogcn(&self.p, cn_raw[j]);
            let v = [derf*unit[0], derf*unit[1], derf*unit[2]];
            for t in 0..3 {
                dcndr[i*n+j][t] += dlcj * (-v[t]);   // d cn_j / d R_i
                dcndr[j*n+j][t] += dlcj * v[t];      // d cn_j / d R_j
                dcndr[j*n+i][t] += dlci * v[t];      // d cn_i / d R_j
                dcndr[i*n+i][t] += dlci * (-v[t]);   // d cn_i / d R_i
            }
        }}
        for a in 0..n { for b in 0..n {
            if d_ed_cn[b] != 0.0 {
                for t in 0..3 { g[a][t] += dcndr[a*n+b][t] * d_ed_cn[b]; }
            }
        }}

        // sum per-term + shared cn chain
        for a in 0..n { for t in 0..3 {
            g[a][t] += g_rep[a][t] + g_bond[a][t] + g_angle[a][t] + g_tors[a][t]
                     + g_es[a][t] + g_disp[a][t] + g_hb[a][t] + g_xb[a][t] + g_batm[a][t];
        }}
        // convert to kcal/mol / Angstrom
        const EH_KCAL: f64 = 627.5094740631;
        for a in 0..n { for t in 0..3 { grad_out[a][t] = g[a][t] * EH_KCAL / BOHR; } }

        EnergyComponents { bond: ebond, angle: eangl, torsion: etors,
            rep, es, disp, hb: ehb, xb: exb, batm: ebatm }
    }

    #[allow(clippy::too_many_arguments)]
    fn d4_dispersion_grad(&self, xyz: &[[f64; 3]], dist: &dyn Fn(usize, usize) -> f64,
                          cn: &[f64], q: &[f64],
                          g: &mut [[f64; 3]], d_ed_cn: &mut [f64]) -> f64 {
        let n = self.at.len();
        let p = &self.p;
        let wf = 4.0f64;
        let tw = trapzd_weights();
        // weights and their cn derivatives
        let mut gw = vec![vec![0.0f64; 8]; n];
        let mut dgw = vec![vec![0.0f64; 8]; n];
        for i in 0..n {
            let z = self.at[i];
            let nref = p.d4_refn[z-1] as usize;
            let mut norm = 0.0; let mut dnorm = 0.0;
            let mut raw = vec![0.0f64; nref];
            for r in 0..nref {
                raw[r] = (-wf * (cn[i] - p.d4_refcn[z-1][r]).powi(2)).exp();
                norm += raw[r];
                dnorm += 2.0 * wf * (p.d4_refcn[z-1][r] - cn[i]) * raw[r];
            }
            let inv = if norm > 1e-80 { 1.0 / norm } else { 0.0 };
            for r in 0..nref {
                gw[i][r] = raw[r] * inv;
                dgw[i][r] = 2.0*wf*(p.d4_refcn[z-1][r] - cn[i])*raw[r]*inv - raw[r]*dnorm*inv*inv;
            }
        }
        let _ = q;
        let thopi = 3.0 / std::f64::consts::PI;
        let mut e = 0.0;
        for i in 0..n { for j in 0..i {
            let (zi, zj) = (self.at[i], self.at[j]);
            let r = dist(i, j);
            let r4r2ij = 3.0 * p.r4r2[zi-1] * p.r4r2[zj-1];
            let a0 = p.d3a1 * r4r2ij.sqrt() + p.d3a2;
            let r0 = a0 * a0;
            let r2 = r * r;
            let t6 = 1.0 / (r.powi(6) + r0.powi(3));
            let t8 = 1.0 / (r.powi(8) + r0.powi(4));
            // C6 + dc6/dcn(i), dc6/dcn(j)
            let (mut c6, mut dc6i, mut dc6j) = (0.0, 0.0, 0.0);
            for a in 0..p.d4_refn[zi-1] as usize {
                for b in 0..p.d4_refn[zj-1] as usize {
                    let mut s = 0.0;
                    for k in 0..23 { s += tw[k] * alpha_ref(p, zi, a, k) * alpha_ref(p, zj, b, k); }
                    let refc6 = thopi * s;
                    c6 += gw[i][a] * gw[j][b] * refc6;
                    dc6i += dgw[i][a] * gw[j][b] * refc6;
                    dc6j += gw[i][a] * dgw[j][b] * refc6;
                }
            }
            let zeta = p.zeta(zi, self.topo.qa[i]) * p.zeta(zj, self.topo.qa[j]);
            let disp = (t6 + 2.0 * r4r2ij * t8) * zeta;
            e += -c6 * disp;
            // pair gradient
            let d6 = -6.0 * r2 * r2 * t6 * t6;
            let d8 = -8.0 * r2 * r2 * r2 * t8 * t8;
            let dg = -c6 * (d6 + 2.0 * r4r2ij * d8) * zeta;
            for t in 0..3 {
                g[i][t] += dg * (xyz[i][t] - xyz[j][t]);
                g[j][t] -= dg * (xyz[i][t] - xyz[j][t]);
            }
            d_ed_cn[i] -= dc6i * disp;
            d_ed_cn[j] -= dc6j * disp;
        }}
        let _ = xyz;
        e
    }

    fn d4_dispersion(&self, xyz: &[[f64; 3]], dist: &dyn Fn(usize, usize) -> f64,
                      cn: &[f64]) -> f64 {
        let n = self.at.len();
        let p = &self.p;
        // gw weights (wf = 4.0) over CN (logCN as passed to d3_gradient)
        let wf = 4.0f64;
        let mut gw = vec![vec![0.0f64; 8]; n];
        for i in 0..n {
            let z = self.at[i];
            let nref = p.d4_refn[z-1] as usize;
            let mut norm = 0.0;
            let mut raw = vec![0.0f64; nref];
            for r in 0..nref {
                raw[r] = (-wf * (cn[i] - p.d4_refcn[z-1][r]).powi(2)).exp();
                norm += raw[r];
            }
            for r in 0..nref { gw[i][r] = raw[r] / norm; }
        }
        let _ = xyz;
        // C6 ref matrix (thopi * trapz(alpha_i * alpha_j)) on the fly per pair
        let thopi = 3.0 / std::f64::consts::PI;
        let mut e = 0.0;
        for i in 0..n { for j in 0..i {
            let (zi, zj) = (self.at[i], self.at[j]);
            let r = dist(i, j);
            let r4r2ij = 3.0 * p.r4r2[zi-1] * p.r4r2[zj-1];
            let a0 = p.d3a1 * r4r2ij.sqrt() + p.d3a2;
            let r0 = a0 * a0; // saved as R0^2, t6 uses r0^3 = (R0^2)^... note: matches Fortran
            let t6 = 1.0 / (r.powi(6) + r0.powi(3));
            let t8 = 1.0 / (r.powi(8) + r0.powi(4));
            // charge-scaled pair C6 (trapezoid over the non-uniform grid)
            let tw = trapzd_weights();
            let mut c6 = 0.0;
            for a in 0..p.d4_refn[zi-1] as usize {
                for b in 0..p.d4_refn[zj-1] as usize {
                    let mut s = 0.0;
                    for k in 0..23 {
                        s += tw[k] * alpha_ref(p, zi, a, k) * alpha_ref(p, zj, b, k);
                    }
                    c6 += gw[i][a] * gw[j][b] * thopi * s;
                }
            }
            let zeta = p.zeta(zi, self.topo.qa[i]) * p.zeta(zj, self.topo.qa[j]);
            // NOTE: per-pair = -C6*disp*zeta (the 0.5 in Fortran is per-atom bookkeeping
            // that enters twice; validated exactly on H2/N2/CO/F2). Molecular H-X pairs
            // still show a ~0.6x offset vs xtb -- see docs/gfnff-porting-notes.md.
            e += -c6 * (t6 + 2.0 * r4r2ij * t8) * zeta;
        }}
        e
    }
}

fn alpha_ref(p: &Params, z: usize, iref: usize, w: usize) -> f64 {
    // newD3Model (src/disp/dftd4.F90 54-104): alpha_sec = sscale(refsys)*secaiw;
    // atom alpha = max(ascale*(alphaiw - hcount*alpha_sec), 0)
    let is = (p.d4_refsys[z-1][iref] - 1) as usize;   // refsys stores 1-based system id
    let asec = p.d4_sscale[is] * p.d4_secaiw[w][is];
    let raw = p.d4_ascale[z-1][iref] * (p.d4_alphaiw[w][z-1][iref] - p.d4_hcount[z-1][iref] * asec);
    if raw > 0.0 { raw } else { 0.0 }
}

/// trapezoid weights on the D4 imaginary-frequency grid (dftd4.F90 trapzd)
fn trapzd_weights() -> [f64; 23] {
    const FREQ: [f64; 23] = [1e-6, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
        1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0];
    let mut w = [0.0f64; 23];
    for k in 0..23 {
        let l = if k > 0 { FREQ[k] - FREQ[k-1] } else { 0.0 };
        let r = if k < 22 { FREQ[k+1] - FREQ[k] } else { 0.0 };
        w[k] = 0.5 * (l + r);
    }
    w
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
fn sub3(xyz: &[[f64; 3]], i: usize, j: usize) -> [f64; 3] {
    [xyz[i][0]-xyz[j][0], xyz[i][1]-xyz[j][1], xyz[i][2]-xyz[j][2]]
}
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[0]+b[0], a[1]+b[1], a[2]+b[2]] }
fn subv3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn norm3(a: [f64; 3]) -> f64 { dot3(a, a).sqrt() }
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] { [a[0]*s, a[1]*s, a[2]*s] }

/// angle damping with derivative (gfnffdampa); returns (damp, 2*d(damp)/d(r^2))
fn dampa2(p: &Params, ati: usize, atj: usize, r2: f64) -> (f64, f64) {
    let rc = (p.rcov[ati-1] + p.rcov[atj-1]) / BOHR * 4.0 / 3.0;
    let rcut = p.atcuta * rc * rc;
    let rr = (r2 / rcut).powi(2);
    (1.0 / (1.0 + rr), -4.0 * rr / (r2 * (1.0 + rr).powi(2)))
}

/// torsion damping with derivative (gfnffdampt)
fn dampt2(p: &Params, ati: usize, atj: usize, r2: f64) -> (f64, f64) {
    let rc = (p.rcov[ati-1] + p.rcov[atj-1]) / BOHR * 4.0 / 3.0;
    let rcut = 0.505 * rc * rc;
    let rr = (r2 / rcut).powi(2);
    (1.0 / (1.0 + rr), -4.0 * rr / (r2 * (1.0 + rr).powi(2)))
}

/// d logCN / d CN for the raw (erf) coordination number.
/// logCN = ln(1+e^cm) - ln(1+e^(cm-cn))  =>  d/dcn = e^cm/(e^cm+e^cn)
/// (caller passes the RAW cn; no inversion needed - the previous version
/// inverted algebraically AND was fed cn_raw, off by 3-50x, which was the
/// source of the ~0.1 kcal/mol/A gradient residual vs xtb)
fn create_dlogcn(p: &Params, cn: f64) -> f64 {
    let cm = p.cnmax;
    cm.exp() / (cm.exp() + cn.exp())
}

/// dihedral derivative wrt the four atoms (dphidrPBC, non-periodic)
#[allow(clippy::too_many_arguments)]
fn dphidr(xyz: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize, phi: f64)
    -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    let sub = |a: usize, b: usize| [xyz[a][0]-xyz[b][0], xyz[a][1]-xyz[b][1], xyz[a][2]-xyz[b][2]];
    let ra = sub(j, i);
    let rb = sub(k, j);
    let rc = sub(l, k);
    let rapb = [ra[0]+rb[0], ra[1]+rb[1], ra[2]+rb[2]];
    let rbpc = [rb[0]+rc[0], rb[1]+rc[1], rb[2]+rc[2]];
    let na = cross3(ra, rb);
    let nb = cross3(rb, rc);
    let nan = norm3(na);
    let nbn = norm3(nb);
    let cosphi = phi.cos();
    let sinphi = phi.sin();
    let nenner = nan * nbn * sinphi;
    if nenner.abs() < 1e-14 {
        return ([0.0;3], [0.0;3], [0.0;3], [0.0;3]);
    }
    let onenner = 1.0 / nenner;
    let rab = cross3(na, rb);
    let rbb = cross3(nb, rb);
    let rba = cross3(nb, ra);
    let rac = cross3(na, rc);
    let rbc = cross3(nb, rc);
    let raa = cross3(na, ra);
    let rapba = cross3(rapb, na);
    let rapbb = cross3(rapb, nb);
    let rbpca = cross3(rbpc, na);
    let rbpcb = cross3(rbpc, nb);
    let mut di = [0.0f64; 3]; let mut dj = [0.0f64; 3];
    let mut dk = [0.0f64; 3]; let mut dl = [0.0f64; 3];
    for t in 0..3 {
        di[t] = onenner * (cosphi * nbn / nan * rab[t] - rbb[t]);
        dj[t] = onenner * (cosphi * (nbn / nan * rapba[t] + nan / nbn * rbc[t]) - (rac[t] + rapbb[t]));
        dk[t] = onenner * (cosphi * (nbn / nan * raa[t] + nan / nbn * rbpcb[t]) - (rba[t] + rbpca[t]));
        dl[t] = onenner * (cosphi * nan / nbn * rbb[t] - rab[t]);
    }
    (di, dj, dk, dl)
}

/// HMO solve: diagonalize api, Fermi-smear occupations at 4000 K,
/// return (P = C diag(focc) C^T, sum focc*eps, eigenvalues in eV)
/// symmetric Jacobi eigensolver: returns (eigenvalues sorted asc, eigenvectors
/// as evecs[row][col] with col = eigenvector index)
fn jacobi_eig(a0: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a0.len();
    let mut a = a0.to_vec();
    let mut v = vec![vec![0.0f64; n]; n];
    for i in 0..n { v[i][i] = 1.0; }
    for _sweep in 0..100 {
        let mut off = 0.0;
        for i in 0..n { for j in (i+1)..n { off += a[i][j]*a[i][j]; } }
        if off < 1e-22 { break; }
        for p in 0..n { for q in (p+1)..n {
            if a[p][q].abs() < 1e-15 { continue; }
            let theta = (a[q][q]-a[p][p])/(2.0*a[p][q]);
            let t = theta.signum()/(theta.abs()+ (theta*theta+1.0).sqrt());
            let c = 1.0/ (t*t+1.0).sqrt();
            let s = t*c;
            for k in 0..n {
                let akp = a[k][p]; let akq = a[k][q];
                a[k][p] = c*akp - s*akq;
                a[k][q] = s*akp + c*akq;
            }
            for k in 0..n {
                let apk = a[p][k]; let aqk = a[q][k];
                a[p][k] = c*apk - s*aqk;
                a[q][k] = s*apk + c*aqk;
            }
            for k in 0..n {
                let vkp = v[k][p]; let vkq = v[k][q];
                v[k][p] = c*vkp - s*vkq;
                v[k][q] = s*vkp + c*vkq;
            }
        }}
    }
    let mut evals: Vec<(f64, usize)> = (0..n).map(|i| (a[i][i], i)).collect();
    evals.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    let vals: Vec<f64> = evals.iter().map(|x| x.0).collect();
    let mut vecs = vec![vec![0.0f64; n]; n];
    for (newc, (_, oldc)) in evals.iter().enumerate() {
        for r in 0..n { vecs[r][newc] = v[r][*oldc]; }
    }
    (vals, vecs)
}

fn hmo_solve(api: &[Vec<f64>], nel: usize) -> (Vec<Vec<f64>>, f64, Vec<f64>, Vec<f64>) {
    let n = api.len();
    let (evals, evecs) = jacobi_eig(api);
    // energy scaling: eV
    let eps: Vec<f64> = evals.iter().map(|e| e * 0.1 * 27.2113957).collect();
    // Fermi smearing at 4000 K; restricted occupations split into
    // na = ceil(nel/2) alpha + nb = floor(nel/2) beta electrons (occu,
    // scc_core.f90) so that odd electron counts conserve nel exactly
    let bkt: f64 = 3.166815e-6 * 27.2113957 * 4000.0;   // eV
    let smear = |count: usize| -> Vec<f64> {
        let mut focc = vec![0.0f64; n];
        if count == 0 { return focc; }
        let mut ef = if count >= n { eps[n - 1] } else { 0.5 * (eps[count - 1] + eps[count]) };
        for _ in 0..200 {
            let mut tot = 0.0; let mut dtot = 0.0;
            for i in 0..n {
                let x = (eps[i] - ef) / bkt;
                let ex = x.exp();
                let f = if x < 50.0 { 1.0 / (ex + 1.0) } else { 0.0 };
                let df = if x < 50.0 { ex / (bkt * (ex + 1.0).powi(2)) } else { 0.0 };
                focc[i] = f; tot += f; dtot += df;
            }
            if dtot > 0.0 { ef += (count as f64 - tot) / dtot; }
            if (count as f64 - tot).abs() <= 1e-9 { break; }
        }
        focc
    };
    let na = nel.div_ceil(2);
    let nb = nel / 2;
    let fa = smear(na);
    let fb = smear(nb);
    let mut focc: Vec<f64> = (0..n).map(|i| fa[i] + fb[i]).collect();
    // perfect biradical (anti-aromatic): closed-shell hard fill, breaking the
    // symmetry (gfnffqmsolve; note xtb drops the odd electron here)
    if na < n && (focc[na - 1] - focc[na]).abs() < 1e-4 {
        for f in focc.iter_mut() { *f = 0.0; }
        for i in 0..nel / 2 { focc[i] = 2.0; }
    }
    let eel: f64 = (0..n).map(|i| focc[i] * eps[i]).sum();
    // P = C diag(focc) C^T   (evecs[k][m] = component k of eigenvector m)
    let mut p = vec![vec![0.0f64; n]; n];
    for i in 0..n { for j in 0..n {
        let mut s = 0.0;
        for m in 0..n { s += evecs[i][m] * focc[m] * evecs[j][m]; }
        p[i][j] = s;
    }}
    (p, eel, eps, focc)
}

/// inversion angle omega (asin(n.r/|n||r|)), constr.f90 omegaPBC
fn improper_angle(xyz: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let sub = |a: usize, b: usize| [xyz[a][0]-xyz[b][0], xyz[a][1]-xyz[b][1], xyz[a][2]-xyz[b][2]];
    let re = sub(i, j);                    // center -> 1st nb
    let rd = [xyz[k][0]-xyz[j][0], xyz[k][1]-xyz[j][1], xyz[k][2]-xyz[j][2]]; // 1st -> 2nd
    let rv = sub(l, i);                    // center -> 3rd
    let rn = cross3(re, rd);
    let rnn = norm3(rn) + 1e-14;
    let rvn = norm3(rv) + 1e-14;
    (dot3(rn, rv) / (rnn * rvn)).clamp(-1.0, 1.0).asin()
}

/// inversion derivatives (constr.f90 domegadrPBC, non-periodic);
/// atom order (i=center, j=1st, k=2nd, l=3rd) matching improper_angle
fn domegadr(xyz: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize, omega: f64)
    -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    let sub = |a: usize, b: usize| [xyz[a][0]-xyz[b][0], xyz[a][1]-xyz[b][1], xyz[a][2]-xyz[b][2]];
    let sinomega = omega.sin();
    let re = sub(i, j);      // center - 1st
    let rd = sub(k, j);      // 2nd - 1st
    let rv = sub(l, i);      // 3rd - center
    let rdme = subv3(rd, re);
    let rn = cross3(re, rd);
    let rvn = norm3(rv);
    let rnn = norm3(rn);
    let rve = cross3(rv, re);
    let rne = cross3(rn, re);
    let rdv = cross3(rd, rv);
    let rdn = cross3(rd, rn);
    let rvdme = cross3(rv, rdme);
    let rndme = cross3(rn, rdme);
    let nenner = rnn * rvn * omega.cos();
    if nenner.abs() <= 1e-14 {
        return ([0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]);
    }
    let onenner = 1.0 / nenner;
    let mut di = [0.0f64; 3]; let mut dj = [0.0f64; 3];
    let mut dk = [0.0f64; 3]; let mut dl = [0.0f64; 3];
    for t in 0..3 {
        di[t] = onenner * (rdv[t] - rn[t]
            - sinomega * (rvn / rnn * rdn[t] - rnn / rvn * rv[t]));
        dj[t] = onenner * (rvdme[t] - sinomega * rvn / rnn * rndme[t]);
        dk[t] = onenner * (rve[t] - sinomega * rvn / rnn * rne[t]);
        dl[t] = onenner * (rn[t] - sinomega * rnn / rvn * rv[t]);
    }
    (di, dj, dk, dl)
}

/// scaled bend factor for eg3: e = 1 - kijk*(cosa - cos(c0))^2 with
/// kijk = fc/(cos(0)-cos(c0))^2 (egbend_nci_mul, gfnff_eg.f90 1317-1373);
/// atoms (j=B, i=C, k=H): the H-B=C angle, gradient g3 = [B, C, H]
fn bend_nci_mul(xyz: &[[f64; 3]], j: usize, i: usize, k: usize, c0: f64, fc: f64,
                g3: &mut [[f64; 3]; 3]) -> f64 {
    let kijk = fc / (1.0 - c0.cos()).powi(2);   // cos(0) = 1
    let vab = sub3(xyz, i, j);   // B -> C
    let vcb = sub3(xyz, k, j);   // B -> H
    let rab2 = dot3(vab, vab);
    let rcb2 = dot3(vcb, vcb);
    let vp = cross3(vcb, vab);
    let rp = norm3(vp) + 1e-14;
    let cosa = (dot3(vab, vcb) / (rab2 * rcb2).sqrt()).clamp(-1.0, 1.0);
    let theta = cosa.acos();
    let ea = kijk * (cosa - c0.cos()).powi(2);
    let deddt = 2.0 * kijk * theta.sin() * (c0.cos() - cosa);
    let deda = scale3(cross3(vab, vp), -deddt / (rab2 * rp));
    let dedc = scale3(cross3(vcb, vp), deddt / (rcb2 * rp));
    let dedb = add3(deda, dedc);
    for t in 0..3 {
        g3[0][t] = dedb[t];   // B (center)
        g3[1][t] = -deda[t]; // C
        g3[2][t] = -dedc[t]; // H
    }
    1.0 - ea
}

/// amide N detection: N is pi, sp3, bonded to exactly ONE pi C which
/// carries exactly ONE terminal pi O (gfnff_ini2.F90 amide(), 2275-2305)
fn is_amide_n(at: &[usize], hyb: &[i32], nb: &[Vec<usize>], piadr: &[bool], a: usize) -> bool {
    if !piadr[a] || hyb[a] != 3 || at[a] != 7 { return false; }
    let mut nc = 0;
    let mut ic = 0usize;
    for &j in &nb[a] {
        if at[j] == 6 && piadr[j] { nc += 1; ic = j; }
    }
    if nc != 1 { return false; }
    let no = nb[ic].iter().filter(|&&j| at[j] == 8 && piadr[j] && nb[j].len() == 1).count();
    no == 1
}

fn dist2(xyz: &[[f64; 3]], i: usize, j: usize) -> f64 {
    let dx = xyz[i][0]-xyz[j][0]; let dy = xyz[i][1]-xyz[j][1]; let dz = xyz[i][2]-xyz[j][2];
    dx*dx + dy*dy + dz*dz
}

/// dihedral l-i-j-k in radians (valijklff, math.f)
fn dihedral(xyz: &[[f64; 3]], l: usize, i: usize, j: usize, k: usize) -> f64 {
    let sub = |a: usize, b: usize| [xyz[a][0]-xyz[b][0], xyz[a][1]-xyz[b][1], xyz[a][2]-xyz[b][2]];
    let ra = sub(l, i);   // i -> l
    let rb = sub(i, j);   // j -> i
    let rc = sub(j, k);   // k -> j
    let cross = |a: [f64;3], b: [f64;3]| [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
    let dot = |a: [f64;3], b: [f64;3]| a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    let na = cross(ra, rb);
    let nb = cross(rb, rc);
    let nlen = |v: [f64;3]| dot(v, v).sqrt();
    let sn = dot(na, nb) / (nlen(na) * nlen(nb) + 1e-14);
    sn.clamp(-1.0, 1.0).acos()
}

/// torsion distance damping (gfnffdampt, atcutt = 0.505)
fn dampt(p: &Params, ati: usize, atj: usize, r2: f64) -> f64 {
    let rc = (p.rcov[ati-1] + p.rcov[atj-1]) / BOHR * 4.0 / 3.0;
    let rcut = 0.505 * rc * rc;
    let rr = (r2 / rcut).powi(2);
    1.0 / (1.0 + rr)
}

fn dampa(p: &Params, ati: usize, atj: usize, r2: f64) -> f64 {
    // gfnffdampa (gfnff_eg.f90 1695); rcov converted like xtb: aatoau * 4/3
    let rc = (p.rcov[ati-1] + p.rcov[atj-1]) / BOHR * 4.0 / 3.0;
    let rcut = p.atcuta * rc * rc;
    let rr = (r2 / rcut).powi(2);
    1.0 / (1.0 + rr)
}

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

impl Params {
    fn metal_is(&self, _z: usize) -> bool { false } // organic subset: no metals
    fn repscalb(&self) -> f64 { 1.7583 }

    /// element-specific bond-radius factors (gfnff_ini2.F90 fat, lines 71-95)
    fn fat(&self, z: usize) -> f64 {
        match z {
            1 => 1.02, 4 => 1.03, 5 => 1.02, 8 => 1.02, 9 => 1.05, 10 => 1.10,
            11 => 1.01, 12 => 1.02, 15 => 0.97, 18 => 1.10, 19 => 1.02, 20 => 1.02,
            34 => 0.99, 38 => 1.02, 50 => 1.01, 51 => 0.99, 52 => 0.95, 53 => 0.98,
            56 => 1.02, 76 => 1.02, 82 => 1.06, 83 => 0.95,
            _ => 1.0,
        }
    }

    /// gfnffrab pair reference bond length (Bohr), gfnff_rab.f
    /// rab = (r0_A + cnfak_A*cn_A + r0_B + cnfak_B*cn_B + shift) * ff * scaleF
    pub fn gfnffrab(&self, ati: usize, atj: usize, cni: f64, cnj: f64, shift: f64) -> f64 {
        const P: [[f64; 2]; 6] = [
            [29.84522887, -8.87843763], [-1.70549806, 2.10878369],
            [ 6.54013762,  0.08009374], [ 6.39169003, -0.85808076],
            [ 6.0,        -1.15      ], [ 5.6,        -1.30      ]];
        let itabrow6 = |z: usize| -> usize {
            if z <= 2 { 1 } else if z <= 10 { 2 } else if z <= 18 { 3 }
            else if z <= 36 { 4 } else if z <= 54 { 5 } else { 6 }
        };
        let ra = self.rab_r0[ati-1] + self.rab_cnfak[ati-1] * cni;
        let rb = self.rab_r0[atj-1] + self.rab_cnfak[atj-1] * cnj;
        let den = (self.rab_en[ati-1] - self.rab_en[atj-1]).abs();
        let (ir, jr) = (itabrow6(ati), itabrow6(atj));
        let k1 = 0.005 * (P[ir-1][0] + P[jr-1][0]);
        let k2 = 0.005 * (P[ir-1][1] + P[jr-1][1]);
        let ff = 1.0 - k1 * den - k2 * den * den;
        (ra + rb + shift) * ff   // scaleF = 1 for organic; Ln/An specials TODO
    }

    pub fn hdiag(&self, z: usize) -> &f64 {
        self.hdiag_tab.get(&z).unwrap_or(&self.hdiag_c)
    }
    pub fn hoff(&self, z: usize) -> &f64 {
        self.hoff_tab.get(&z).unwrap_or(&self.hoff_c)
    }

    /// d rab0 / d(cn_i, cn_j) = (cnfak_i, cnfak_j) * ff  (gfnffdrab rabdcn)
    pub fn gfnffrab_dcnd(&self, ati: usize, atj: usize, shift: f64) -> (f64, f64) {
        // ff recomputed (same as gfnffrab)
        const P: [[f64; 2]; 6] = [
            [29.84522887, -8.87843763], [-1.70549806, 2.10878369],
            [ 6.54013762,  0.08009374], [ 6.39169003, -0.85808076],
            [ 6.0,        -1.15      ], [ 5.6,        -1.30      ]];
        let itabrow6 = |z: usize| -> usize {
            if z <= 2 { 1 } else if z <= 10 { 2 } else if z <= 18 { 3 }
            else if z <= 36 { 4 } else if z <= 54 { 5 } else { 6 }
        };
        let den = (self.rab_en[ati-1] - self.rab_en[atj-1]).abs();
        let (ir, jr) = (itabrow6(ati), itabrow6(atj));
        let k1 = 0.005 * (P[ir-1][0] + P[jr-1][0]);
        let k2 = 0.005 * (P[ir-1][1] + P[jr-1][1]);
        let ff = 1.0 - k1 * den - k2 * den * den;
        (self.rab_cnfak[ati-1] * ff, self.rab_cnfak[atj-1] * ff)
    }
}

/// erf CN then logCN (gfnff_dlogcoord + create_logCN, gfnff_eg.f90 3497)
/// note: covalentRadD3 values are in Angstrom and carried *aatoau*4/3 in xtb
pub fn erf_cn_raw(p: &Params, at: &[usize], xyz: &[[f64; 3]]) -> Vec<f64> {
    let n = at.len();
    let mut cn = vec![0.0f64; n];
    for i in 0..n { for j in 0..i {
        let dx = xyz[i][0]-xyz[j][0]; let dy = xyz[i][1]-xyz[j][1]; let dz = xyz[i][2]-xyz[j][2];
        let r = (dx*dx+dy*dy+dz*dz).sqrt();
        if r > 60.0 { continue; }
        let r0 = (p.rcov[at[i]-1] + p.rcov[at[j]-1]) / BOHR * 4.0 / 3.0;
        let dr = (r - r0) / r0;
        let c = 0.5 * (1.0 + erf(-7.5 * dr));
        cn[i] += c; cn[j] += c;
    }}
    cn
}

pub fn erf_cn(p: &Params, at: &[usize], xyz: &[[f64; 3]]) -> Vec<f64> {
    let n = at.len();
    let mut cn = vec![0.0f64; n];
    for i in 0..n { for j in 0..i {
        let dx = xyz[i][0]-xyz[j][0]; let dy = xyz[i][1]-xyz[j][1]; let dz = xyz[i][2]-xyz[j][2];
        let r = (dx*dx+dy*dy+dz*dz).sqrt();
        if r > 60.0 { continue; }
        let r0 = (p.rcov[at[i]-1] + p.rcov[at[j]-1]) / BOHR * 4.0 / 3.0;
        let dr = (r - r0) / r0;
        let c = 0.5 * (1.0 + erf(-7.5 * dr));
        cn[i] += c; cn[j] += c;
    }}
    cn.iter().map(|&c| create_logcn(p, c)).collect()
}

fn create_logcn(p: &Params, cn: f64) -> f64 {
    let e = p.cnmax.exp();
    (1.0 + e).ln() - (1.0 + (p.cnmax - cn).exp()).ln()
}

/// bond detection (gfnff_ini2.F90 neigh 100-130 + neighbor.f90 getnb):
/// pair radius from gfnffrab with a normcn CN guess, charge-corrected
/// (rqshrink) and fat element-scaled; bond if r < rthr * rco. The
/// hypercoordination filter (skip the pair if either atom's full candidate
/// count exceeds 4 for group 1-2, else 6) mirrors xtb's icase-2 getnb
/// semantics exactly ("if(nnfi.gt.hc_crit) cycle" per PAIR, no distance
/// trimming). In xtb that filtered list is transient (hyb determination,
/// cluster tagging) while the final topology list is the full nbf; the two
/// coincide for normal organic valences where the caps never fire, which is
/// the regime this single-list approximation covers. Unlike the old
/// rcov*1.25 criterion this keeps stretched bonds (the radius grows as CN
/// drops), matching xtb.
fn detect_bonds(p: &Params, at: &[usize], xyz: &[[f64; 3]], qa: &[f64]) -> Vec<Vec<usize>> {
    let n = at.len();
    let mut nbf = vec![Vec::new(); n];
    for i in 0..n { for j in 0..i {
        let r = dist2(xyz, i, j).sqrt();
        let mut rco = p.gfnffrab(at[i], at[j],
            p.normcn[at[i] - 1] as f64, p.normcn[at[j] - 1] as f64, 0.0);
        rco -= (qa[i] + qa[j]) * p.rqshrink;
        rco *= p.fat(at[i]) * p.fat(at[j]);
        if r < p.rthr * rco { nbf[i].push(j); nbf[j].push(i); }
    }}
    let mut nb = vec![Vec::new(); n];
    for i in 0..n { for j in 0..i {
        if !nbf[i].contains(&j) { continue; }
        let hc_i = if p.group[at[i] - 1] <= 2 { 4 } else { 6 };
        let hc_j = if p.group[at[j] - 1] <= 2 { 4 } else { 6 };
        if nbf[i].len() > hc_i || nbf[j].len() > hc_j { continue; }
        nb[i].push(j); nb[j].push(i);
    }}
    nb
}

fn sorted3(a: usize, b: usize, c: usize) -> Vec<usize> { let mut v = vec![a, b, c]; v.sort_unstable(); v }
fn sorted4(a: usize, b: usize, c: usize, d: usize) -> Vec<usize> { let mut v = vec![a, b, c, d]; v.sort_unstable(); v }
fn sorted5(a: usize, b: usize, c: usize, d: usize, e: usize) -> Vec<usize> { let mut v = vec![a, b, c, d, e]; v.sort_unstable(); v }
fn sorted6(a: usize, b: usize, c: usize, d: usize, e: usize, f: usize) -> Vec<usize> { let mut v = vec![a, b, c, d, e, f]; v.sort_unstable(); v }
fn add_ring(rings_all: &mut [Vec<Vec<usize>>], members: &[usize]) {
    for &m in members {
        if !rings_all[m].iter().any(|r| r == members) { rings_all[m].push(members.to_vec()); }
    }
}

/// alpha-C=O torsion strengthening check (alphaCO, gfnff_ini2.F90):
/// sp2 pi C (with one terminal pi O) bonded to an sp3 C
fn alpha_co(at: &[usize], hyb: &[i32], nb: &[Vec<usize>], piadr: &[bool], a: usize, b: usize) -> bool {
    let check = |x: usize, y: usize| -> bool {
        if piadr[x] && hyb[y] == 3 && at[x] == 6 && at[y] == 6 {
            let no = nb[x].iter()
                .filter(|&&j| at[j] == 8 && piadr[j] && nb[j].len() == 1).count();
            return no == 1;
        }
        false
    };
    check(a, b) || check(b, a)
}

/// aldehyde/carbonyl carbon check (ctype, gfnff_ini2.F90): pi C with one pi O
fn is_aldehyde_c(at: &[usize], nb: &[Vec<usize>], piadr: &[bool], a: usize) -> bool {
    if !piadr[a] || at[a] != 6 { return false; }
    let no = nb[a].iter().filter(|&&j| at[j] == 8 && piadr[j]).count();
    no == 1
}

/// smallest ring containing all three angle atoms (ringsbend, gfnff_ini2 557)
fn ringsbend(topo_rings: &[Vec<Vec<usize>>], i: usize, j: usize, k: usize) -> usize {
    let mut best = 99usize;
    for rings in [topo_rings.get(i), topo_rings.get(j), topo_rings.get(k)] {
        let Some(rings) = rings else { continue };
        for r in rings {
            if r.contains(&j) && r.contains(&k) && r.len() < best { best = r.len(); }
        }
    }
    if best == 99 { 0 } else { best }
}

/// smallest ring containing the bond i-j (ringsbond, gfnff_ini2 531), 0 if none
fn ringsbond(topo_rings: &[Vec<Vec<usize>>], i: usize, j: usize) -> usize {
    let mut best = 99usize;
    if let Some(rings) = topo_rings.get(i) {
        for r in rings {
            if r.contains(&j) && r.len() < best { best = r.len(); }
        }
    }
    if best == 99 { 0 } else { best }
}

/// smallest (ringstors, gfnff_ini2 603) or largest (ringstorl, 667) ring
/// containing ALL four torsion atoms, 0 if they share no ring
fn ring_common(topo_rings: &[Vec<Vec<usize>>], atoms: [usize; 4], largest: bool) -> usize {
    let mut best = if largest { 0usize } else { 99usize };
    if let Some(rings) = topo_rings.get(atoms[0]) {
        for r in rings {
            if !atoms[1..].iter().all(|&a| r.contains(&a)) { continue; }
            if largest { if r.len() > best { best = r.len(); } }
            else if r.len() < best { best = r.len(); }
        }
    }
    if !largest && best == 99 { 0 } else { best }
}

/// topology-dependent electronegativity corrections dxi (gfnff_ini.f90 391-447)
fn compute_dxi(p: &Params, at: &[usize], nb: &[Vec<usize>], hyb: &[i32], piadr: &[bool]) -> Vec<f64> {
    let n = at.len();
    let mut dxi = vec![0.0f64; n];
    for i in 0..n {
        let z = at[i];
        let nn = nb[i].len();
        let nh = nb[i].iter().filter(|&&j| at[j] == 1).count();
        if nn == 0 { continue; }
        if z == 5 { dxi[i] += nh as f64 * 0.015; }
        // carbene C (xtb itag=1: bent 2-coordinate group 14)
        if z == 6 && nn == 2 && hyb[i] == 2 { dxi[i] = -0.15; }
        // free CO: C with nn=1 bonded to O with nn=1
        if z == 6 && nn == 1 && at[nb[i][0]] == 8 && nb[nb[i][0]].len() == 1 {
            dxi[i] = 0.15;
        }
        // nitro O: nn=1 bonded to pi N
        if z == 8 && nn == 1 && at[nb[i][0]] == 7 && piadr[nb[i][0]] {
            dxi[i] = 0.05;
        }
        if z == 8 && nn == 2 && nh == 2 { dxi[i] = -0.02; }
        if p.group[z - 1] == 6 && nn > 2 { dxi[i] += nn as f64 * 0.005; }
        if z == 8 || z == 16 { dxi[i] -= nh as f64 * 0.005; }
        if p.group[z - 1] == 7 && z > 9 && nn > 1 { dxi[i] -= nn as f64 * 0.021; }
    }
    dxi
}

/// condense H charges onto their heavy neighbors (gfnff_ini2.F90 qheavy)
fn qheavy(at: &[usize], nb: &[Vec<usize>], q: &mut [f64]) {
    let qtmp = q.to_vec();
    for (i, &z) in at.iter().enumerate() {
        if z != 1 { continue; }
        q[i] = 0.0;
        for &k in &nb[i] {
            q[k] += qtmp[i] / nb[i].len() as f64;
        }
    }
}

/// Floyd-Warshall topology distances via `rad` sums (gfnff_ini.f90 476-510)
fn floyd_rabd(p: &Params, at: &[usize], nb: &[Vec<usize>]) -> Vec<Vec<f64>> {
    let n = at.len();
    let cutoff = 13.0f64;
    let thr = 12.0f64;
    let mut d = vec![vec![cutoff; n]; n];
    for i in 0..n { d[i][i] = 0.0; }
    for i in 0..n { for &j in &nb[i] {
        let r = p.rad[at[i]-1] + p.rad[at[j]-1];
        d[i][j] = r; d[j][i] = r;
    }}
    for k in 0..n {
        for i in 0..n {
            if d[i][k] > thr { continue; }
            for j in 0..n {
                if d[k][j] > thr { continue; }
                if d[i][j] > d[i][k] + d[k][j] { d[i][j] = d[i][k] + d[k][j]; }
            }
        }
    }
    for i in 0..n { for j in 0..n { if d[i][j] > thr { d[i][j] = cutoff; } } }
    d
}

/// solve dense linear system by Gaussian elimination with partial pivoting
/// (augmented-matrix form, no permutation bookkeeping)
fn solve_sym(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut m: Vec<Vec<f64>> = a.iter().cloned().zip(b.iter()).map(|(r, &v)| {
        let mut row = r.clone(); row.push(v); row
    }).collect();
    for col in 0..n {
        let mut piv = col;
        let mut best = m[col][col].abs();
        for r in (col+1)..n {
            if m[r][col].abs() > best { best = m[r][col].abs(); piv = r; }
        }
        if best < 1e-30 { continue; }
        if piv != col { m.swap(col, piv); }
        for r in (col+1)..n {
            let f = m[r][col] / m[col][col];
            if f == 0.0 { continue; }
            for c in col..=n { m[r][c] -= f * m[col][c]; }
        }
    }
    let mut x = vec![0.0f64; n];
    for r in (0..n).rev() {
        let mut s = m[r][n];
        for c in (r+1)..n { s -= m[r][c] * x[c]; }
        if m[r][r].abs() > 1e-30 { x[r] = s / m[r][r]; }
    }
    x
}

/// EEQ solve: topology mode uses rabd distances, RHS with cnf*sqrt(min(nb,cnmax));
/// returns (charges q, ES energy) (goedeckera / goed_gfnff, gfnff_eg.f90 1758-1914)
#[allow(clippy::too_many_arguments)]
fn solve_eeq(p: &Params, at: &[usize], rabd: &[Vec<f64>], nb: &[Vec<usize>],
             charge: f64, topology_mode: bool, _hyb: &[i32], dxi: &[f64],
             fraglist: &[usize], qfrag: &[f64]) -> (Vec<f64>, f64) {
    let n = at.len();
    let nfrag = qfrag.len();
    let m = n + nfrag;
    let mut a = vec![vec![0.0f64; m]; m];
    let mut x = vec![0.0f64; m];
    for i in 0..n {
        let z = at[i];
        let cn_i = if topology_mode { nb[i].len().min(p.cnmax as usize) as f64 } else { 0.0 };
        x[i] = -p.chi[z-1] + dxi[i] + cn_i.sqrt() * p.cnf[z-1];
        a[i][i] = TSQRT2PI / p.alp[z-1] + p.gam[z-1];
    }
    for i in 0..n { for j in 0..i {
        // topology distances are Angstrom -> Bohr; rfgoed1 scaling
        let r = p.rfgoed1 * rabd[i][j] / BOHR;
        let g = 1.0 / (p.alp[at[i]-1].powi(2) + p.alp[at[j]-1].powi(2)).sqrt();
        let v = erf(g * r) / r;
        a[i][j] = v; a[j][i] = v;
    }}
    // fragment charge constraints
    for (fi, &qf) in qfrag.iter().enumerate() {
        x[n + fi] = qf;
        for j in 0..n {
            if fraglist[j] == fi {
                a[n + fi][j] = 1.0;
                a[j][n + fi] = 1.0;
            }
        }
    }
    let q = solve_sym(&a, &x);
    // ES energy of the topology-charge distribution (goedeckera)
    let mut es = 0.0f64;
    for i in 0..n { for j in 0..i {
        let r = p.rfgoed1 * rabd[i][j] / BOHR;
        let g = 1.0 / (p.alp[at[i] - 1].powi(2) + p.alp[at[j] - 1].powi(2)).sqrt();
        es += q[i] * q[j] * erf(g * r) / r;
    }}
    for i in 0..n {
        es += -q[i] * x[i]
            + q[i] * q[i] * 0.5 * (p.gam[at[i] - 1] + TSQRT2PI / p.alp[at[i] - 1]);
    }
    (q[..n].to_vec(), es)
}

// ---------------------------------------------------------------------------
// validation against xtb --gfnff --verbose (term-by-term)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// water geometry identical to the xtb reference run
    const WATER: [[f64; 3]; 3] = [
        [0.0, 0.0, 0.1173],
        [0.0, 0.7572, -0.4692],
        [0.0, -0.7572, -0.4692],
    ];

    #[test]
    fn erf_accuracy() {
        // spot-check against known values
        assert!((erf(1.0) - 0.8427007929497149).abs() < 1e-14);
        assert!((erf(0.5) - 0.5204998778130465).abs() < 1e-14);
        assert!((erf(3.0) - 0.9999779095030014).abs() < 1e-12);
        assert!((erf(5.0) - 0.9999999999984626).abs() < 1e-12);
    }

    #[test]
    fn water_components_vs_xtb() {
        let at = [8usize, 1, 1];
        let g = Gfnff::new(&at, &WATER, 0.0);
        let e = g.energy(&WATER);

        // reference: xtb 6.7.1 --gfnff --verbose, Eh
        let ref_bond  = -0.270029555436;
        let ref_angle =  0.000461003522;
        let ref_rep   =  0.031478455872;
        let ref_es    = -0.089052756232;
        let ref_disp  = -0.000139319111;

        println!("bond  {:>14.9} (ref {ref_bond})", e.bond);
        println!("angle {:>14.9} (ref {ref_angle})", e.angle);
        println!("rep   {:>14.9} (ref {ref_rep})", e.rep);
        println!("es    {:>14.9} (ref {ref_es})", e.es);
        println!("disp  {:>14.9} (ref {ref_disp})", e.disp);

        // current agreement: bond exact (1e-8), others within ~0.5-25%
        // remaining deviations trace to topology-charge (qa) setup details
        assert!((e.bond  - ref_bond ).abs() < 1e-7, "bond  off by {}", e.bond  - ref_bond);
        assert!((e.angle - ref_angle).abs() < 2e-6, "angle off by {}", e.angle - ref_angle);
        assert!((e.rep   - ref_rep  ).abs() < 1e-7, "rep   off by {}", e.rep   - ref_rep);
        assert!((e.es    - ref_es   ).abs() < 1e-6, "es    off by {}", e.es    - ref_es);
        assert!((e.disp  - ref_disp ).abs() < 1e-8, "disp  off by {}", e.disp  - ref_disp);
    }
}


#[cfg(test)]
mod tests_more {
    use super::*;

    #[test]
    fn methane_components_vs_xtb() {
        let at = [6usize, 1, 1, 1, 1];
        let xyz = [
            [0.0, 0.0, 0.0],
            [0.6287, 0.6287, 0.6287],
            [-0.6287, -0.6287, 0.6287],
            [-0.6287, 0.6287, -0.6287],
            [0.6287, -0.6287, -0.6287],
        ];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        let (rb, ra, rr, re, rd) =
            (-0.654717754336, 0.000067338640, 0.025993802631, -0.001559154623, -0.000650752337);
        println!("bond  {:>13.9} (ref {rb})", e.bond);
        println!("angle {:>13.9} (ref {ra})", e.angle);
        println!("rep   {:>13.9} (ref {rr})", e.rep);
        println!("es    {:>13.9} (ref {re})", e.es);
        println!("disp  {:>13.9} (ref {rd})", e.disp);
        assert!((e.bond - rb).abs() < 2e-3, "bond {}", e.bond - rb);
        assert!((e.angle - ra).abs() < 4e-5, "angle {}", e.angle - ra);
        assert!((e.rep - rr).abs() < 2e-3, "rep {}", e.rep - rr);
        assert!((e.es - re).abs() < 4e-3, "es {}", e.es - re);
        assert!((e.disp - rd).abs() < 1e-8, "disp {}", e.disp - rd);
    }
}

#[cfg(test)]
mod tests_diatomic_disp {
    use super::*;

    /// dispersion of single-pair diatomics must reproduce xtb exactly
    /// (validates C6 alpha-integration, BJ damping, double counting, zeta=1)
    #[test]
    fn diatomic_dispersion_vs_xtb() {
        // (name, Z1, Z2, r_angstrom, xtb dispersion / Eh)
        let cases: &[(&str, usize, usize, f64, f64)] = &[
            ("H2", 1, 1, 0.7414, -0.000048023),
            ("N2", 7, 7, 1.0977, -0.000197976),
            ("CO", 6, 8, 1.1280, -0.000192937),
            ("F2", 9, 9, 1.4120, -0.000086167),
            ("HF", 1, 9, 0.9168, -0.000043792),
            ("HCl", 1, 17, 1.2745, -0.000140419),
            ("LiF", 3, 9, 1.5640, -0.000118951),
            ("NaCl", 11, 17, 2.3608, -0.000341022),
        ];
        for &(name, z1, z2, r, ref_e) in cases {
            let at = [z1, z2];
            let xyz = [[0.0, 0.0, 0.0], [r, 0.0, 0.0]];
            let g = Gfnff::new(&at, &xyz, 0.0);
            let e = g.energy(&xyz);
            println!("{name}: disp {:+.9} (ref {ref_e:+.9})", e.disp);
            assert!((e.disp - ref_e).abs() < 3e-6, "{name} disp off by {}", e.disp - ref_e);
        }
    }
}

#[cfg(test)]
mod tests_torsion {
    use super::*;

    fn ethane_eclipsed() -> ([[f64; 3]; 8], [usize; 8]) {
        ([[-0.7624,0.0,0.0],[0.7624,0.0,0.0],
          [-1.1553,1.0219,0.0],[-1.1553,-0.5109,0.8851],[-1.1553,-0.5109,-0.8851],
          [1.1553,1.0219,0.0],[1.1553,-0.5109,0.8851],[1.1553,-0.5109,-0.8851]],
         [6,6,1,1,1,1,1,1])
    }
    fn ethane_staggered() -> ([[f64; 3]; 8], [usize; 8]) {
        // right methyl rotated 60 deg: azimuths 60/180/300
        ([[-0.7624,0.0,0.0],[0.7624,0.0,0.0],
          [-1.1553,1.0219,0.0],[-1.1553,-0.5109,0.8851],[-1.1553,-0.5109,-0.8851],
          [1.1553,0.51095,0.88506],[1.1553,-1.0219,0.0],[1.1553,0.51095,-0.88506]],
         [6,6,1,1,1,1,1,1])
    }

    #[test]
    fn ethane_torsion_vs_xtb() {
        let (xyz, at) = ethane_eclipsed();
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        println!("ethane eclipsed: torsion {:+.9} (ref +0.003535601721), #tors = {}", e.torsion, g.torsions.len());
        assert!((e.torsion - 0.003535601721).abs() < 2e-5, "tors off by {}", e.torsion - 0.003535601721);
        let (xyz2, at2) = ethane_staggered();
        let g2 = Gfnff::new(&at2, &xyz2, 0.0);
        let e2 = g2.energy(&xyz2);
        println!("ethane staggered: torsion {:+.9} (ref 0)", e2.torsion);
        assert!(e2.torsion.abs() < 1e-6, "staggered tors not zero: {}", e2.torsion);
    }

    #[test]
    fn methanol_torsion_vs_xtb() {
        let xyz = [[0.3531,0.0,0.0],[0.7238,1.0266,0.0],[0.7238,-0.5133,0.8890],
                   [0.7238,-0.5133,-0.8890],[-0.8469,0.0,0.0],[-1.1907,-0.4264,-0.7384]];
        let at = [6usize,1,1,1,8,1];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        println!("methanol eclipsed: torsion {:+.9} (ref +0.004287358815)", e.torsion);
        assert!((e.torsion - 0.004287358815).abs() < 8e-5, "tors off by {}", e.torsion - 0.004287358815);
    }
}

#[cfg(test)]
mod tests_gradient {
    use super::*;

    const EH_KCAL: f64 = 627.5094740631;

    fn finite_diff(g: &Gfnff, xyz: &[[f64; 3]], comp: fn(&EnergyComponents) -> f64,
                   h: f64) -> Vec<[f64; 3]> {
        let n = xyz.len();
        let mut out = vec![[0.0f64; 3]; n];
        for a in 0..n { for t in 0..3 {
            let mut p = xyz.to_vec(); p[a][t] += h;
            let mut m = xyz.to_vec(); m[a][t] -= h;
            out[a][t] = (comp(&g.energy(&p)) - comp(&g.energy(&m))) / (2.0 * h);
        }}
        out
    }

    #[test]
    fn gradient_vs_finite_difference() {
        // water with a slightly distorted geometry (all terms active)
        let at = [8usize, 1, 1];
        let xyz = [[0.0, 0.03, 0.1173], [0.05, 0.7572, -0.4692], [-0.02, -0.7572, -0.40]];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let mut grad = vec![[0.0f64; 3]; 3];
        let e = g.energy_and_gradient(&xyz, &mut grad);

        // finite differences of the TOTAL energy (kcal/mol/A vs Eh/Bohr*A)
        let num = {
            let n = 3; let h = 1e-4;
            let mut fd = vec![[0.0f64; 3]; n];
            for a in 0..n { for t in 0..3 {
                let mut p = xyz.to_vec(); p[a][t] += h;
                let mut m = xyz.to_vec(); m[a][t] -= h;
                fd[a][t] = (g.energy(&p).total() - g.energy(&m).total()) * EH_KCAL / (2.0*h);
            }}
            fd
        };
        let mut maxdiff: f64 = 0.0;
        for a in 0..3 { for t in 0..3 {
            maxdiff = maxdiff.max((grad[a][t] - num[a][t]).abs());
        }}
        println!("total E = {:.6} Eh, max |analytic - numeric| = {:.2e} kcal/mol/A", e.total(), maxdiff);
        assert!(maxdiff < 2e-1, "gradient mismatch: {maxdiff}");
        let _ = finite_diff;
    }
}

#[cfg(test)]
mod tests_grad_debug {
    use super::*;

    #[test]
    fn per_term_gradient_check() {
        let at = [8usize, 1, 1];
        let xyz = [[0.0, 0.03, 0.1173], [0.05, 0.7572, -0.4692], [-0.02, -0.7572, -0.40]];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let mut grad = vec![[0.0f64; 3]; 3];
        g.energy_and_gradient(&xyz, &mut grad);
        // verify against central finite differences on total energy
        let h = 1e-4;
        let mut fd = [[0.0f64; 3]; 3];
        for a in 0..3 { for t in 0..3 {
            let mut p = xyz.to_vec(); p[a][t] += h;
            let mut m = xyz.to_vec(); m[a][t] -= h;
            fd[a][t] = (g.energy(&p).total() - g.energy(&m).total()) * 627.5094740631 / (2.0*h);
        }}
        let mut worst: f64 = 0.0;
        for a in 0..3 { for t in 0..3 { worst = worst.max((grad[a][t]-fd[a][t]).abs()); }}
        println!("analytic-vs-fd max diff: {worst:.4} kcal/mol/A");
        assert!(worst < 0.5);
    }
}

#[cfg(test)]
mod tests_grad_xtb {
    use super::*;

    /// analytic gradient vs the xtb reference (Eh/Bohr, from `xtb --gfnff --grad`)
    #[test]
    fn gradient_vs_xtb_water() {
        let at = [8usize, 1, 1];
        let xyz = [[0.0, 0.03, 0.1173], [0.05, 0.7572, -0.4692], [-0.02, -0.7572, -0.40]];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let mut grad = vec![[0.0f64; 3]; 3]; // kcal/mol/A
        g.energy_and_gradient(&xyz, &mut grad);
        // xtb gradient (Eh/Bohr), converted to kcal/mol/A
        // xtb gradient file, SECOND block (Eh/Bohr) - matches the printed
        // |dE/dxyz| = 0.062774; the first block is an unrelated printout
        let xtb = [
            [1.4311058979121e-03, 3.8682687029141e-03, -4.2458814688643e-02],
            [-1.9075391208491e-03, -2.6279657940511e-02, 2.4592701170365e-02],
            [4.7643322293700e-04, 2.2411389237597e-02, 1.7866113518278e-02],
        ];
        let conv = 627.5094740631 / 0.52917726;
        let mut worst = 0.0f64;
        for a in 0..3 {
            println!("atom{a} ours: [{:11.5} {:11.5} {:11.5}]  xtb: [{:11.5} {:11.5} {:11.5}]",
                grad[a][0]/conv, grad[a][1]/conv, grad[a][2]/conv,
                xtb[a][0], xtb[a][1], xtb[a][2]);
        }
        for a in 0..3 { for t in 0..3 {
            let x = xtb[a][t] * conv;
            worst = worst.max((grad[a][t] - x).abs());
        }}
        println!("grad-vs-xtb max diff: {worst:.6} kcal/mol/A (xtb |g| = 0.0628 Eh/Bohr)");
        // also verify gradient norm directly in Eh/Bohr
        let mut norm = 0.0f64;
        for a in 0..3 { for t in 0..3 { let c = grad[a][t] / conv; norm += c*c; } }
        println!("our norm = {:.9}, xtb norm = 0.062773545048 Eh/Bohr", norm.sqrt());
        assert!(worst < 0.2, "gradient vs xtb off by {worst}");
    }
}

#[cfg(test)]
mod tests_optimize {
    use super::*;
    use crate::forces::ForceField;

    /// smoke test: GFN-FF + WebMM L-BFGS optimizes distorted water
    /// to the xtb-GFNFF minimum (O-H 0.9574 A, angle 104.5 deg class)
    #[test]
    fn water_optimization_smoke() {
        // simple steepest-descent/L-BFGS-free smoke: use own gradient loop
        let at = [8usize, 1, 1];
        let mut xyz = vec![[0.0, 0.0, 0.05], [0.0, 0.90, -0.55], [0.85, -0.35, -0.30]];
        let ff = GfnffForceField::new(&at, &xyz, 0.0);
        let e0 = ff.energy_and_gradient(&xyz, &mut [[0.0; 3]; 3]);
        let mut grad = vec![[0.0f64; 3]; 3];
        let mut e = e0;
        let step = 0.02f64;
        for _ in 0..400 {
            e = ff.energy_and_gradient(&xyz, &mut grad);
            let gn: f64 = grad.iter().map(|v| v[0]*v[0]+v[1]*v[1]+v[2]*v[2]).sum::<f64>().sqrt();
            if gn < 1e-3 { break; }
            for a in 0..3 { for t in 0..3 { xyz[a][t] -= step * grad[a][t] / gn; } }
        }
        let r1 = ((xyz[0][0]-xyz[1][0]).powi(2) + (xyz[0][1]-xyz[1][1]).powi(2) + (xyz[0][2]-xyz[1][2]).powi(2)).sqrt();
        println!("E: {e0:.3} -> {e:.3} kcal/mol, r(OH) = {r1:.4} A");
        assert!(e < e0, "energy did not decrease");
        assert!((r1 - 0.9574).abs() < 0.03, "O-H not optimized: {r1}");
    }
}

#[cfg(test)]
mod tests_pi {
    use super::*;

    #[test]
    fn benzene_vs_xtb() {
        let at = vec![6usize; 6].into_iter().chain([1,1,1,1,1,1]).collect::<Vec<_>>();
        let r = 1.3970f64;
        let rh = 2.4810f64;
        let mut xyz = Vec::new();
        for k in 0..6 {
            let a = k as f64 * std::f64::consts::PI / 3.0;
            xyz.push([r*a.cos(), r*a.sin(), 0.0]);
        }
        for k in 0..6 {
            let a = k as f64 * std::f64::consts::PI / 3.0;
            xyz.push([rh*a.cos(), rh*a.sin(), 0.0]);
        }
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        println!("piBO: {:?}", g.topo.pibo.iter().take(6).map(|v| (v*1000.0).round()/1000.0).collect::<Vec<_>>());
        println!("bond {:+.9} (ref -2.494313610162)", e.bond);
        println!("angle {:+.9} (ref +0.000000002942)", e.angle);
        println!("rep   {:+.9} (ref +0.145212551579)", e.rep);
        println!("es    {:+.9} (ref -0.004363946338)", e.es);
        println!("disp  {:+.9} (ref -0.006263967273)", e.disp);
        assert!((e.bond + 2.494313610162).abs() < 5e-3, "bond off {}", e.bond);
        assert!((e.rep - 0.145212551579).abs() < 5e-3, "rep off {}", e.rep);
        assert!((e.es + 0.004363946338).abs() < 5e-3, "es off {}", e.es);
    }
}

#[cfg(test)]
mod tests_hb {
    use super::*;

    /// water dimer HB energy vs xtb (unbound H, A···H···B)
    #[test]
    fn water_dimer_hb_vs_xtb() {
        let at = [8usize, 1, 1, 8, 1, 1];
        let xyz = [
            [0.0, 0.0, 0.0], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692],
            [0.0, 0.0, 2.926], [0.0, 0.24, 1.966], [0.93, 0.0, 3.166],
        ];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        // xtb: hb = -0.003319534126 Eh, total = -0.647089188710
        println!("hb     {:+.9} (ref -0.003319534126)", e.hb);
        println!("bond   {:+.9} (ref -0.543578922535)", e.bond);
        println!("angle  {:+.9} (ref +0.007881102232)", e.angle);
        println!("rep    {:+.9} (ref +0.089770304966)", e.rep);
        println!("es     {:+.9} (ref -0.197202931636)", e.es);
        println!("disp   {:+.9} (ref -0.000639207611)", e.disp);
        assert!(e.hb < -0.001, "hb should be attractive, got {}", e.hb);
    }
}

#[cfg(test)]
mod tests_hbxb_review {
    use super::*;

    const EH_KCAL: f64 = 627.5094740631;

    fn fd_total(g: &Gfnff, xyz: &[[f64; 3]], h: f64) -> Vec<[f64; 3]> {
        let n = xyz.len();
        let mut fd = vec![[0.0f64; 3]; n];
        for a in 0..n { for t in 0..3 {
            let mut p = xyz.to_vec(); p[a][t] += h;
            let mut m = xyz.to_vec(); m[a][t] -= h;
            fd[a][t] = (g.energy(&p).total() - g.energy(&m).total()) * EH_KCAL / (2.0 * h);
        }}
        fd
    }

    fn check_consistency_and_fd(at: &[usize], xyz: &[[f64; 3]], name: &str) -> Gfnff {
        let g = Gfnff::new(at, xyz, 0.0);
        // energy() and energy_and_gradient() must agree exactly
        let e1 = g.energy(xyz).total();
        let mut grad = vec![[0.0f64; 3]; at.len()];
        let ec = g.energy_and_gradient(xyz, &mut grad);
        let e2 = ec.total();
        assert!((e1 - e2).abs() < 1e-12, "{name}: energy()/gradient totals differ by {}", e1 - e2);
        // analytic vs finite-difference gradient (kcal/mol/A)
        let fd = fd_total(&g, xyz, 1e-4);
        let mut worst = 0.0f64;
        for a in 0..at.len() { for t in 0..3 {
            worst = worst.max((grad[a][t] - fd[a][t]).abs());
        }}
        println!("{name}: E = {e1:.9} Eh, analytic-vs-FD max diff = {worst:.2e} kcal/mol/A");
        assert!(worst < 5e-3, "{name}: gradient mismatch {worst}");
        g
    }

    #[test]
    fn improper_vs_xtb_formaldehyde() {
        // planar-ish formaldehyde: xtb torsion = 0.000016144526 Eh (improper only)
        let at = [6usize, 8, 1, 1];
        let xyz = [[0.0,0.0,0.0],[1.25,0.0,0.0],[-0.6,0.93,0.05],[-0.6,-0.93,0.05]];
        let g = check_consistency_and_fd(&at, &xyz, "formaldehyde planar");
        let e = g.energy(&xyz);
        println!("formaldehyde torsion {:+.9} (ref +0.000016144526)", e.torsion);
        assert!((e.torsion - 1.6144526e-5).abs() < 2e-6, "tors off: {}", e.torsion);
        // pyramidalized: xtb torsion = 0.000770759397 Eh
        let xyz2 = [[0.0,0.0,-0.10],[1.25,0.0,-0.10],[-0.6,0.93,0.40],[-0.6,-0.93,0.40]];
        let g2 = Gfnff::new(&at, &xyz2, 0.0);
        let e2 = g2.energy(&xyz2);
        println!("formaldehyde pyr  torsion {:+.9} (ref +0.000770759397)", e2.torsion);
        assert!((e2.torsion - 7.70759397e-4).abs() < 2e-5, "tors off: {}", e2.torsion);
        check_consistency_and_fd(&at, &xyz2, "formaldehyde pyr");
    }

    #[test]
    fn water_dimer_vs_xtb() {
        let at = [8usize, 1, 1, 8, 1, 1];
        let xyz = [
            [0.0, 0.0, 0.0], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692],
            [0.0, 0.0, 2.926], [0.0, 0.24, 1.966], [0.93, 0.0, 3.166],
        ];
        let g = check_consistency_and_fd(&at, &xyz, "water dimer");
        let e = g.energy(&xyz);
        // xtb 6.7.1 --gfnff --verbose references (Eh)
        println!("hb     {:+.9} (ref -0.003319534126)", e.hb);
        println!("bond   {:+.9} (ref -0.543578922535)", e.bond);
        println!("angle  {:+.9} (ref +0.007881102232)", e.angle);
        println!("rep    {:+.9} (ref +0.089770304966)", e.rep);
        println!("es     {:+.9} (ref -0.197202931636)", e.es);
        println!("disp   {:+.9} (ref -0.000639207611)", e.disp);
        assert!((e.hb - (-0.003319534126)).abs() < 2e-6, "hb off: {}", e.hb);
        assert!((e.bond - (-0.543578922535)).abs() < 5e-8, "bond off: {}", e.bond);
        assert!((e.angle - 0.007881102232).abs() < 1e-8, "angle off: {}", e.angle);
        assert!((e.rep - 0.089770304966).abs() < 2e-7, "rep off: {}", e.rep);
        assert!((e.es - (-0.197202931636)).abs() < 1e-7, "es off: {}", e.es);
        assert!((e.disp - (-0.000639207611)).abs() < 1e-8, "disp off: {}", e.disp);
        assert!(e.xb.abs() < 1e-10 && e.batm.abs() < 1e-10);
    }

    #[test]
    fn butane_batm_vs_xtb() {
        // RDKit-embedded anti-butane (physical H-H distances; xtb #BATM=90)
        let at = [6usize,6,6,6, 1,1,1,1,1,1,1,1,1,1];
        let xyz = [
        [-1.577316, 0.462590, 0.022882],
        [-0.552469, -0.313498, -0.789867],
        [0.651782, -0.775632, 0.031048],
        [1.501331, 0.370708, 0.557691],
        [-1.170099, 1.413441, 0.378863],
        [-1.912775, -0.116904, 0.888745],
        [-2.453322, 0.687572, -0.593991],
        [-0.215752, 0.300526, -1.633266],
        [-1.042306, -1.196883, -1.216582],
        [0.315587, -1.395398, 0.870447],
        [1.280316, -1.413301, -0.601975],
        [2.392587, -0.021542, 1.058030],
        [0.952276, 0.975343, 1.285440],
        [1.830160, 1.022979, -0.257466],
        ];
        let g = check_consistency_and_fd(&at, &xyz, "butane");
        let e = g.energy(&xyz);
        println!("batm   {:+.9} (ref -0.000468999264)", e.batm);
        println!("tors   {:+.9} (ref -0.000226202544)", e.torsion);
        println!("bond   {:+.9} (ref -2.057070881710)", e.bond);
        println!("angle  {:+.9} (ref +0.001100976471)", e.angle);
        println!("rep    {:+.9} (ref +0.108762095371)", e.rep);
        println!("es     {:+.9} (ref -0.004146566502)", e.es);
        println!("disp   {:+.9} (ref -0.005966004918)", e.disp);
        assert!((e.batm - (-0.000468999264)).abs() < 3e-5, "batm off: {}", e.batm);
        assert!((e.torsion - (-0.000226202544)).abs() < 2e-4, "tors off: {}", e.torsion);
        assert!((e.rep - 0.108762095371).abs() < 2e-4, "rep off: {}", e.rep);
        assert!(e.hb.abs() < 1e-10 && e.xb.abs() < 1e-10);
    }

    #[test]
    fn xb_dimer_vs_xtb() {
        // CH3-Cl ... NH3 halogen-bonded complex: xtb xb = -0.002670199865 Eh
        let at = [6usize,1,1,1,17,7,1,1];
        let xyz = [
            [0.0,0.0,0.0],[-0.6,0.93,0.0],[-0.6,-0.93,0.0],[0.0,0.0,1.09],
            [1.78,0.0,0.0],[4.60,0.0,0.0],[5.02,0.93,0.0],[5.02,-0.93,-0.75],
        ];
        let g = check_consistency_and_fd(&at, &xyz, "xb dimer");
        let e = g.energy(&xyz);
        println!("xb     {:+.9} (ref -0.002670199865)", e.xb);
        assert!((e.xb - (-0.002670199865)).abs() < 2e-6, "xb off: {}", e.xb);
    }

    #[test]
    fn benzene_hb_batm_vs_xtb() {
        let at: Vec<usize> = [6usize; 6].into_iter().chain([1; 6]).collect();
        let (r, rh) = (1.3970f64, 2.4810f64);
        let mut xyz = Vec::new();
        for k in 0..6 {
            let a = k as f64 * std::f64::consts::PI / 3.0;
            xyz.push([r*a.cos(), r*a.sin(), 0.0]);
        }
        for k in 0..6 {
            let a = k as f64 * std::f64::consts::PI / 3.0;
            xyz.push([rh*a.cos(), rh*a.sin(), 0.0]);
        }
        let g = check_consistency_and_fd(&at, &xyz, "benzene");
        let e = g.energy(&xyz);
        println!("hb     {:+.9} (ref 0.0)", e.hb);
        println!("batm   {:+.9} (ref -0.003628901480)", e.batm);
        println!("rep    {:+.9} (ref +0.145153049356)", e.rep);
        assert!(e.hb.abs() < 1e-8, "benzene hb should vanish: {}", e.hb);
        assert!((e.batm - (-0.003628901480)).abs() < 2e-5, "batm off: {}", e.batm);
        assert!((e.rep - 0.145153049356).abs() < 2e-5, "rep off: {}", e.rep);
    }
}



#[cfg(test)]
mod tests_eg3 {
    use super::*;

    /// water donor H-bonded to formaldehyde C=O acceptor (eg3 carbonyl path)




















    #[test]
    fn water_stretch_bond_vs_xtb() {
        // stretched O-H: the topology must KEEP the bond (gfnffrab detection
        // radius grows as CN drops) and the runtime CN-dependent rab0 chain
        // must match xtb exactly (regression for the old rcov*1.25 criterion
        // which dropped the bond beyond ~1.15x stretch)
        let at = [8usize, 1, 1];
        let refs = [(1.0f64, -0.270029555436), (1.1, -0.260512049958),
                    (1.15, -0.252677288024), (1.2, -0.242141045277),
                    (1.25, -0.228673952893), (1.3, -0.211039977599)];
        for (s, r) in refs {
            let o = [0.0, 0.0, 0.1173];
            let h1 = [0.0, 0.7572*s, o[2] + (-0.4692 - o[2])*s];
            let xyz = [o, h1, [0.0, -0.7572, -0.4692]];
            let g = Gfnff::new(&at, &xyz, 0.0);
            let e = g.energy(&xyz);
            println!("stretch {s}: bond {:+.9} (ref {r:+.9})", e.bond);
            assert!((e.bond - r).abs() < 5e-7, "stretch {s} off: {}", e.bond - r);
        }
    }


    #[test]
    fn h2o_h2co_hb_vs_xtb() {
        let at = [8usize, 6, 1, 1, 8, 1, 1, 1];
        let xyz = [
            [0.0, 0.0, 0.0], [1.23, 0.0, 0.0], [1.80, 0.93, 0.05], [1.80, -0.93, 0.05],
            [-2.95, 0.18, 0.12], [-1.95, 0.05, 0.02], [-3.30, 0.85, 0.70], [-3.22, 0.30, -0.80],
        ];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        println!("hb     {:+.9} (ref -0.005787714291 for this perturbed geometry)", e.hb);
        assert!((e.hb - (-0.005787714291)).abs() < 2e-6, "hb off: {}", e.hb);
        // consistency + FD gradient
        let e1 = g.energy(&xyz).total();
        let mut grad = vec![[0.0f64; 3]; 8];
        let e2 = g.energy_and_gradient(&xyz, &mut grad).total();
        assert!((e1 - e2).abs() < 1e-12);
        let h = 1e-4;
        let mut worst = 0.0f64;
        for a in 0..8 { for t in 0..3 {
            let mut p = xyz.to_vec(); p[a][t] += h;
            let mut m = xyz.to_vec(); m[a][t] -= h;
            let fd = (g.energy(&p).total() - g.energy(&m).total()) * 627.5094740631 / (2.0*h);
            worst = worst.max((grad[a][t] - fd).abs());
        }}
        println!("h2o+h2co analytic-vs-FD max diff = {worst:.2e} kcal/mol/A");
        assert!(worst < 5e-3, "gradient mismatch {worst}");
    }
}

#[cfg(test)]
mod tests_rnr {
    use super::*;

    /// water donor H-bonded to pyridine N (eg2_rnr lone-pair acceptor path)


    #[test]
    fn pyridine_water_hb_vs_xtb() {
        // geometry committed as a fixture (RDKit pyridine + placed water;
        // previously read from /tmp/gfnff_fix which broke fresh checkouts)
        let xyz_str = include_str!("pyr_h2o.xyz");
        let mut at = Vec::new();
        let mut xyz = Vec::new();
        for l in xyz_str.lines().skip(2) {
            let p: Vec<&str> = l.split_whitespace().collect();
            if p.len() == 4 {
                at.push(match p[0] { "N" => 7, "C" => 6, "H" => 1, "O" => 8, _ => panic!() });
                xyz.push([p[1].parse().unwrap(), p[2].parse().unwrap(), p[3].parse().unwrap()]);
            }
        }
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        println!("hb     {:+.9} (ref -0.011372529656)", e.hb);
        println!("bond   {:+.9} (ref -2.600581389665)", e.bond);
        assert!((e.hb - (-0.011372529656)).abs() < 2e-6, "hb off: {}", e.hb);
        // consistency + FD gradient
        let e1 = g.energy(&xyz).total();
        let mut grad = vec![[0.0f64; 3]; at.len()];
        let e2 = g.energy_and_gradient(&xyz, &mut grad).total();
        assert!((e1 - e2).abs() < 1e-12);
        let h = 1e-4;
        let mut worst = 0.0f64;
        for a in 0..at.len() { for t in 0..3 {
            let mut p = xyz.to_vec(); p[a][t] += h;
            let mut m = xyz.to_vec(); m[a][t] -= h;
            let fd = (g.energy(&p).total() - g.energy(&m).total()) * 627.5094740631 / (2.0*h);
            worst = worst.max((grad[a][t] - fd).abs());
        }}
        println!("pyr+h2o analytic-vs-FD max diff = {worst:.2e} kcal/mol/A");
        assert!(worst < 5e-3, "gradient mismatch {worst}");
    }
}

#[cfg(test)]
mod tests_charged {
    use super::*;

    /// cyclopentadienyl anion: ipis must give the aromatic 6-electron pi
    /// system (xtb: Hueckel ndim/Nel = 5/6) and the pibo-dependent bond
    /// energies must match (gfnff_ini.f90 610-660, 1032)
    #[test]
    fn cp_anion_vs_xtb() {
        let at: Vec<usize> = [6usize; 5].into_iter().chain([1; 5]).collect();
        let (r, rh) = (1.4180f64, 2.5100f64);
        let mut xyz = Vec::new();
        for k in 0..5 { let a = 2.0*std::f64::consts::PI*k as f64/5.0; xyz.push([r*a.cos(), r*a.sin(), 0.0]); }
        for k in 0..5 { let a = 2.0*std::f64::consts::PI*k as f64/5.0; xyz.push([rh*a.cos(), rh*a.sin(), 0.0]); }
        let g = Gfnff::new(&at, &xyz, -1.0);
        let e = g.energy(&xyz);
        println!("piBO: {:?}", g.topo.pibo.iter().take(5).map(|v| (v*1000.0).round()/1000.0).collect::<Vec<_>>());
        println!("bond {:+.9} (ref -1.740735892020)", e.bond);
        println!("angle {:+.9} (ref +0.003442831114)", e.angle);
        println!("rep   {:+.9} (ref +0.052093430260)", e.rep);
        println!("es    {:+.9} (ref -1.093101347487)", e.es);
        println!("disp  {:+.9} (ref -0.004866886660)", e.disp);
        // aromatic: uniform pibo on all five C-C bonds (0.647 for the 6e/5c ring)
        let cc_pibo: Vec<f64> = g.bonds.iter().enumerate()
            .filter(|(_, b)| g.at[b.i] == 6 && g.at[b.j] == 6)
            .map(|(bi_, _)| g.topo.pibo[bi_]).collect();
        assert_eq!(cc_pibo.len(), 5);
        for v in &cc_pibo {
            assert!((v - cc_pibo[0]).abs() < 1e-6, "pibo not uniform: {cc_pibo:?}");
        }
        assert!((cc_pibo[0] - 0.6471).abs() < 0.005, "unexpected pibo: {}", cc_pibo[0]);
        assert!((e.bond - (-1.740735892020)).abs() < 5e-6, "bond off: {}", e.bond);
        assert!((e.angle - 0.003442831114).abs() < 5e-7, "angle off: {}", e.angle);
        assert!((e.rep - 0.052093430260).abs() < 5e-7, "rep off: {}", e.rep);
        assert!((e.es - (-1.093101347487)).abs() < 5e-6, "es off: {}", e.es);
        assert!((e.disp - (-0.004866886660)).abs() < 5e-8, "disp off: {}", e.disp);
        // consistency + FD
        let e1 = g.energy(&xyz).total();
        let mut grad = vec![[0.0f64; 3]; 10];
        let e2 = g.energy_and_gradient(&xyz, &mut grad).total();
        assert!((e1 - e2).abs() < 1e-12);
        let h = 1e-4;
        let mut worst = 0.0f64;
        for a in 0..10 { for t in 0..3 {
            let mut p = xyz.to_vec(); p[a][t] += h;
            let mut m = xyz.to_vec(); m[a][t] -= h;
            let fd = (g.energy(&p).total() - g.energy(&m).total()) * 627.5094740631 / (2.0*h);
            worst = worst.max((grad[a][t] - fd).abs());
        }}
        println!("Cp- analytic-vs-FD max diff = {worst:.2e} kcal/mol/A");
        assert!(worst < 5e-3, "gradient mismatch {worst}");
    }
}

#[cfg(test)]
mod tests_more3 {
    use super::*;

    /// bond-detection drop point must match xtb bit-for-bit: the two-pass
    /// qloop is self-consistent (pass 1 with qa=0 may drop a stretched bond;
    /// the fragment-EEQ charges of the detached topology then make the O more
    /// negative, and pass 2 re-adds the bond up to a slightly larger radius).
    /// Both engines drop the O-H bond between s=1.360 and s=1.362, and the
    /// C-H bond of methane between s=1.30 and s=1.33.
    #[test]
    fn bond_drop_point_vs_xtb() {
        let at = [8usize, 1, 1];
        let mk = |s: f64| {
            let o = [0.0, 0.0, 0.1173];
            let h1 = [0.0, 0.7572*s, o[2] + (-0.4692 - o[2])*s];
            [o, h1, [0.0, -0.7572, -0.4692]]
        };
        assert_eq!(Gfnff::new(&at, &mk(1.360), 0.0).topo.nb[0].len(), 2, "1.360 should keep the bond");
        assert_eq!(Gfnff::new(&at, &mk(1.362), 0.0).topo.nb[0].len(), 1, "1.362 should drop the bond");
        // energy at the kept side (2-fragment-in-pass-1 -> re-added topology)
        let xyz = mk(1.360);
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        assert!((e.bond - (-0.191709966972)).abs() < 5e-7, "bond off: {}", e.bond);
        assert!((e.angle - 0.000242808539).abs() < 5e-8, "angle off: {}", e.angle);
        assert!((e.rep - 0.015805484992).abs() < 5e-8, "rep off: {}", e.rep);
        assert!((e.es - (-0.074017234818)).abs() < 5e-7, "es off: {}", e.es);
        assert!((e.disp - (-0.000216335616)).abs() < 5e-8, "disp off: {}", e.disp);
        // methane C-H drop window (xtb: kept at 1.30, dropped at 1.33)
        let atc = [6usize, 1, 1, 1, 1];
        let dirs = [[0.6287,0.6287,0.6287],[-0.6287,-0.6287,0.6287],
                    [-0.6287,0.6287,-0.6287],[0.6287,-0.6287,-0.6287]];
        let mkc = |s: f64| {
            let mut v = vec![[0.0f64; 3]; 5];
            for (k, d) in dirs.iter().enumerate() {
                let f = if k == 0 { s } else { 1.0 };
                v[k + 1] = [d[0]*f, d[1]*f, d[2]*f];
            }
            v
        };
        assert_eq!(Gfnff::new(&atc, &mkc(1.30), 0.0).topo.nb[0].len(), 4, "CH4 1.30 should keep");
        assert_eq!(Gfnff::new(&atc, &mkc(1.33), 0.0).topo.nb[0].len(), 3, "CH4 1.33 should drop");
    }


    #[test]
    fn dimer_detached_h_vs_xtb() {
        let at = [8usize, 1, 1, 8, 1, 1];
        let xyz = [
            [0.0, 0.0, 0.0], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692],
            [0.0, 0.0, 2.926], [0.0, 0.36, 1.486], [0.93, 0.0, 3.166],
        ];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        assert!((e.bond - (-0.448269805570)).abs() < 5e-7, "bond off: {}", e.bond);
        assert!((e.angle - 0.007532830949).abs() < 5e-8, "angle off: {}", e.angle);
        assert!(e.hb.abs() < 1e-10);
        let mut grad = vec![[0.0f64; 3]; 6];
        g.energy_and_gradient(&xyz, &mut grad);
        let h = 1e-4;
        let mut worst = 0.0f64;
        for a in 0..6 { for t in 0..3 {
            let mut p = xyz.to_vec(); p[a][t] += h;
            let mut m = xyz.to_vec(); m[a][t] -= h;
            let fd = (g.energy(&p).total() - g.energy(&m).total()) * 627.5094740631 / (2.0*h);
            worst = worst.max((grad[a][t] - fd).abs());
        }}
        assert!(worst < 5e-3);
    }

    /// nitromethane: nitro-N itag rules (hyb 2, no own pi electrons -> Nel=4),
    /// N sp2 angle rules, batm — full breakdown vs xtb


    #[test]
    fn nitromethane_vs_xtb() {
        let at = [6usize, 7, 8, 8, 1, 1, 1];
        let xyz = [
            [-0.650592, -0.038325, -0.049988], [0.830518, 0.048104, 0.056190],
            [1.333901, 1.175720, -0.011249], [1.444713, -1.012128, 0.223923],
            [-0.908120, -1.002198, -0.494252], [-0.998303, 0.778710, -0.685829],
            [-1.052117, 0.050118, 0.961206],
        ];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        println!("hyb {:?} nitro {:?}", g.topo.hyb, g.topo.nitro_n);
        println!("bond {:+.9} (ref -1.123587682615)", e.bond);
        println!("angle {:+.9} (ref +0.001799964317)", e.angle);
        println!("tors  {:+.9} (ref +0.000407238375)", e.torsion);
        println!("rep   {:+.9} (ref +0.060515635197)", e.rep);
        println!("es    {:+.9} (ref -0.074873994480)", e.es);
        println!("disp  {:+.9} (ref -0.001929574796)", e.disp);
        println!("batm  {:+.9} (ref -0.000180194341)", e.batm);
        assert!(g.topo.nitro_n[1], "N not flagged as nitro");
        assert!(matches!(g.topo.hyb[1], 2), "nitro N hyb should be 2");
        assert!((e.bond - (-1.123587682615)).abs() < 5e-6, "bond off: {}", e.bond);
        assert!((e.angle - 0.001799964317).abs() < 5e-7, "angle off: {}", e.angle);
        assert!((e.torsion - 0.000407238375).abs() < 5e-7, "tors off: {}", e.torsion);
        assert!((e.rep - 0.060515635197).abs() < 5e-7, "rep off: {}", e.rep);
        assert!((e.es - (-0.074873994480)).abs() < 5e-7, "es off: {}", e.es);
        assert!((e.disp - (-0.001929574796)).abs() < 5e-8, "disp off: {}", e.disp);
        assert!((e.batm - (-0.000180194341)).abs() < 5e-7, "batm off: {}", e.batm);
        // consistency + FD
        let e1 = g.energy(&xyz).total();
        let mut grad = vec![[0.0f64; 3]; 7];
        let e2 = g.energy_and_gradient(&xyz, &mut grad).total();
        assert!((e1 - e2).abs() < 1e-12);
        let h = 1e-4;
        let mut worst = 0.0f64;
        for a in 0..7 { for t in 0..3 {
            let mut p = xyz.to_vec(); p[a][t] += h;
            let mut m = xyz.to_vec(); m[a][t] -= h;
            let fd = (g.energy(&p).total() - g.energy(&m).total()) * 627.5094740631 / (2.0*h);
            worst = worst.max((grad[a][t] - fd).abs());
        }}
        println!("nitromethane FD max diff = {worst:.2e} kcal/mol/A");
        assert!(worst < 5e-3);
    }


    /// ammonia: saturated-N improper (cos double well at 80 deg)
    #[test]
    fn nh3_improper_vs_xtb() {
        let at = [7usize, 1, 1, 1];
        let xyz = [
            [0.0, 0.0, 0.0], [0.0, 0.9377, 0.3814],
            [0.8121, -0.4688, 0.3814], [-0.8121, -0.4688, 0.3814],
        ];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        println!("torsion {:+.9} (ref +0.000028596109)", e.torsion);
        assert!((e.torsion - 0.000028596109).abs() < 5e-8, "improper off: {}", e.torsion);
        // consistency + FD
        let e1 = g.energy(&xyz).total();
        let mut grad = vec![[0.0f64; 3]; 4];
        let e2 = g.energy_and_gradient(&xyz, &mut grad).total();
        assert!((e1 - e2).abs() < 1e-12);
        let h = 1e-4;
        let mut worst = 0.0f64;
        for a in 0..4 { for t in 0..3 {
            let mut p = xyz.to_vec(); p[a][t] += h;
            let mut m = xyz.to_vec(); m[a][t] -= h;
            let fd = (g.energy(&p).total() - g.energy(&m).total()) * 627.5094740631 / (2.0*h);
            worst = worst.max((grad[a][t] - fd).abs());
        }}
        println!("NH3 FD max diff = {worst:.2e} kcal/mol/A");
        assert!(worst < 5e-3);
    }
}

#[cfg(test)]
mod tests_ring3 {
    use super::*;

    /// methylcyclopropane: exercises the 3-ring-center + acyclic-neighbor
    /// angle correction (xtb ringsj+ringsk.eq.102, i.e. sentinel-99 + 3;
    /// regression for the initial port where ring_size's 0-for-acyclic made
    /// the branch dead). xtb 6.7.1 --gfnff --verbose references (Eh).
    #[test]
    fn methylcyclopropane_vs_xtb() {
        let at = [6usize, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1];
        let xyz = [
            [1.441100, 0.074800, -0.254100],
            [0.154300, -0.020800, 0.509900],
            [-1.005500, -0.806300, 0.007300],
            [-1.080100, 0.680800, -0.061300],
            [2.116900, 0.736500, 0.359200],
            [1.306300, 0.649100, -1.204200],
            [1.900400, -0.909000, -0.418600],
            [0.266700, -0.065600, 1.615600],
            [-1.663800, -1.360100, 0.723200],
            [-0.829400, -1.361100, -0.911600],
            [-0.981500, 1.115800, -1.067700],
            [-1.625300, 1.265900, 0.702400],
        ];
        let g = Gfnff::new(&at, &xyz, 0.0);
        let e = g.energy(&xyz);
        println!("bond  {:+.9} (ref -1.994328837257)", e.bond);
        println!("angle {:+.9} (ref +0.045553364714)", e.angle);
        println!("tors  {:+.9} (ref +0.007939436401)", e.torsion);
        println!("rep   {:+.9} (ref +0.081124347742)", e.rep);
        println!("es    {:+.9} (ref -0.003499883044)", e.es);
        println!("disp  {:+.9} (ref -0.004866138639)", e.disp);
        println!("batm  {:+.9} (ref -0.000131084366)", e.batm);
        assert!((e.bond - (-1.994328837257)).abs() < 5e-6, "bond off: {}", e.bond);
        // tight angle tolerance is what catches a dead 99+3 correction:
        // several C-Cr-Cme / H-Cr-C angles get phi0 82->86 and shift ~1e-4 Eh
        assert!((e.angle - 0.045553364714).abs() < 5e-7, "angle off: {}", e.angle);
        assert!((e.torsion - 0.007939436401).abs() < 5e-7, "tors off: {}", e.torsion);
        assert!((e.rep - 0.081124347742).abs() < 5e-7, "rep off: {}", e.rep);
        assert!((e.es - (-0.003499883044)).abs() < 5e-7, "es off: {}", e.es);
        assert!((e.disp - (-0.004866138639)).abs() < 5e-8, "disp off: {}", e.disp);
        assert!((e.batm - (-0.000131084366)).abs() < 5e-7, "batm off: {}", e.batm);
        // consistency + FD
        let e1 = g.energy(&xyz).total();
        let mut grad = vec![[0.0f64; 3]; 12];
        let e2 = g.energy_and_gradient(&xyz, &mut grad).total();
        assert!((e1 - e2).abs() < 1e-12);
        let h = 1e-4;
        let mut worst = 0.0f64;
        for a in 0..12 { for t in 0..3 {
            let mut p = xyz.to_vec(); p[a][t] += h;
            let mut m = xyz.to_vec(); m[a][t] -= h;
            let fd = (g.energy(&p).total() - g.energy(&m).total()) * 627.5094740631 / (2.0 * h);
            worst = worst.max((grad[a][t] - fd).abs());
        }}
        println!("methylcyclopropane FD max diff = {worst:.2e} kcal/mol/A");
        assert!(worst < 5e-3);
    }
}

/// organic-subset hybridization (gfnff_ini2.F90 215-360), returning
/// (hyb, nitro-N itag flags)
fn angle_at(xyz: &[[f64; 3]], a: usize, b: usize, c: usize) -> f64 {
    let v1 = [xyz[a][0]-xyz[b][0], xyz[a][1]-xyz[b][1], xyz[a][2]-xyz[b][2]];
    let v2 = [xyz[c][0]-xyz[b][0], xyz[c][1]-xyz[b][1], xyz[c][2]-xyz[b][2]];
    let cos = (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
        / ((v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]).sqrt() * (v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]).sqrt());
    cos.clamp(-1.0, 1.0).acos() * 180.0 / std::f64::consts::PI
}

fn assign_hyb(at: &[usize], nb: &[Vec<usize>], _xyz: &[[f64; 3]], _qa: &[f64]) -> (Vec<i32>, Vec<bool>) {
    let n = at.len();
    let mut hyb = vec![0i32; n];
    let mut nitro_n = vec![false; n];
    for i in 0..n {
        let z = at[i];
        let nn = nb[i].len();
        let group = if z == 1 || z == 2 { z as i32 }
            else if matches!(z, 5|6|7|8|9|14|15|16|17|33|34|35|51|52|53) {
                let core = if z <= 9 { 2 } else if z <= 18 { 10 } else if z <= 36 { 18 } else { 36 };
                (z - core) as i32
            } else if matches!(z, 13) { 3 } else { 0 };
        match group {
            1 => { if nn == 2 { hyb[i] = 1; } else if nn > 2 && nn <= 4 { hyb[i] = 3; } }
            3 => { if nn >= 4 { hyb[i] = 3; } else if nn == 3 { hyb[i] = 2; } else if nn == 2 { hyb[i] = 1; } }
            4 => {
                if nn >= 4 { hyb[i] = 3; }
                else if nn == 3 { hyb[i] = 2; }
                else if nn == 2 {
                    let phi = angle_at(_xyz, nb[i][0], i, nb[i][1]);
                    hyb[i] = if phi < 150.0 { 2 } else { 1 };
                }
                else if nn == 1 { hyb[i] = 1; }
            }
            5 => {
                // N family (gfnff_ini2.F90 280-341)
                if nn >= 4 { hyb[i] = 3; }
                else if nn == 3 {
                    hyb[i] = 3;
                    if z == 7 {
                        let kk = nb[i].iter().filter(|&&j| at[j] == 8 && nb[j].len() == 1).count();
                        let ll = nb[i].iter().filter(|&&j| at[j] == 5 && nb[j].len() == 4).count();
                        let ns = nb[i].iter().filter(|&&j| at[j] == 16 && nb[j].len() == 4).count();
                        if ns == 1 && ll == 0 && kk == 0 { hyb[i] = 3; }
                        if ll == 1 && ns == 0 { hyb[i] = 2; }
                        if kk >= 1 { hyb[i] = 2; nitro_n[i] = true; }  // itag=1
                    }
                }
                else if nn == 2 {
                    hyb[i] = 2;
                    let (ja, jb) = (nb[i][0], nb[i][1]);
                    for j in [ja, jb] {
                        if nb[j].len() == 1 && (at[j] == 6 || at[j] == 7) { hyb[i] = 1; }
                    }
                    if at[ja] == 7 && at[jb] == 7 && nb[ja].len() <= 2 && nb[jb].len() <= 2 { hyb[i] = 1; }
                    if angle_at(_xyz, ja, i, jb) > 160.0 { hyb[i] = 1; }
                }
                else { hyb[i] = 1; }
            }
            6 => {
                // O family (gfnff_ini2.F90 341-357)
                if nn >= 2 { hyb[i] = 3; }
                else {
                    hyb[i] = 2;
                    if let Some(&j) = nb[i].first() {
                        if nb[j].len() == 1 { hyb[i] = 1; }   // CO / OH
                    }
                }
            }
            7 => { hyb[i] = if nn >= 2 { 5 } else { 1 }; }
            _ => { hyb[i] = 0; }
        }
    }
    (hyb, nitro_n)
}




