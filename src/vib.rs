//! Harmonic vibrational analysis: Hessian, frequencies, free-energy correction.
//!
//! Gives the iMTD-GC landscape its free-energy tier: conformer free energies
//! G_i(T) = E_i + G_corr,i(T) with the harmonic ZPE and vibrational entropy,
//! so Boltzmann populations can be reconstructed from G (not just E). For
//! conformers of the SAME molecule the rotational and translational entropy
//! cancel, so only the vibrational part is needed.
//!
//! Hessian: central finite differences of the ANALYTIC force-field gradient
//! (h = 0.005 Å), mass-weighted, diagonalized with the shared Jacobi solver.
//! All deterministic.

use crate::forces::ForceField;

/// Central-difference step for the Hessian (Å).
const FD_H: f64 = 0.005;

/// Frequency (cm⁻¹) per sqrt of a mass-weighted eigenvalue in
/// kcal·mol⁻¹·Å⁻²·amu⁻¹ — derived from fundamental constants:
/// 1 kcal/mol = 6.9477e-21 J per molecule, 1 amu = 1.66054e-27 kg,
/// 1 Å = 1e-10 m, ν̃ = ω/(2πc).
const WAVENUMBER_PER_SQRT: f64 = 108.6;

/// Signed harmonic frequencies (cm⁻¹) at a (near-)minimum geometry:
/// negative values mark imaginary modes. Translation/rotation leave ~0
/// eigenvalues; callers filter them (see [`harmonic_correction`]).
pub fn frequencies(ff: &dyn ForceField, coords: &[[f64; 3]], masses: &[f64]) -> Vec<f64> {
    let n = coords.len();
    let dim = 3 * n;
    // mass-weighted Hessian rows by central differences of the gradient
    let mut h = vec![0.0f64; dim * dim];
    let mut grad_plus = vec![[0.0f64; 3]; n];
    let mut grad_minus = vec![[0.0f64; 3]; n];
    for c in 0..dim {
        let atom = c / 3;
        let xyz = c % 3;
        let mut p = coords.to_vec();
        p[atom][xyz] += FD_H;
        ff.energy_and_gradient(&p, &mut grad_plus);
        p[atom][xyz] -= 2.0 * FD_H;
        ff.energy_and_gradient(&p, &mut grad_minus);
        for r in 0..dim {
            let ra = r / 3;
            let rx = r % 3;
            // M^-1/2 H M^-1/2: divide both rows and columns by sqrt(m)
            let w = 1.0 / (masses[atom] * masses[ra]).sqrt();
            h[r * dim + c] = w * (grad_plus[ra][rx] - grad_minus[ra][rx]) / (2.0 * FD_H);
        }
    }
    // symmetrize (kills FD asymmetry noise)
    for r in 0..dim {
        for c in (r + 1)..dim {
            let m = 0.5 * (h[r * dim + c] + h[c * dim + r]);
            h[r * dim + c] = m;
            h[c * dim + r] = m;
        }
    }
    let mut eig = h;
    let mut _vecs = vec![0.0f64; dim * dim];
    crate::optimizer::jacobi::sym_jacobi(&mut eig, &mut _vecs, dim);
    let mut freqs: Vec<f64> = eig
        .iter()
        .step_by(dim + 1)
        .map(|&l| WAVENUMBER_PER_SQRT * l.abs().sqrt() * if l < 0.0 { -1.0 } else { 1.0 })
        .collect();
    freqs.sort_by(|a, b| a.total_cmp(b));
    freqs
}

/// Modes below this wavenumber are treated as translation/rotation residue
/// and excluded from the harmonic sum (soft torsions are well above).
pub const TR_ROT_CUTOFF_CM: f64 = 25.0;

/// Harmonic free-energy correction at temperature `t` (K), in kcal/mol:
/// Σ over real modes of [½·kB·θ + kB·T·ln(1 − exp(−θ/T))] with θ = (hc/kB)·ν̃;
/// imaginary modes are counted (`n_imag`) and included with |ν̃| — the
/// conformer is then flagged as not a true minimum by the caller.
/// Returns (g_corr, zpe, n_imag).
pub fn harmonic_correction(freqs: &[f64], t: f64) -> (f64, f64, usize) {
    let kb = crate::md::KB; // kcal/(mol·K)
    let t = if t > 1.0 { t } else { 298.15 };
    let hc_over_kb = 1.43877; // K·cm
    let mut g_corr = 0.0;
    let mut zpe = 0.0;
    let mut n_imag = 0;
    for &nu in freqs {
        if nu.abs() < TR_ROT_CUTOFF_CM {
            continue; // translation/rotation residue
        }
        if nu < 0.0 {
            n_imag += 1;
        }
        let theta = hc_over_kb * nu.abs(); // K
        let zpe_m = 0.5 * kb * theta;
        zpe += zpe_m;
        let ln_term = (1.0 - (-theta / t).exp()).ln(); // ln(1 − e^{−θ/T})
        g_corr += zpe_m + kb * t * ln_term;
    }
    (g_corr, zpe, n_imag)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::md::atomic_mass;
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;

    const METHANE_SDF: &str = "Methane\n     RDKit          3D\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6294    0.6294    0.6294 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6294   -0.6294    0.6294 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6294    0.6294   -0.6294 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6294   -0.6294   -0.6294 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\n  1  4  1  0\n  1  5  1  0\nM  END\n";

    const WATER_SDF: &str = "Water\n     RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.1173 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7572    0.0000   -0.4694 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7572    0.0000   -0.4694 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\nM  END\n";

    fn quenched(
        sdf: &str,
    ) -> (
        crate::molecule::Molecule,
        Vec<[f64; 3]>,
        Vec<f64>,
        MMFFForceField,
    ) {
        let mol = parse_sdf(sdf).expect("parse");
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let start: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let r = crate::optimizer::optimize(
            &ff,
            &start,
            &crate::ConvergenceOptions {
                max_iterations: 1000,
                ..Default::default()
            },
        );
        let masses: Vec<f64> = mol
            .atoms
            .iter()
            .map(|a| atomic_mass(a.atomic_number))
            .collect();
        (mol, r.optimized_coords, masses, ff)
    }

    #[test]
    fn methane_frequencies_in_physical_bands() {
        let (_mol, coords, masses, ff) = quenched(METHANE_SDF);
        let freqs = frequencies(&ff, &coords, &masses);
        assert_eq!(freqs.len(), 15);
        // 6 translation/rotation residue modes near zero
        let near_zero = freqs
            .iter()
            .filter(|&&f| f.abs() < TR_ROT_CUTOFF_CM)
            .count();
        assert!(near_zero >= 5, "tr/rot residue near zero, got {freqs:?}");
        // highest C–H stretch in the 2800–3300 cm⁻¹ band (validates the
        // kcal/Å/amu → cm⁻¹ conversion constant end to end)
        let hi = *freqs.last().expect("nonempty");
        assert!(
            (2800.0..3300.0).contains(&hi),
            "methane top MMFF stretch {hi} cm⁻¹"
        );
        // all real above the tr/rot residue (tiny ±near-zero eigenvalues
        // from translation/rotation are expected)
        assert!(
            freqs.iter().all(|f| *f >= -TR_ROT_CUTOFF_CM),
            "no significant imaginary modes, got {freqs:?}"
        );
        // determinism
        let again = frequencies(&ff, &coords, &masses);
        assert_eq!(freqs, again);
    }

    #[test]
    fn water_frequencies_in_physical_bands() {
        let (_mol, coords, masses, ff) = quenched(WATER_SDF);
        let freqs = frequencies(&ff, &coords, &masses);
        let real: Vec<f64> = freqs
            .iter()
            .copied()
            .filter(|f| f.abs() >= TR_ROT_CUTOFF_CM)
            .collect();
        assert_eq!(
            real.len(),
            3,
            "9 DOF − 6 tr/rot = 3 vibrations, got {freqs:?}"
        );
        // bend < stretches; stretches in a sane MMFF band
        assert!(real[0] < real[1] && real[1] < real[2]);
        assert!(real[0] > 1200.0 && real[0] < 2000.0, "bend {}", real[0]);
        assert!(real[2] > 3300.0 && real[2] < 4200.0, "stretch {}", real[2]);
    }

    #[test]
    fn harmonic_correction_limits() {
        // single mode 1000 cm⁻¹
        let f = [1000.0];
        let (g298, zpe, n_imag) = harmonic_correction(&f, 298.15);
        assert_eq!(n_imag, 0);
        let zpe_expect = 0.5 * crate::md::KB * 1.43877 * 1000.0;
        assert!(
            (zpe - zpe_expect).abs() < 1e-12,
            "zpe {zpe} vs {zpe_expect}"
        );
        // T→0: G_corr → ZPE (entropy term vanishes)
        let (g0, _, _) = harmonic_correction(&f, 1.5);
        assert!((g0 - zpe).abs() < 1e-9, "T→0 gives {g0}, zpe {zpe}");
        // high T: entropy dominates, G_corr < ZPE eventually
        let (ghot, _, _) = harmonic_correction(&f, 5000.0);
        assert!(ghot < g298, "entropy grows with T");
        // imaginary mode counted, included with |ν|
        let (gi, _, ni) = harmonic_correction(&[-1000.0, 1000.0], 298.15);
        assert_eq!(ni, 1);
        assert!((gi - 2.0 * g298).abs() < 1e-12);
        // tr/rot residue excluded
        let (gn, _, nn) = harmonic_correction(&[0.0, 5.0, 24.9, 1000.0], 298.15);
        assert_eq!(nn, 0);
        assert!((gn - g298).abs() < 1e-12);
    }

    #[test]
    fn ethanol_thermo_is_physical() {
        let sdf = "ethanol\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n    0.8867    0.1097    0.1398 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4649   -0.1632   -0.4871 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4908    0.4341    0.2924 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0600    1.1872    0.2278 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6897   -0.3301   -0.4581 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9343   -0.3042    1.1524 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5119    0.2599   -1.4947 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6562   -1.2384   -0.5550 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4469    0.0448    1.1825 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  1  4  1  0\n  1  5  1  0\n  1  6  1  0\n  2  7  1  0\n  2  8  1  0\n  3  9  1  0\nM  END\n";
        let (_mol, coords, masses, ff) = quenched(sdf);
        let freqs = frequencies(&ff, &coords, &masses);
        let real: Vec<&f64> = freqs
            .iter()
            .filter(|f| f.abs() >= TR_ROT_CUTOFF_CM)
            .collect();
        assert_eq!(real.len(), 21, "27 DOF − 6 tr/rot, got {freqs:?}");
        assert!(
            freqs.iter().all(|f| *f >= -TR_ROT_CUTOFF_CM),
            "minimum is real above residue, got {freqs:?}"
        );
        let (g_corr, zpe, n_imag) = harmonic_correction(&freqs, 298.15);
        assert_eq!(n_imag, 0);
        assert!(
            (42.0..=62.0).contains(&zpe),
            "ethanol MMFF ZPE ~48-52 kcal/mol, got {zpe}"
        );
        // the entropy term kT·ln(1−e^{−θ/T}) is always negative: G_corr < ZPE
        assert!(
            g_corr < zpe && g_corr > zpe - 8.0,
            "G_corr below ZPE by the entropy term: {g_corr} vs {zpe}"
        );
    }
}
