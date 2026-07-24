use webmm::mmff::{MMFFForceField, MMFFVariant};
use webmm::molecule::parser::parse_sdf;

fn main() {
    let sdf = "Aniline\n     RDKit          3D\n\n 14 14  0  0  0  0  0  0  0  0999 V2000\n   -1.8551    0.3019   -0.2147 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9433    1.3121   -0.5108 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4265    1.0872   -0.3490 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9000   -0.1487    0.0976 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2537   -0.3576    0.2878 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0248   -1.1486    0.4072 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3958   -0.9291    0.2472 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9206    0.4752   -0.3382 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2957    2.2773   -0.8642 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1231    1.8892   -0.5767 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5964   -1.2716    0.5480 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9224    0.3435    0.0017 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3154   -2.1119    0.7767 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1023   -1.7188    0.4874 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  4  6  1  0\n  6  7  2  0\n  7  1  1  0\n  1  8  1  0\n  2  9  1  0\n  3 10  1  0\n  5 11  1  0\n  5 12  1  0\n  6 13  1  0\n  7 14  1  0\nM  END";
    let mol = parse_sdf(sdf).expect("parse");
    let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
    
    // RDKit MMFF94s coords (planar NH2)
    let rdkit: Vec<[f64; 3]> = vec![
        [-1.8620, 0.3170, -0.1937],
        [-0.9752, 1.3728, -0.2018],
        [0.3772, 1.1349, -0.0621],
        [0.8601, -0.1480, 0.0862],
        [2.2659, -0.3629, 0.2287],
        [-0.0179, -1.2087, 0.0954],
        [-1.3883, -0.9698, -0.0459],
        [-2.9257, 0.4377, -0.2981],
        [-1.3210, 2.3802, -0.3153],
        [1.0621, 1.9746, -0.0704],
        [2.6208, -0.8508, 1.0875],
        [2.9644, -0.0541, -0.4891],
        [0.3930, -2.1883, 0.2124],
        [-2.0534, -1.8347, -0.0338],
    ];
    
    // Our optimized coords (slightly pyramidal NH2)
    let ours: Vec<[f64; 3]> = vec![
        [-1.0728, -0.2855, 0.3027],
        [-0.1631, 0.7657, 0.4288],
        [1.2004, 0.5712, 0.1328],
        [1.6290, -0.6628, -0.2814],
        [1.0061, -1.6993, -0.3964],
        [-0.3542, -1.5046, -0.1109],
        [-2.4340, -0.5147, 0.5788],
        [-0.7337, -1.7341, -0.9433],
        [-0.4371, -0.1619, 0.9202],
        [-0.3991, 1.6413, 0.8846],
        [1.7302, 1.2980, -0.0285],
        [0.1728, -2.3116, 0.3510],
        [-0.5946, -2.3347, -0.5800],
        [-2.7920, 0.2068, -0.1796],
    ];
    
    println!("Energy at RDKit MMFF94s geometry:");
    let bd_rdkit = ff.calculate_energy_breakdown(&rdkit);
    println!("  bond={:.4} angle={:.4} sb={:.4} tor={:.4} oop={:.4} vdw={:.4} elec={:.4} total={:.4}",
        bd_rdkit.bond, bd_rdkit.angle, bd_rdkit.stretch_bend, bd_rdkit.torsion, bd_rdkit.oop, bd_rdkit.vdw, bd_rdkit.electrostatic,
        bd_rdkit.bond+bd_rdkit.angle+bd_rdkit.stretch_bend+bd_rdkit.torsion+bd_rdkit.oop+bd_rdkit.vdw+bd_rdkit.electrostatic);
    
    println!("\nEnergy at our optimized geometry:");
    let bd_ours = ff.calculate_energy_breakdown(&ours);
    println!("  bond={:.4} angle={:.4} sb={:.4} tor={:.4} oop={:.4} vdw={:.4} elec={:.4} total={:.4}",
        bd_ours.bond, bd_ours.angle, bd_ours.stretch_bend, bd_ours.torsion, bd_ours.oop, bd_ours.vdw, bd_ours.electrostatic,
        bd_ours.bond+bd_ours.angle+bd_ours.stretch_bend+bd_ours.torsion+bd_ours.oop+bd_ours.vdw+bd_ours.electrostatic);
    
    // Also check the SDF geometry energy  
    let sdf_coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
    println!("\nEnergy at SDF geometry:");
    let bd_sdf = ff.calculate_energy_breakdown(&sdf_coords);
    println!("  bond={:.4} angle={:.4} sb={:.4} tor={:.4} oop={:.4} vdw={:.4} elec={:.4} total={:.4}",
        bd_sdf.bond, bd_sdf.angle, bd_sdf.stretch_bend, bd_sdf.torsion, bd_sdf.oop, bd_sdf.vdw, bd_sdf.electrostatic,
        bd_sdf.bond+bd_sdf.angle+bd_sdf.stretch_bend+bd_sdf.torsion+bd_sdf.oop+bd_sdf.vdw+bd_sdf.electrostatic);
}
