use std::io::{self, Read};

fn main() {
    let mut sdf = String::new();
    io::stdin().read_to_string(&mut sdf).expect("Failed to read stdin");
    
    let mol = match webmm::molecule::parser::parse_sdf(&sdf) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("Parse error: {}", e);
            std::process::exit(1);
        }
    };
    
    let coords = webmm::etkdg::generate_initial_coords(&mol);
    let ff = webmm::mmff::MMFFForceField::new(&mol, webmm::MMFFVariant::MMFF94s);
    let energy = ff.calculate_energy(&coords);
    let (_, grad) = ff.calculate_energy_and_gradient(&coords);
    
    println!("energy: {}", energy);
    println!("atoms: {}", coords.len());
    for i in 0..coords.len() {
        let g_norm = (grad[i][0]*grad[i][0] + grad[i][1]*grad[i][1] + grad[i][2]*grad[i][2]).sqrt();
        println!("{} {:.6} {:.6} {:.6} | grad: {:.6} {:.6} {:.6} | norm: {:.6}", 
            mol.atoms[i].symbol, coords[i][0], coords[i][1], coords[i][2],
            grad[i][0], grad[i][1], grad[i][2], g_norm);
    }
}
