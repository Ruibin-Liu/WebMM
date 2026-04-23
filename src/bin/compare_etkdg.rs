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
    
    println!("{} atoms", coords.len());
    for (i, c) in coords.iter().enumerate() {
        println!("{} {:.6} {:.6} {:.6}", mol.atoms[i].symbol, c[0], c[1], c[2]);
    }
}
