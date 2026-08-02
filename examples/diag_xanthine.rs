use webmm::etkdg::{generate_initial_coords_with_config, ETKDGConfig};
use webmm::mmff::MMFFForceField;
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;
fn bad(c:&[[f64;3]],i:usize,j:usize,k:usize)->f64{
    let v1=[c[i][0]-c[j][0],c[i][1]-c[j][1],c[i][2]-c[j][2]];
    let v2=[c[k][0]-c[j][0],c[k][1]-c[j][1],c[k][2]-c[j][2]];
    let n1=(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]).sqrt();
    let n2=(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]).sqrt();
    ((v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/(n1*n2)).clamp(-1.0,1.0).acos()*180.0/std::f64::consts::PI}
fn main(){
    let mol=parse_sdf(&std::fs::read_to_string("scripts/val_set/xanthine.sdf").unwrap()).unwrap();
    let ff=MMFFForceField::new(&mol,MMFFVariant::MMFF94s);
    let cfg=ETKDGConfig{random_seed:42,..Default::default()};
    let c=generate_initial_coords_with_config(&mol,&cfg);
    let(e,_)=ff.calculate_energy_and_gradient(&c);
    println!("xanthine: energy={:.3} (was -129.6, RDKit -142.8)  O-C-C@C3={:.1}deg (was 112.5, eq 126.5)",e,bad(&c,4,3,5));
}
