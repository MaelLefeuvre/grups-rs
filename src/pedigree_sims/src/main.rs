use std::{
    path::Path,
};
use pedigree_sims::io;


fn main() {
    println!("Hello, world!");
    let file = Path::new("tests/test-data/vcf/g1k-phase3-v5a/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz");
    //let _vcf = io::VCFReader::new(file).unwrap();
    //println!("{:?}", vcf.samples());
    //for line in vcf.lines() {
    //    println!("{:?}", line);
    //}

    let ped_path = Path::new("tests/test-data/peds/example_pedigree.txt");
    let pedigree= io::pedigree_parser(ped_path).unwrap();
    pedigree.individuals.get(&"mother".to_string()).unwrap().borrow_mut().set_id("Hello!");
    //for individual in pedigree.0 {
    //    println!("{}", individual);
    //}
    //println!("{:?}", pedigree.0.get(&String::from("cousin")).unwrap().borrow());
    //let mut father = pedigree.0.get_mut(&String::from("father")).unwrap();
    //father.borrow_mut().set_id("Hello!");
    //pedigree.0.insert(father);
    //println!("{:?}", pedigree.0.get(&String::from("father")).unwrap().borrow());
    //for comparison in pedigree.1 {
    //    println!("{}", comparison);
    //}
 
    //let mut pedigree = Pedigree::new();
    //pedigree.add_individual(&"father".to_string(), None);
    //pedigree.add_individual(&"mother".to_string(), None);
    //pedigree.add_individual(&"son".to_string(), Some((&"father".to_string(), &"mother".to_string())));
    //let father = Rc::new(RefCell::new(Individual::new("father".to_string(), None)));
    //let mother = Rc::new(RefCell::new(Individual::new("mother".to_string(), None)));
    //let son = Rc::new(RefCell::new(Individual::new("son".to_string(), None)));
    //son.borrow_mut().set_parents((&father, &mother));
    //father.borrow_mut().id = Some("Hello".to_string());
    //println!("{:#?}", son);
    //pedigree.individuals.get_mut(&"father".to_string()).unwrap().borrow_mut().set_id("Hello!");

    //pedigree.add_comparison(&"paternal".to_string(), (&"father".to_string(), &"son".to_string()));
    println!("{:#?}", pedigree);
    
}
