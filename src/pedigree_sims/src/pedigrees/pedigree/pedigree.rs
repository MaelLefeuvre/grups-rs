use std::{collections::{BTreeMap},  error::Error};

use super::{
    Contaminant,
    Individual,
    PedigreeParams,
    PedComparisons,
    PedComparison
};

use crate::io::{
    genotype_reader::GenotypeReader,
    vcf::reader::VCFPanelReader,
};

use genome::Genome;


use rand::{Rng, prelude::ThreadRng};




use std::{

    cell::{RefCell, RefMut},
    rc::Rc,

};



#[derive(Debug, Clone)]
pub struct Pedigree {
    pub individuals: BTreeMap<String, Rc<RefCell<Individual>>>,
    pub comparisons: PedComparisons,
    params: Option<PedigreeParams>,
    pop : Option<String>,
}

impl Pedigree {
    pub fn new() -> Pedigree {
        Pedigree { individuals: BTreeMap::new(), comparisons: PedComparisons::new(), params: None, pop: None}
    }

    pub fn compare_alleles(&mut self, cont_af: [f64; 2]) -> Result<(), Box<dyn Error>> {
        // --------------------- Compare genomes.
        let contam_rate = self.get_params()?.contam_rate;
        let seq_error_rate = self.get_params()?.seq_error_rate;
        for comparison in &mut self.comparisons.iter_mut() {
            comparison.compare_alleles(contam_rate, cont_af, seq_error_rate)?;
        }
        Ok(())
    }

    pub fn update_founder_alleles(&mut self, reader: &dyn GenotypeReader, rng: &mut ThreadRng) -> Result<(), Box<dyn Error>> {
        let af_downsampling_rate = self.get_params()?.af_downsampling_rate;
        for mut founder in self.founders_mut() {
            founder.alleles = match rng.gen::<f64>() < af_downsampling_rate { 
                false => reader.get_alleles(founder.get_tag().unwrap()),
                true  => Some([0, 0]),
            };
        }
        Ok(())
    }

    pub fn compute_offspring_alleles(&mut self, interval_prob_recomb: f64, pedigree_index: usize) -> Result<(), Box<dyn Error>> {
        for mut offspring in self.offsprings_mut() {
            offspring.assign_alleles(interval_prob_recomb, pedigree_index)?;
        }
        Ok(())
    }



    pub fn set_params(&mut self, snp_downsampling_rate: f64, af_downsampling_rate: f64, seq_error_rate: [f64; 2], contam_rate: [f64; 2]) {
        //println!("error_rate: {seq_error_rate} | contam_rate: {contam_rate}");
        self.params = Some(
            PedigreeParams::new(snp_downsampling_rate, af_downsampling_rate, seq_error_rate, contam_rate)
        );
    }

    pub fn get_params(&self) -> Result<&PedigreeParams, &str> {
        self.params.as_ref().ok_or("Attempt to access empty pedigree params.")
    }

    pub fn assign_offspring_strands(&mut self) -> Result<(), String> {
        for mut offspring in self.offsprings_mut() {
            offspring.assign_strands()?;
        }
        Ok(())
    }

    pub fn add_individual(&mut self, label: &str, parents: Option<(&String, &String)>, genome: Genome) -> std::io::Result<()>{
        use std::io::ErrorKind::InvalidInput;
        let parents = match parents {
            None => None,
            Some((parent1, parent2)) => {
                let parent1 = self.individuals.get(parent1).ok_or(InvalidInput)?;
                let parent2 = self.individuals.get(parent2).ok_or(InvalidInput)?;
                Some([parent1, parent2])
            },
        };
        let ind = Rc::new(RefCell::new(Individual::new(&label, parents, genome)));
        self.individuals.insert(label.to_owned(), ind);
        Ok(())
    }

    pub fn add_comparison(&mut self, label: &str, pair: (&String, &String)) -> std::io::Result<()> {
        use std::io::ErrorKind::InvalidInput;
        let pair0 = self.individuals.get(pair.0).ok_or(InvalidInput)?;
        let pair1 = self.individuals.get(pair.1).ok_or(InvalidInput)?;
        self.comparisons.push(PedComparison::new(label.to_owned(), (pair0, pair1), pair0==pair1));
        Ok(())
    }

    pub fn set_relationship(&mut self, ind: &String, parents: (&String, &String)) -> std::io::Result<()>{
        use std::io::ErrorKind::InvalidInput;
        let parent0 = &self.individuals.get(parents.0).ok_or(InvalidInput)?.clone();
        let parent1 = &self.individuals.get(parents.1).ok_or(InvalidInput)?.clone();

        self.individuals.get_mut(ind)
            .ok_or(InvalidInput)?
            .borrow_mut()
            .set_parents([parent0, parent1]);
            Ok(())

    }

    #[cfg(test)]
    pub fn get_mutind(&mut self, label: &String) -> Result<std::cell::RefMut<Individual>, std::io::Error>{
        use std::io::ErrorKind::InvalidInput;
        Ok(self.individuals.get_mut(label)
            .ok_or(InvalidInput)?
            .borrow_mut())
    }

    //#[cfg(test)]
    //pub fn get_ind(&mut self, label: &String) -> Result<std::cell::Ref<Individual>, std::io::Error>{
    //    use std::io::ErrorKind::InvalidInput;
    //    Ok(RefCell::borrow(self.individuals.get(label).ok_or(InvalidInput)?))
    //}

    /// TODO: founders_mut and offsprings_mut have the same code, except the predicate:
    ///  - founters_mut    : ind.is_founder()
    ///  - offsprings_mut : !ind.is_founder()
    pub fn founders_mut(&mut self) -> Vec<RefMut<Individual>> {
        self.individuals
            .values_mut()
            .filter(|ind| RefCell::borrow(ind).is_founder())
            .map(|x| x.borrow_mut())
            .collect()
    }

    pub fn offsprings_mut(&mut self) -> Vec<RefMut<Individual>> {
        self.individuals
            .values_mut()
            .filter(|ind| !RefCell::borrow(ind).is_founder())
            .map(|x| x.borrow_mut())
            .collect()
    }


    pub fn set_tags(&mut self, panel: &VCFPanelReader, pop: &String, contaminants: Option<&Contaminant>) {
        self.pop = Some(pop.to_owned());

        let contam_tags = contaminants.map(|cont| cont.as_flat_list());

        for mut founder in self.founders_mut() {
            founder.set_tag(panel.random_sample(pop, contam_tags.as_ref()).unwrap().clone());
        }
    }

    pub fn clear_alleles(&mut self){
        for ind in self.individuals.values_mut(){
            ind.borrow_mut().clear_alleles()
        }
    }

}

impl Default for Pedigree {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_pedigree() -> std::io::Result<Pedigree> {
        let mut pedigree = Pedigree::new();
        pedigree.add_individual("father", None, Genome::default())?;
        pedigree.add_individual("mother", None, Genome::default())?;
        pedigree.add_individual("offspr", Some((&"father".to_string(), &"mother".to_string())), Genome::default())?;

        let mut father = pedigree.get_mutind(&"father".to_string()).expect("Cannot extract father");
        father.set_alleles([0, 1]);
        drop(father);

        let mut mother = pedigree.get_mutind(&"mother".to_string()).expect("Cannot extract mother");
        mother.set_alleles([1, 0]);
        drop(mother);


        Ok(pedigree)
    }

    #[test]
    #[should_panic]
    fn meiosis_assign_alleles_empty_strands(){
        let mut pedigree = test_pedigree().expect("Cannot generate test pedigree");
        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.assign_alleles(0.0, 0).expect("Failed to assign alleles");
    }

    #[test]
    fn meiosis_assign_alleles_filled_strands(){
        let mut pedigree = test_pedigree().expect("Cannot generate test pedigree");
        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands= Some([0,0]);

        let output = offspr.assign_alleles(0.0, 0).expect("Failed to assign alleles");
        assert_eq!(output, true);

        let output = offspr.assign_alleles(0.0, 0).expect("Failed to assign alleles");
        assert_eq!(output, false);
    }

    #[test]
    fn meiosis_check_strands_00() {
        let mut pedigree = test_pedigree().expect("Cannot generate test pedigree");

        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands= Some([0,0]);
        offspr.assign_alleles(0.0, 0).expect("Failed to assign alleles");
        assert_eq!(offspr.alleles, Some([0, 1]))
    }

    #[test]
    fn meiosis_check_strands_01() {
        let mut pedigree = test_pedigree().expect("Cannot generate test pedigree");
        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands= Some([1,1]);

        offspr.assign_alleles(0.0, 0).expect("Failed to assign alleles");
        assert_eq!(offspr.alleles, Some([1, 0]));
        assert_eq!(offspr.currently_recombining, [false, false]);

    }

    #[test]
    fn meiosis_check_recombination() {
        let mut pedigree = test_pedigree().expect("Cannot generate test pedigree");
        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands= Some([0,1]);
        offspr.assign_alleles(1.0, 0).expect("Failed to assign alleles");
        assert_eq!(offspr.alleles, Some([1, 1]));
        assert_eq!(offspr.currently_recombining, [true, true]);

    }
}