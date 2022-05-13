use rand::Rng;

use crate::{
    locus::Locus,
    geneticmap::GeneticMap,
    chromatid::Chromatid,
    alleles::{Allele, Alleles},
};

#[derive(Debug, Clone)]
/// A simple struct representing a chromosome index. This is mainly used to compute Jackknife Blocks
pub struct Chromosome {
    pub index  : usize,
    pub name   : u8,
    pub length : u32,
    loci       : Vec<Locus>,
}

impl Chromosome {
    pub fn new(index: usize, name: u8, length: u32) -> Chromosome{
        Chromosome{index, name, length, loci: Vec::new()}
    }

    pub fn snp_len(&self) -> usize {
        self.loci.len()
    }

    pub fn loci(&self) -> impl Iterator<Item = &Locus> {
        self.loci.iter()
    }

    pub fn add_locus(&mut self, pos: u32, alleles: (u8, u8), af: f64) {
        self.loci.push(Locus::new(pos, alleles, af))
    }

    pub fn is_empty(&self) -> bool {
        self.loci.is_empty()
    }

    pub fn meiosis(&self, genetic_map: &GeneticMap) -> Chromatid {
        println!("Performing meosis in chr {}", self.name);

        let mut rng = rand::thread_rng();

        let mut previous_position = 0;
        let mut currently_recombining = false;

        let mut gamete = self.loci.clone();
        
        for locus in gamete.iter_mut() {
            let mut interval_prob_recomb: f64 = 0.0;

            for recombination_range in genetic_map[&self.name].find(previous_position, locus.get_pos()) {
                let real_start = if previous_position < recombination_range.start {recombination_range.start} else {previous_position};
                let real_stop  = if locus.get_pos()   > recombination_range.stop  {recombination_range.stop } else {locus.get_pos()  };

                interval_prob_recomb += recombination_range.val.prob() * (real_stop as f64 - real_start as f64 + 1.0);
            }

            if rng.gen::<f64>() < interval_prob_recomb {
                println!("Switch!");
                currently_recombining = ! currently_recombining;
            }

            if currently_recombining {
                locus.crossover();
            }

            previous_position = locus.get_pos();
        }

        let alleles: (Alleles, Alleles) = gamete.into_iter()
            .map(|locus| {
                let (pos, alleles, af) = (locus.get_pos(), locus.alleles(), locus.get_af());
                (Allele::new(pos, alleles[0], af), Allele::new(pos, alleles[1], af))
            })
            .unzip();

        let random_chromatid = if rng.gen::<f64>() < 0.5 {alleles.0} else {alleles.1};
        Chromatid::new(self.index, self.name, self.length, random_chromatid)
    }
}

impl PartialEq for Chromosome {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}
impl Eq for Chromosome {}

impl PartialOrd for Chromosome {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Chromosome {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.name).cmp(&other.name)
    }
}

impl Extend<Locus> for Chromosome {
    fn extend<T: IntoIterator<Item=Locus>>(&mut self, alleles: T) {
        for locus in alleles {
            self.add_locus(locus.get_pos(), locus.alleles_tuple(), locus.get_af());
        }
    }
}