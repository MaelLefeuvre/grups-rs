use crate::genome::Genome;
use crate::chromatid::Chromatid;

use std::{
    collections::BTreeMap,
    ops::{Deref, DerefMut},
};

#[derive(Default)]
pub struct Gamete(BTreeMap<u8, Chromatid>);

impl Deref for Gamete {
    type Target = BTreeMap<u8, Chromatid>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Gamete {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Gamete {

    pub fn add_chromatid(&mut self, chr: u8, chromatid: Chromatid) -> bool {
        self.insert(chr,  chromatid).is_none()
    }

    pub fn fertilize(&self, other: &Self) -> Genome {
        let mut genome = Genome::new();
        for ((_, chromatid1), (_, chromatid2)) in self.iter().zip(other.iter()) {
            assert_eq!(chromatid1.index, chromatid2.index);
            assert_eq!(chromatid1.name, chromatid2.name);
            assert_eq!(chromatid1.length, chromatid2.length);
                        assert_eq!(chromatid1.length, chromatid2.length);


            genome.add_chromosome(chromatid1.index, chromatid1.name, chromatid1.length);

            //let mut chromosome = Chromosome::new(chromatid1.index, chromatid1.name, chromatid1.length);
            for (allele1, allele2) in chromatid1.alleles.iter().zip(chromatid2.alleles.iter()) {
                assert_eq!(allele1.get_pos(), allele2.get_pos());
                assert_eq!(allele1.get_af(),  allele2.get_af() );

                genome.get_chr_mut(&chromatid1.name).unwrap().add_locus(allele1.get_pos(), (allele1.get_allele(), allele2.get_allele()), allele1.get_af())
            }
        }
        genome
    }
}
