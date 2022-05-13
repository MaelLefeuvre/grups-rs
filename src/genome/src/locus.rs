#[derive(Debug, Clone)]
pub struct Locus {
    pos: u32,
    alleles: [u8; 2],
    af: f64,
}

impl Locus {

    pub fn new(pos: u32, alleles: (u8, u8), af: f64) -> Self {
        Self{pos, alleles: [alleles.0, alleles.1], af}
    }

    pub fn get_pos(&self) -> u32 {
        self.pos
    }

    pub fn alleles_tuple(&self) -> (u8, u8) {
        (self.alleles[0], self.alleles[1])
    }

    pub fn get_af(&self) -> f64 {
        self.af
    }

    pub fn crossover(&mut self) {
        self.alleles.reverse();
    }

    pub fn alleles(&self) -> [u8; 2] {
        self.alleles
    }
}

impl PartialEq<Locus> for Locus {
    fn eq(&self, other: &Self) -> bool { 
        self.pos == other.pos && self.alleles == other.alleles
    }
}

impl Eq for Locus {}

impl PartialOrd for Locus {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Locus {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.pos).cmp(&(other.pos))
    }
}
