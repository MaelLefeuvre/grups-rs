#[derive(Debug, Clone)]
/// A simple struct representing a chromosome index. This is mainly used to compute Jackknife Blocks
/// # Fields:
/// - `index` : 0-based index of the chromosome
/// - `name`  : 1-based name of the chromosome (i.e. 1-22)
/// - `length`: length of the chromosome (bp)
pub struct Chromosome {
    pub index  : usize,
    pub name   : u8,
    pub length : u32,
}

impl Chromosome {
    /// Instantiate a new chromosome
    /// # Fields:
    /// - `index` : 0-based index of the chromosome
    /// - `name`  : 1-based name of the chromosome (i.e. 1-22)
    /// - `length`: length of the chromosome (bp)
    #[must_use]
    pub fn new(index: usize, name: u8, length: u32) -> Chromosome{
        Chromosome{index, name, length}
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

#[cfg(test)]
mod tests {
    use super::*;

    fn get_two_chromosome(name_first: u8, name_second: u8) -> [Chromosome; 2] {
        let chr1 = Chromosome::new(0, name_first, 1000);
        let chr2 = Chromosome::new(1, name_second, 2000);
        [chr1, chr2]
    }

    #[test]
    fn chromosome_name_equality() {
        let test_chr = get_two_chromosome(1, 1);
        assert_eq!(test_chr[0], test_chr[1]);
    }

    #[test]
    fn chromosome_name_inequality() {
        let test_chr = get_two_chromosome(1, 22);
        assert!(test_chr[0] != test_chr[1]);
    }

    #[test]
    fn chromosome_name_lge(){
        let test_chr = get_two_chromosome(1, 1);
        assert!(test_chr[0] <= test_chr[1]);
        assert!(test_chr[0] >= test_chr[1]);

    }

    #[test]
    fn chromosome_name_lgt(){
        let test_chr = get_two_chromosome(1, 2);
        assert!(test_chr[0] < test_chr[1]);
        assert!(test_chr[1] > test_chr[0]);

    }
}