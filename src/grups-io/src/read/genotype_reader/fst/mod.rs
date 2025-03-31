use std::{fs::File, path::{Path, PathBuf}, io::Read};

use crate::{
    parse,
    read::SampleTag,
    read::genotype_reader::{GenotypeReader, GenotypeReaderError},
};


use genome::coordinate::Coordinate;
use located_error::prelude::*;

use memmap2::Mmap;
use ahash::AHashMap;
use log::{info, debug};
use fst::{Set, IntoStreamer, automaton::{Automaton, Str, StartsWith}, Streamer};

mod error;
use error::FSTReaderError;

pub const FST_EXT: &str = "fst";
pub const FRQ_EXT: &str = "fst.frq";


/// Generic trait to either obtain a Set<Mmap> or a Set<Vec<u8>>.
pub trait SetRead<T: AsRef<[u8]>> {
    fn get_fst_memory(path: &str) -> Result<Set<T>>;
}

/// Implementation for the Mmap memory model of FSTReader
impl SetRead<Mmap> for Mmap {
    fn get_fst_memory(path: &str) -> Result<Set<Mmap>> {
        use FSTReaderError::{OpenFST, LoadFST, BuildFST};
        let loc_msg = || format!("While attempting to set the memory of {path}");
        
        // ---- Instantiate Memory Map
        info!("Loading Memory-mapped file: {path}");
        let file_handle = File::open(path).map_err(OpenFST).with_loc(loc_msg)?;
        let mmap = unsafe { Mmap::map(&file_handle).map_err(LoadFST).with_loc(loc_msg)? };
        Set::new(mmap).map_err(BuildFST).with_loc(loc_msg)
    }
}

/// Implementation for the in-RAM memory model of FSTReader
impl SetRead<Vec<u8>> for Vec<u8> {
    fn get_fst_memory(path: &str) -> Result<Set<Vec<u8>>> {
        let loc_msg = || format!("While attempting to set the memory of {path}");
        use FSTReaderError::{OpenFST, LoadFST, BuildFST};
        info!("Loading in memory : {path}");
        // ---- Get file size
        let size = std::fs::metadata(path).map_err(OpenFST).with_loc(loc_msg)?.len();
        debug!("File size (bytes): {size}");
        // ---- Open FST file in ro mode.
        let mut file_handle = File::open(path).map_err(OpenFST).with_loc(loc_msg)?;

        // ---- Load FST set into memory.
        let mut bytes = Vec::with_capacity(size as usize);
        Read::read_to_end(&mut file_handle, &mut bytes).map_err(LoadFST).with_loc(loc_msg)?;
        Set::new(bytes).map_err(BuildFST).with_loc(loc_msg)
    }
}


/// GenotypeReader from a pair of genotypes (`.fst`) and frequencies (`.fst.frq`) FST-Sets
/// #### Implemented Traits: GenotypeReader
pub struct FSTReader<T: SetRead<T> + AsRef<[u8]>> {
    genotypes_set: Set<T>,
    frequency_set: Set<T>,
    genotypes    : AHashMap<u128, [u8; 2]>,
    frequencies  : AHashMap<String, f32>,
}

impl<T: AsRef<[u8]> + SetRead<T>> GenotypeReader for FSTReader<T> {
    // Return the alleles for a given SampleTag. Search is performed using `sample_tag.id()`;
    fn get_alleles(&self, sample_tag: &SampleTag ) -> Result<[u8; 2]> {
        use GenotypeReaderError::MissingAlleles;
        let loc_msg = || format!("While retrieving alleles of {}", sample_tag.id());

        match self.genotypes.get(&sample_tag.hashed_id()) {
            Some(alleles) => Ok(alleles.map(|all| all - 48)),
            None          => Err(MissingAlleles).with_loc(loc_msg)
        }
    }

    // Return the alleles frequencies for a given population id.
    fn get_pop_allele_frequency(&self, pop: &str) -> Result<f32> {
        use GenotypeReaderError::MissingFreq;
        match self.frequencies.get(pop) {
            Some(freq) => Ok(*freq),
            None       => Err(MissingFreq(pop.to_string())).loc("while parsing coordinate")
        }
    }

    fn fetch_input_files(input_dir: &Path) -> Result<Vec<PathBuf>> {
        let mut fsts = parse::fetch_input_files(input_dir, &[FST_EXT]).loc("While searching for candidate .fst files.")?;
        fsts.dedup();
        debug!("Found the following fst file candidates: {:#?}", fsts);
        // Ensure there's a matching 'fst.frq' file for each found '.fst' file
        for fst in fsts.iter() {
            Self::match_input_frq(fst).with_loc(|| format!("While fetching for .fst files in {}", input_dir.display()))?;
        }
        Ok(fsts)
    }

}

impl<T: AsRef<[u8]> + SetRead<T>>FSTReader<T> {
    /// Instantiate a new reader from a `.fst` file.
    /// The companion `.fst.frq` file must be located within the same directory and display the same file-stem.
    /// # Arguments:
    /// - `path`: path leading to the targeted `.fst` file.
    pub fn new(path: &str) -> Result<FSTReader<T>> {
        let loc_msg = || format!("While attempting to create FSTReader from {path}");
        let genotypes_set = <T as SetRead<T>>::get_fst_memory(path).with_loc(loc_msg)?;
        let frequency_set = <T as SetRead<T>>::get_fst_memory(&format!("{path}.frq")).with_loc(loc_msg)?; //@TODO: const FRQ_EXT should be used in place.
        Ok(Self{genotypes_set, frequency_set, genotypes: AHashMap::new(), frequencies: AHashMap::new()})
    }

    /// Public wrapper for `find_chromosome()`. returns a sorted, unduplicated list of chromosomes contained within the set.
    pub fn find_chromosomes(&self) -> Result<Vec<u8>> {
        let root = self.genotypes_set.as_fst().root();
        let mut chromosomes = Vec::from_iter( root.transitions().map(|transition| transition.inp) );
        chromosomes.sort_unstable();
        chromosomes.dedup();
        Ok(chromosomes)
    }

    /// @TODO: This should be a macro ?
    #[inline]
    fn format_coordinate_pattern(coord_bytes: &[u8; 5]) -> StartsWith<Str> {
        let regex = unsafe { std::str::from_utf8_unchecked(coord_bytes) };
        Str::new(regex).starts_with()
    }
    
    /// Populate the `genotypes` field using genome coordinates.
    /// `self.genotypes` is thus a HashMap, with keys==<sample-id>, and values==<alleles>
    #[inline]
    pub fn search_coordinate_genotypes(&mut self, coordinate: &Coordinate) {
        // ---- Create a new fst-matcher, matching any entry starting with our chromosome coordinates.
        // ---- genotype-fst index fields are '{chr} {pos} {id} {alleles}'
        let coord_bytes: [u8; 5] = (*coordinate).into();
        let matcher = Self::format_coordinate_pattern(&coord_bytes);

        // ---- Search through the set and format each match within `self.genotypes` (key=<sample-id>, val=<alleles>)
        let mut stream = self.genotypes_set.search(&matcher).into_stream();
        while let Some(key) = stream.next() {
            // ---- Retrieve alleles and sample id.
            let key = &mut key[5..].to_vec();                // Skip coordinate bytes.
            let (id, alleles) = key.split_at(key.len() - 2); // alleles are the last two u8.

            // ---- Hash sample id. to u128
            let panic_msg        = |_| panic!("Invalid allele length. expected 2, got {}", alleles.len());
            let hashed_tag       = SampleTag::hash_id_u128(id);
            let alleles: [u8; 2] = alleles.try_into().unwrap_or_else(panic_msg);

            self.genotypes.insert(hashed_tag, alleles);
        }
    }

    /// Populate the `frequencies` field using genome coordinates.
    /// `self.frequencies` is thus a HashMap, with keys==<sample-id>, and values==<allele-frequency>
    /// # Arguments:
    /// - `chr`: chromosome name
    /// - `pos`: 0-based chromosome position
    #[inline]
    pub fn search_coordinate_frequencies(&mut self, coordinate: &Coordinate) {
        // ---- Create a new fst-matcher, matching any entry starting with our chromosome coordinates.
        //----- genotype-fst index fields are '{chr(u8)}{pos(u32_be)}{pop(chars)}{freq(f32_be)}'
        let coord_bytes: [u8; 5] = (*coordinate).into();
        let matcher = Self::format_coordinate_pattern(&coord_bytes);

        // ---- Search through the set and format each match within `self.frequencies` (key=<sample-id>, val=<allele-frequency>)
        let mut stream = self.frequency_set.search(&matcher).into_stream();
        while let Some(key) = stream.next() {
            // ---- Retrieve population tag and frequency.
            let key              = &mut key[5..].to_vec(); // Go from the end of the coordinate to the end.
            let (tag, freq_be)   = key.split_at(key.len()-4); 

            // ---- Parse population tag
            let pop              = unsafe {String::from_utf8_unchecked(tag.to_vec())};

            // ---- Convert frequency to f32
            let panic_msg        = |_| panic!("Invalid freq bytes length. expected 4, got {}", freq_be.len());
            let freq_be: [u8; 4] = freq_be.try_into().unwrap_or_else(panic_msg);
            let freq: f32        = f32::from_be_bytes(freq_be);

            self.frequencies.insert(pop, freq);
        }
    }


    /// Clear `self.genotypes` and `self.frequencies` buffer.
    /// Generally called before skipping to the next coordinate.
    pub fn clear_buffers(&mut self) {
        self.genotypes.clear();
        self.frequencies.clear();
    }

    /// Check if the `self.genotypes` field has been filled.
    pub fn has_genotypes(&self) -> bool {
       ! self.genotypes.is_empty()
    }

    /// Search for a companion `.fst.frq` file next to the provided `.fst`. Returns an error if there are none.
    /// # Arguments: 
    /// - `fst_file`: path leading to a `.fst` file.
    /// # Errors:
    /// - returns an error if there are no `.fst.frq` companion file for the provided `fst_file`
    fn match_input_frq(fst: &Path) -> Result<()> {
        let frq = fst.with_extension(FRQ_EXT);
        if ! frq.exists() {
            return Err(FSTReaderError::MatchFrqFile{path: fst.to_path_buf()})
                .loc("Note that .fst files must have a matching .fst.frq file within the same location.")
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use genome::coordinate::ChrIdx;
    
    fn fst_line(chr: u8, coord: u32, id: &str, alleles: [u8; 2]) -> Vec<u8> {
        // ---- {chr(u8)}{pos(u32_be)}{sample_id(chars)}{allele(u8)}{allele(u8)}   
        let mut line = vec![ChrIdx(chr).into()];
        line.extend_from_slice(&coord.to_be_bytes());
        line.extend_from_slice(id.as_bytes());
        line.extend_from_slice(&alleles.iter().map(|a| a + 48).collect::<Vec<u8>>());
        line
    }

    fn frq_line(chr: u8, coord: u32, pop: &str, freq: f32) -> Vec<u8> {
        // ---- {chr(u8)}{pos(u32_be)}{pop(chars)}{freq(f32_be)}   
        let mut line = vec![ChrIdx(chr).into()];
        line.extend_from_slice(&coord.to_be_bytes());
        line.extend_from_slice(pop.as_bytes());
        line.extend_from_slice(&freq.to_be_bytes());
        line
    }

    fn fake_fst() -> FSTReader<Vec<u8>> {
        let genotypes_set = Set::from_iter(vec![
            fst_line(1, 50000, "HG00096", [0, 0]),
            fst_line(1, 50000, "HG00097", [0, 1]),
            fst_line(1, 50000, "HG00098", [1, 0]),
            fst_line(1, 50000, "HG00099", [1, 1]),
            fst_line(1, 60000, "HG00096", [0, 0]),
            fst_line(1, 60000, "HG00097", [0, 0]),
            fst_line(1, 60000, "HG00098", [0, 0]),
            fst_line(1, 60000, "HG00099", [0, 1]),
        ]).unwrap();

        let frequency_set = Set::from_iter(vec![
            frq_line(1, 50000, "AFR", 0.0),
            frq_line(1, 50000, "EUR", 0.5),
            frq_line(1, 60000, "AFR", 1.0),
            frq_line(1, 60000, "EUR", 0.125),
        ]).unwrap();

        FSTReader::<Vec<u8>>{genotypes_set, frequency_set, genotypes: AHashMap::new(), frequencies: AHashMap::new()}
    }

    #[test]
    fn test_get_pop_allele_frequency() {
        let mut reader = fake_fst();
        
        reader.search_coordinate_frequencies(&Coordinate::new(1, 50_000));
        assert!(reader.get_pop_allele_frequency("EUR").is_ok_and(|freq| freq == 0.5));
        assert!(reader.get_pop_allele_frequency("AFR").is_ok_and(|freq| freq == 0.0));
        
        reader.search_coordinate_frequencies(&Coordinate::new(1, 60_000));
        assert!(reader.get_pop_allele_frequency("EUR").is_ok_and(|freq| freq == 0.125));
        assert!(reader.get_pop_allele_frequency("AFR").is_ok_and(|freq| freq == 1.0));
    }

    #[test]
    fn test_get_pop_allele_frequency_missing() {
        let mut reader = fake_fst();
        reader.search_coordinate_frequencies(&Coordinate::new(1, 50_000));
        let alleles = reader.get_pop_allele_frequency("SAS"); // Pop not found 
        assert!(alleles.is_err_and(|e| {
            matches!(e.downcast_ref::<GenotypeReaderError>(), Some(GenotypeReaderError::MissingFreq(_)))
        }));
    }
    #[test]
    fn test_get_allele() {
        let mut reader = fake_fst();
    
        reader.search_coordinate_genotypes(&Coordinate::new(1, 50_000));
        for (sample, want) in [("HG00096", [0, 0]), ("HG00097", [0, 1]), ("HG00098", [1, 0]), ("HG00099", [1, 1])] {
            let tag = SampleTag::new(sample, None, None);
            let alleles = reader.get_alleles(&tag);
            assert!(alleles.is_ok_and(|alleles| alleles == want));
        }
        
        reader.search_coordinate_genotypes(&Coordinate::new(1, 60_000));
        for (sample, want) in [("HG00096", [0, 0]), ("HG00097", [0, 0]), ("HG00098", [0, 0]), ("HG00099", [0, 1])] {
            let tag = SampleTag::new(sample, None, None);
            let alleles = reader.get_alleles(&tag);
            assert!(alleles.is_ok_and(|alleles| alleles == want));
        }
    }

    #[test]
    fn test_get_allele_missing() {
        let mut reader = fake_fst();
        reader.search_coordinate_genotypes(&Coordinate::new(1, 50_000));
        let alleles = reader.get_alleles(&SampleTag::new("NA0206", None, None)); // Sample not found 
        assert!(alleles.is_err_and(|e| {
            matches!(e.downcast_ref::<GenotypeReaderError>(), Some(GenotypeReaderError::MissingAlleles))
        }));
    }

    #[test]
    fn test_clear_buffers() {
        let mut reader = fake_fst();

        assert!(reader.frequencies.is_empty());
        assert!(reader.genotypes.is_empty());

        reader.search_coordinate_genotypes(&Coordinate::new(1, 50_000));
        reader.search_coordinate_frequencies(&Coordinate::new(1, 50_000));
        assert!(!reader.frequencies.is_empty());
        assert!(!reader.genotypes.is_empty());

        reader.clear_buffers();
        assert!(reader.frequencies.is_empty());
        assert!(reader.genotypes.is_empty());
    }

    #[test]
    fn test_has_genotypes() {
        let mut reader = fake_fst();
        assert!(!reader.has_genotypes());
        
        reader.search_coordinate_frequencies(&Coordinate::new(1, 50_000));
        assert!(!reader.has_genotypes());

        reader.search_coordinate_genotypes(&Coordinate::new(1, 50_000));
        assert!(reader.has_genotypes());
        
        reader.clear_buffers();
        assert!(!reader.has_genotypes());
    }

    #[test]
    fn test_find_chromosomes() {
        let genotypes_set = Set::from_iter(vec![
            fst_line(1, 50000, "HG00096", [0, 0]),
            fst_line(3, 50000, "HG00098", [1, 0]),
            fst_line(5, 60000, "HG00096", [0, 0]),
            fst_line(7, 60000, "HG00097", [0, 0]),
            fst_line(9, 60000, "HG00098", [0, 0]),
            fst_line(21, 60000, "HG00099", [0, 1]),
        ]).unwrap();

        let frequency_set = Set::from_iter(vec![frq_line(1, 50000, "AFR", 0.0)]).unwrap();

        let reader = FSTReader::<Vec<u8>>{genotypes_set, frequency_set, genotypes: AHashMap::new(), frequencies: AHashMap::new()};

        assert!(reader.find_chromosomes().is_ok_and(|v| v == vec![1,3,5,7,9,21]));
    }

    #[test]
    fn test_match_input_frq() -> anyhow::Result<()> {
        let tmpdir = tempfile::tempdir()?;
        let fst_path = tmpdir.path().join(format!("genome.{FST_EXT}"));
        let _ = File::create(&fst_path)?;
        assert!(FSTReader::<Vec<u8>>::match_input_frq(&fst_path).is_err_and(|e| {
            matches!(e.downcast_ref::<FSTReaderError>(), Some(FSTReaderError::MatchFrqFile { path: _}))
        }));

        let frq_path = tmpdir.path().join(format!("genome.{FRQ_EXT}"));
        let _ = File::create(frq_path)?;
        assert!(FSTReader::<Vec<u8>>::match_input_frq(&fst_path).is_ok());

        Ok(())
    }
}
