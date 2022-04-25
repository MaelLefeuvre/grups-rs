use pedigree_sims::io;
use pwd_from_stdin::genome::SNPCoord;
use std::path::{PathBuf};
use std::error::Error;

use std::fs::File;
use fst::{IntoStreamer, Map, MapBuilder};
//use fst::Streamer;
use fst::automaton::{Automaton, Str};
//use fst::automaton::Subsequence;
//use memmap::Mmap;

use std::rc::Rc;
use std::collections::BTreeSet;
// [POP][SAMPLES][Nucleotide]

#[derive(Clone, Debug)]
pub struct VCFCoord {
    coordinate: SNPCoord,
    _af: f32,
}

impl VCFCoord {
    pub fn new(chromosome: u8, position: u32, reference: Option<char>, alternate: Option<char>, af: f32) -> VCFCoord {
        let coordinate = SNPCoord{chromosome, position, reference, alternate };
        VCFCoord{coordinate, _af: af}
    }
}

impl std::borrow::Borrow<u32> for VCFCoord {
    fn borrow(&self) -> &u32 {
        self.coordinate.position.borrow()
    }
}


impl PartialEq for VCFCoord {
    fn eq(&self, other: &VCFCoord) -> bool {
        self.coordinate == other.coordinate
    }
}

impl Eq for VCFCoord {}

impl std::cmp::Ord for VCFCoord {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.coordinate).cmp(&(other.coordinate))
    }
}

impl PartialOrd for VCFCoord {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug)]
pub struct Allele {
    coordinate: Rc<VCFCoord>,
    _haplos: (u8,u8)
}

impl Allele {
    pub fn new(coordinate: &Rc<VCFCoord>, haplos: (u8,u8)) -> Allele {
        Allele {coordinate: Rc::clone(coordinate), _haplos: haplos}
    }
}

impl PartialEq for Allele {
    fn eq(&self, other: &Allele) -> bool {
        self.coordinate == other.coordinate
    }
}

impl Eq for Allele {}

impl Ord for Allele {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.coordinate).cmp(&(other.coordinate))
    }
}

impl PartialOrd for Allele {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Debug)]
pub struct VCFSample {
    id: String,
    idx: usize,
    pub genome: BTreeSet<Allele>
}

impl VCFSample {
    pub fn new(id: String, idx: usize) -> VCFSample {
        VCFSample{ id, idx, genome: BTreeSet::new()}
    }

    pub fn insert_allele(&mut self, coordinate:  &Rc<VCFCoord>, haplos: (u8,u8)) {
        self.genome.insert(Allele::new(coordinate, haplos));
    }

    pub fn id(&self) -> &String {
        &self.id
    }
    pub fn idx(&self) -> &usize {
        &self.idx
    }
}

impl<'a> std::borrow::Borrow<String> for VCFSample {
    fn borrow(&self) -> &String {
        self.id.borrow()
    }
}

impl<'a> std::borrow::Borrow<usize> for VCFSample {
    fn borrow(&self) -> &usize {
        self.idx.borrow()
    }
}

#[derive(Debug, Default)]
pub struct VCF {
    pub coordinates: BTreeSet<Rc<VCFCoord>>,
    pub samples : Vec<VCFSample>
}

impl VCF {

    pub fn insert_coordinate(&mut self, chromosome: u8, position: u32, reference: char, alternate: char, af: f32) -> Rc<VCFCoord> {
        let coordinate = Rc::new(VCFCoord::new(chromosome, position, Some(reference), Some(alternate), af));
        let rcx = Rc::clone(&coordinate);

        self.coordinates.insert(coordinate);
        rcx

    }

    pub fn insert_sample(&mut self, id: String, idx: usize) {
        let sample = VCFSample::new(id, idx);
        self.samples.push(sample);
    }
}

fn _serialize(vcf_paths: &[PathBuf]) -> Result<(), Box<dyn Error>> {

    for vcf in vcf_paths {
        println!("{:?}", vcf);


        let mut vcf_hash = VCF::default();

        //Start a reader.
        let mut reader = io::VCFReader::new(vcf)?;

        //Extract samples
        let samples = reader.samples();
        let samples = &samples[9..];

        
        //Convert to struct to keep id.
        let mut struct_samples = Vec::new();
        for (i, sample) in samples.iter().enumerate(){
            let sample = io::SampleTag::new(sample, i);
            struct_samples.push(sample);
        }

        //Sort samples
        struct_samples.sort();

        for sample in struct_samples.iter() {
            vcf_hash.insert_sample(sample.id().clone(), *sample.idx())
        }

        let mut i = 0;
        while reader.has_data_left()? {
            let chr: u8  = reader.next_field()?.parse()?;
            reader.clear_buffer();
            let pos: u32 = reader.next_field()?.parse()?;
            reader.clear_buffer();
            reader.skip(1).unwrap();
            let reference: char = reader.next_field()?.parse()?;
            reader.clear_buffer();
            let alternate: char = reader.next_field()?.parse()?;
            reader.clear_buffer();

            let coord = vcf_hash.insert_coordinate(chr, pos, reference, alternate, 0.0);

            //let coord: &Rc<VCFCoord> = vcf_hash.coordinates.get().unwrap();

            reader.fill_genotypes()?;
            for sample in vcf_hash.samples.iter_mut() {
                let haplos = reader.get_alleles(sample.idx())?;
                //sample.insert_allele(&Rc::clone(coord), haplos);
                sample.genome.insert(Allele::new(&Rc::clone(&coord), haplos));
                //if i % 500000 == 0 {
                //    println!("{:?}", sample);
                    //println!("{} {: >2} {: >9} {} {}",
                    //    sample.idx(), chr, pos, haplos.0, haplos.1
                    //);
                //}
            }
            

            reader.clear_buffer();
            if i % 1000 == 0 {
                println!("{}", i);
            }
            i+=1;
        }


    }
    Ok(())
}

fn _sort_vcf_map(vcf_paths: &[PathBuf]) -> Result<(), Box<dyn Error>> {

    let wtr = std::io::BufWriter::new(File::create("tests/test-data/fst/g1k-map.fst")?);
    let mut build= MapBuilder::new(wtr)?;

    for vcf in vcf_paths {
        println!("{:?}", vcf);

        //Start a reader.
        let mut reader = io::VCFReader::new(vcf)?;

        //Extract samples
        let samples = reader.samples();
        let samples = &samples[9..];

        //Convert to struct to keep id.
        let mut struct_samples = Vec::new();
        for (i, sample) in samples.iter().enumerate(){
            let sample = io::SampleTag::new(sample, i);
            struct_samples.push(sample);
        }

        //Sort
        struct_samples.sort();

        //Read vcf
        let (mut buf, mut genotypes) = (Vec::new(), Vec::new());
        //for (i, line) in reader.lines().enumerate() {
            let mut i = 0;
        loop {
            let chr_bytes = reader.source.read_until(b'\t', &mut buf)?; // 1
            if chr_bytes == 0 { break };
            buf.pop(); buf.push(b' ');
            let pos_bytes = reader.source.read_until(b'\t', &mut buf)?; // 2
            buf.pop(); buf.push(b' ');
            for _ in 0..(10-pos_bytes) {                 // Add 9 leading zeroes.
                buf.insert(chr_bytes, b'0')   
            }

            for _ in 3..10 {
                reader.source.read_until(b'\t', &mut Vec::new())?;
            }
            
            reader.source.read_until(b'\n', &mut genotypes)?;

            if i % 50000 == 0 {
                println!("{i: >9} {:?}",std::str::from_utf8(&buf)?);
            }
            unsafe{
            for sample in struct_samples.iter() {
                let geno_idx=sample.idx()*4;
                let genotype = genotypes[geno_idx..geno_idx+3].to_owned();
                let mut key = buf.clone();
                key.append(sample.id().clone().as_mut_vec()); //Unsafe! 

                let _haplo1 = genotypes[geno_idx]   - 64;
                let _haplo2 = genotypes[geno_idx+2] - 64;

                let val = std::str::from_utf8(&[genotype[0]+1, genotype[2]+1])?.parse::<u64>()?;
                //println!("{:?} {:?}", std::str::from_utf8(&key), val);
                build.insert(key, val).unwrap();

            }}
            genotypes.clear();
            buf.clear();
            i+=1;

        }
    }
    build.finish()?;
    Ok(())
}

fn main() {

    let data_dir = PathBuf::from("tests/test-data/vcf/g1k-phase3-v5b-first-50k-filtered/");
    println!("Fetching input VCF files in {}", &data_dir.to_str().unwrap());
    let mut input_vcf_paths = io::get_input_vcfs(&data_dir).unwrap();
    input_vcf_paths.sort();


    let panel = PathBuf::from("tests/test-data/vcf/g1k-phase3-v5b/integrated_call_samples_v3.20130502.ALL.panel");
    let _panel = io::VCFPanelReader::new(panel.as_path(), input_vcf_paths[0].as_path()).unwrap();

    //for pop in panel.samples.into_keys(){
    //    println!("{}", pop); 
    //}




    // Bincode strategy

    //serialize(&input_vcf_paths, &tpool).unwrap();

    // FST Strategy
    println!("Generating FST index...");
    //sort_vcf_map(&input_vcf_paths, &tpool).unwrap();
    println!("Done!");

    //// Memory Map strategy.
    //println!("Opening memory map");
    //let mmap = unsafe { Mmap::map(&File::open("tests/fst/g1k-map.fst").unwrap()).unwrap() };
    //let map = Map::new(mmap).unwrap();

    //// In RAM Strategy
    println!("Reading in memory.");
    let mut file_handle = File::open("tests/test-data/fst/g1k-map.fst").unwrap();
    let mut bytes = vec![];
    std::io::Read::read_to_end(&mut file_handle, &mut bytes).unwrap();
    let map = Map::new(bytes).unwrap();


    let samples = ["HG01067", "NA20885", "HG00267", "HG00236", "NA20356", "NA19346"];

    for sample in samples.iter() {
        println!("Searching {sample}");
        let regex = format!("1 000055164 {sample}");
        //let matcher = Subsequence::new(&regex);
        let matcher = Str::new(&regex)
            .starts_with();
            //.intersection(Subsequence::new(&sample));
        let stream = map.search(&matcher).into_stream();

        println!("  - Populating vector");
        //let mut kvs = vec![];
        
        let kvs = stream.into_str_vec().unwrap();
        println!("{kvs:?}");
        //while let Some((k, v)) = stream.next() {
        //    println!("{:?} {}", k, v);
        //    let record = (String::from_utf8(k.to_vec()).unwrap(), v);
        //    println!("{:?}", &record);
        //    kvs.push(record);
        //    break
        //}
        //println!("{kvs:?}\n");
    }
}
