use pedigree_sims::io;
use std::path::{PathBuf};
use std::error::Error;

use std::fs::File;
use fst::{IntoStreamer, Streamer, Map, MapBuilder};
use fst::automaton::Subsequence;
use memmap::Mmap;

use rust_htslib::tpool::ThreadPool;



fn sort_vcf_map(vcf_paths: &[PathBuf], tpool: &ThreadPool) -> Result<(), Box<dyn Error>> {

    let wtr = std::io::BufWriter::new(File::create("tests/fst/g1k-map.fst")?);
    let mut build= MapBuilder::new(wtr)?;


    for vcf in vcf_paths {
        println!("{:?}", vcf);

        //Start a reader.
        let mut reader = io::VCFReader::new(vcf, tpool)?;

        //Extract samples
        let samples = reader.samples();
        let samples = &samples[9..];

        //Convert to struct to keep id.
        let mut struct_samples = Vec::new();
        for (i, sample) in samples.iter().enumerate(){
            let sample = io::Sample::new(sample, i);
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
            //let line = line?;
            //let line = line.split('\t');

            //match i % 50000 {
            //    0 => println!("{}", line[1]),
            //    _ => ()
            //}
            //let chr = line.next(); // 1
            //let pos = line.next(); // 2

            //line.skip(7);                       // 9
            //let line = line.to_str();

            if i % 50000 == 0 {
                println!("{i: >9} {:?}",std::str::from_utf8(&buf)?);
            }
            unsafe{
            for sample in struct_samples.iter() {
                let geno_idx=sample.idx()*4;
                let genotype = genotypes[geno_idx..geno_idx+3].to_owned();
                let mut key = buf.clone();
                key.append(sample.id().clone().as_mut_vec()); //Unsafe! 

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

    let tpool = rust_htslib::tpool::ThreadPool::new(16).unwrap();



    let data_dir = PathBuf::from("tests/test-data/vcf/g1k-phase3-v5b-first-5000-filtered/");
    println!("Fetching input VCF files in {}", &data_dir.to_str().unwrap());
    let mut input_vcf_paths = io::get_input_vcfs(&data_dir).unwrap();
    input_vcf_paths.sort();


    let panel = PathBuf::from("tests/test-data/vcf/g1k-phase3-v5b/integrated_call_samples_v3.20130502.ALL.panel");
    let _panel = io::VCFPanelReader::new(panel.as_path(), input_vcf_paths[0].as_path(), &tpool).unwrap();

    //for pop in panel.samples.into_keys(){
    //    println!("{}", pop); 
    //}

    println!("Generating FST index...");
    sort_vcf_map(&input_vcf_paths, &tpool).unwrap();
    println!("Done!");


    // Memory Map strategy.
    println!("Opening memory map");
    let mmap = unsafe { Mmap::map(&File::open("tests/fst/g1k-map.fst").unwrap()).unwrap() };
    let map = Map::new(mmap).unwrap();

    //// In RAM Strategy
    //println!("Reading in memory.");
    //let mut file_handle = File::open("tests/fst/g1k-map.fst").unwrap();
    //let mut bytes = vec![];
    //std::io::Read::read_to_end(&mut file_handle, &mut bytes).unwrap();
    //let map = Map::new(bytes).unwrap();


    println!("Searching HG02067");
    let matcher = Subsequence::new("HG01067");
    let mut stream = map.search(&matcher).into_stream();

    println!("Populating vector");
    let mut kvs = vec![];
    while let Some((k, v)) = stream.next() {
        let record = (String::from_utf8(k.to_vec()).unwrap(), v);
        //println!("{:?}", &record);
        kvs.push(record);
    }
    //println!("{kvs:#?}");
}
