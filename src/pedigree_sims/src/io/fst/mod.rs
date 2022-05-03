mod fstreader;
pub use fstreader::FSTReader;

use std::{
    path::PathBuf,
};

use log::{info, debug};

pub fn get_input_fst(input_dir: &PathBuf) -> std::io::Result<Vec<PathBuf>>{
    let paths = std::fs::read_dir(input_dir)?;

    let mut fsts = paths.filter_map(Result::ok)
        .filter_map(|d| d.path()
            .to_str()
            .and_then(|f|
                if f.ends_with(".fst") { 
                    debug!("Found: {}", f);
                    Some(d) 
                } 
                else { 
                    debug!("Skipping: {}", f);
                    None
                }
            )
            .map(|f| f.path())
        )
    .collect::<Vec<PathBuf>>();
    fsts.sort();
    fsts.dedup();
    info!("Found input vcf file candidates: {:#?}", fsts);
    Ok(fsts)
}