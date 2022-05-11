mod fstreader;
pub use fstreader::FSTReader;

use std::{
    path::{PathBuf, Path},
};

use log::{info, debug};

fn match_input_frq(fst_file: &Path) -> std::io::Result<()> {
    let frq = fst_file.with_extension("fst.frq");
    if ! frq.exists() {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("Could not find any valid paired '.fst.frq' file for {}\n
        Note that .fst files must have a matching .fst.frq file within the same location.", fst_file.to_str().unwrap_or("None"))))
    }
    Ok(())
}

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
    info!("Found input fst file candidates: {:#?}", fsts);

    // Ensure there's a matching 'fst.frq' file for each found '.fst' file
    for fst in fsts.iter() {
        match_input_frq(fst)?;
    }

    if fsts.is_empty() {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("Could not find any valid input fst file within {}\n
        Please specify an appropriate directory using `--data-dir`.
        Note that files are searched by matching files ending with '.fst', and '.fst.frq'.", input_dir.to_str().unwrap_or("None"))))
    }

    Ok(fsts)
}