pub mod reader;

mod sampletag;
pub use sampletag::SampleTag;

use log::{info, debug};
use std::{
    path::{PathBuf}
};

pub fn get_input_vcfs(input_dir: &PathBuf) -> std::io::Result<Vec<PathBuf>>{
    let paths = std::fs::read_dir(input_dir)?;

    let mut vcfs = paths.filter_map(Result::ok)
        .filter_map(|d| d.path()
            .to_str()
            .and_then(|f|
                if f.ends_with(".vcf") || f.ends_with(".vcf.gz") { 
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
    vcfs.sort();
    info!("Found input vcf file candidates: {:#?}", vcfs);
    Ok(vcfs)
}

pub fn fetch_input_panel(input_dir: &PathBuf) -> std::io::Result<PathBuf>{

    let paths = std::fs::read_dir(input_dir)?;
    let panel = paths.filter_map(Result::ok)
        .filter_map(|d| d.path()
            .to_str()
            .and_then(|f|
                if f.ends_with(".panel") { 
                    Some(d) 
                } 
                else { 
                    None
                }
            )
            .map(|f| f.path())
        )
    .collect::<Vec<PathBuf>>();

    if panel.len() != 1 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, 
            format!("Found multiple candidate Panel definition files: {:#?}\n
            Please specify the relevant file using '--panel'.
            Exiting.", panel)
        ))
    }

    info!("Found: {}", panel[0].to_str().unwrap_or("None"));
    Ok(panel[0].clone())
}