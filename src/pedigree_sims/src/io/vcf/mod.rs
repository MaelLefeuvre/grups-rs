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

    if vcfs.is_empty() {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("Could not find any valid input vcf/vcf.gz files within {}\n
        Please specify an appropriate directory using `--data-dir`.
        Note that files are searched by matching files ending with '.vcf', or '.vcf.gz'.", input_dir.to_str().unwrap_or("None"))))
    }


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

    match panel.len() {
        1 => (),
        0 => { return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("Could not find a candidate panel definition file within {}\n
                Please specify the relevant file using '--panel.
                Exiting.", input_dir.to_str().unwrap_or("None"))))
        },
        _ =>  { return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("Found multiple candidate Panel definition files: {:#?}\n
            Please specify the relevant file using '--panel'.
            Exiting.", panel)
        ))},
    }

    info!("Found: {}", panel[0].to_str().unwrap_or("None"));
    Ok(panel[0].clone())
}