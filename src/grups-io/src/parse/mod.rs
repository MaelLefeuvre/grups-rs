use std::{collections::HashMap, fs, path::{Path, PathBuf}};

use located_error::LocatedError;

use log::trace;
use anyhow::Result;

mod error;
use error::ParseError;

/// Attempt to create the parent directories of a path (if needed) and return an error if it failed.
pub fn create_parent_directory(path: &Path) -> Result<()> {
    use ParseError::CreateParentDirectory;
    let parent_dir = path.parent().unwrap_or(path);
    let loc_msg = || format!("While attempting to create output directory '{}'", path.display());
    fs::create_dir_all(parent_dir).map_err(CreateParentDirectory).with_loc(loc_msg)?;
    Ok(())
}

/// Attempt to convert a path to string, and return an error if it failed.
fn maybe_to_str(path: &Path) -> Result<&str> {
    use ParseError::InvalidFilename;
    path.to_str().ok_or(InvalidFilename).loc("While converting path to string")
}

/// Simple enum for `get_output_files()`
///  Suffix => `HashMap` will use provided file suffixes as keys
///  Key    => `HashMap` will use provided file extensions as keys 
#[derive(Clone, Copy)]
pub enum FileKey{Suffix, Ext}

/// Obtain predefined filenames for jackknife blocks from the given output directory and pileup file. 
/// Return a `HashMap` with (K (pair): "{ind1}-{ind2}", V (File): "{out dir}/blocks/{file prefix}.block"
/// ==> A new file is generated for each pair of individuals. 
/// 
/// # Errors
/// 
/// If creating the parent directory of `file_prefix` is required, will throw a `PermissionDenied` if
/// the user does not have the proper UNIX permissions 
/// 
/// # Panics
/// 
/// when failing to convert `file_prefix` from `PathBuf` to `&str`
/// 
pub fn get_output_files<'a>(
    file_prefix    : &'a mut Path,
    allow_overwrite: bool,
    sort           : FileKey,
    suffixes       : &[String],
    file_ext       : &[&'a str]
) -> Result<HashMap<String, String>> {
    let err_context = "While attempting to format the name of the output_files";
    // Create output directory. Early return if we can't create it
    create_parent_directory(file_prefix)?;
    // Generate a HashMap of filepaths from the file_prefix
    let mut outfiles_hash = HashMap::with_capacity(suffixes.len() * file_ext.len());
    for suffix in suffixes {
        let mut file = PathBuf::from({
            /*final dot to fake an extension*/
            let prefix = maybe_to_str(file_prefix).loc(err_context)?;
            let sep = if suffix.is_empty() {""} else {"-"};
            format!("{prefix}{sep}{suffix}.") // Final dot to fake an extension.
        });
        for ext in file_ext {
            file.set_extension(ext);
            can_write_file(allow_overwrite, &file)?;
            let key = match &sort {
                FileKey::Suffix => suffix.clone(),
                FileKey::Ext    => (*ext).to_string()
            };
            let outfile = maybe_to_str(&file).loc(err_context)?;
            outfiles_hash.insert(key, outfile.to_string());
        }
    }

    trace!("Output File(s): {:#?}", outfiles_hash.values());
    Ok(outfiles_hash)
}

/// Check if a given file already exists ; raise an error if such is the case, and the user did not explicitly 
/// allow file overwriting.
/// # Errors
/// - If the provided `pathbuf` already exists and the user did not specifically allow for file
///   overwrite using the `--overwrite` argument
/// 
/// # Panics
/// - if the provided `pathbuf` fails to get parsed as a string.
/// 
/// # @TODO:
/// - This method is a duplicate of `parser::can_write_file`. I/O functions, methods and structs should be contained
///   is their own package.
pub fn can_write_file(overwrite: bool, path: &Path) -> Result<bool> {
    let loc_msg = "While ensuring that file permissions were appropriate";
    if !overwrite && path.exists() {   // Check if this file already exists and/or if overwrite is allowed.
        return Err(ParseError::OverwriteDisallowed{path: path.to_path_buf()}).loc(loc_msg)
    }
    Ok(true)
}


/// Iterate over the contents of a OS-directory and search for all files matching a given list of extensions.
/// Return these matches as a vector of PathBuf.
/// 
/// # Arguments:
/// - `input_dir`: path of the targeted directory, were search should be performed.
/// 
/// # Errors:
/// - Returns an error if the output Vec<PathBuf> is empty after iterating over all the contents of `input_dir`
pub fn fetch_input_files(input_dir: &Path, extensions: &[&str] ) -> Result<Vec<PathBuf>>{
    use ParseError::MissingInput;
    // Fetch any file matching the provided extensions. Note that we don't use the '.extension()'
    // method, since we want to target chained-ext file, as in '.vcf.gz'.
    let mut files: Vec<PathBuf> = fs::read_dir(input_dir)?.filter_map(|file| {
        let Ok(file) = file else {
            return None
        };
        let file     = file.path();
        let str_file = maybe_to_str(&file).unwrap_or("");
        match extensions.iter().any(|ext| str_file.ends_with(ext)) {
            true  => Some(file),
            false => None,
        }
    }).collect::<Vec<PathBuf>>();

    files.sort();

    // There should be at least one file within that directory. If not, we must notify the user.
    let loc_msg = "Please ensure the requested directory contains the appropriate files, and that they all end with a valid file extension";
    match files.is_empty() { 
        true  => Err(MissingInput{dir: input_dir.to_path_buf()}).loc(loc_msg),
        false => Ok(files)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    
    #[test]
    fn test_fetch_input_files_ok() -> anyhow::Result<()> {
        // Create a directory inside of std::env::temp_dir()
        let tmpdir = tempfile::tempdir()?;

        let want = "targets.vcf.gz";
        for file in ["README.txt", want, "targets.vcf.gz.tbi"] {
            let path = tmpdir.path().join(file);
            let _ = File::create(path)?;
        }
        let provided_path = tmpdir.path();
        let files = fetch_input_files(provided_path, &[".vcf.gz"]).expect("Failed to fetch files");
        assert_eq!(files.len(), 1);
        assert_eq!(files[0].file_name().map(Path::new), Some(Path::new(want)));

        Ok(())
    }

    #[test]
    fn test_fetch_input_files_missing_input() -> anyhow::Result<()> {
        // Create a directory inside of `std::env::temp_dir()`
        let tmpdir = tempfile::tempdir()?;

        for file in ["README.txt", "targets.vcf.gz", "targets.vcf.gz.tbi"] {
            let path = tmpdir.path().join(file);
            let _ = File::create(path)?;
        }
        let provided_path = tmpdir.path();
        let files = fetch_input_files(provided_path, &[".vcf"]);
        assert!(files.is_err_and(|e| matches!(e.downcast_ref::<ParseError>(), Some(ParseError::MissingInput{dir: _}))));

        Ok(())
    }

    #[test]
    fn test_can_write_file() -> anyhow::Result<()> {
        let tmpdir = tempfile::tempdir()?;

        let path   = tmpdir.path().join("README.md");
        assert!(can_write_file(false, &path).is_ok_and(|x| x)); // No overwrite, no file => should return true
        assert!(can_write_file(true, &path).is_ok_and(|x| x));  // Overwrite, no file    => should return true

        let _   = File::create(&path)?;
        assert!(can_write_file(true, &path).is_ok_and(|x| x));  // Overwrite, file       => should return true
        assert!(can_write_file(false, &path).is_err_and(|e| {   // No overwrite, file       => should error
            matches!(e.downcast_ref::<ParseError>(), Some(ParseError::OverwriteDisallowed{path: _}))
        })); 


        Ok(())
    }
}