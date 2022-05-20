use std::path::Path;
use std::fs;
use std::io;

use tempfile;
use tempfile::TempDir;

pub const TEST_DATA_DIR: &str = "./tests/test-data";

pub struct Fixture {
    path: std::path::PathBuf,
    source: std::path::PathBuf,
    _tempdir: TempDir,
}

impl Fixture {
    pub fn blank(fixture_filename: &str) -> Self {
        // First, figure out the right file in `tests/fixtures/`:
        let root_dir = &std::env::var("CARGO_MANIFEST_DIR").expect("$CARGO_MANIFEST_DIR");
        let mut source = std::path::PathBuf::from(root_dir);
        source.push(TEST_DATA_DIR);
        source.push(&fixture_filename);

        // The "real" path of the file is going to be under a temporary directory:
        let tempdir = tempfile::tempdir().unwrap();

        // Check if the given path is absolute... match only the filename if so. (allows recursion.)
        let mut path = std::path::PathBuf::from(&tempdir.path());
        let fixture_path = Path::new(&fixture_filename);
        match fixture_path.is_absolute() {
            true => path.push(&fixture_path.file_name().unwrap()), // PathBuf overwrites the buffer if the given str is absolute.
            false => path.push(&fixture_filename), 
        };

        Fixture { _tempdir: tempdir, source, path }
    }
}

impl Fixture {
    pub fn copy(fixture_filename: &str) -> Self {
        let fixture = Fixture::blank(fixture_filename);
        if fixture.source.is_dir() {
            copy_dir_all(&fixture.source, &fixture.path).unwrap();
        } else {
            std::fs::create_dir_all(&fixture.path.parent().unwrap()).unwrap();
            std::fs::copy(&fixture.source, &fixture.path).unwrap();
        }
        fixture
    }
}

fn copy_dir_all(src: impl AsRef<Path>, dst: impl AsRef<Path>) -> io::Result<()> {
    fs::create_dir_all(&dst)?;
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let ty = entry.file_type()?;
        if ty.is_dir() {
            copy_dir_all(entry.path(), dst.as_ref().join(entry.file_name()))?;
        } else {
            fs::copy(entry.path(), dst.as_ref().join(entry.file_name()))?;
        }
    }
    Ok(())
}

impl std::ops::Deref for Fixture {
    type Target = std::path::Path;

    fn deref(&self) -> &Self::Target {
        self.path.deref()
    }
}

impl std::fmt::Display for Fixture {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.path.to_str().unwrap())
    }
}