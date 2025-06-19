use std::{io, env, fs, path::{Path, PathBuf}, ops::Deref, fmt::{self, Formatter, Display}};
use tempfile::{self, TempDir};

pub const TEST_DATA_DIR: &str = "./tests/test-data";

pub struct Fixture {
    path: std::path::PathBuf,
    source: std::path::PathBuf,
    _tempdir: TempDir,
}

impl Fixture {
    pub fn blank(fixture_filename: &str) -> Self {
        // First, figure out the right file in `tests/fixtures/`:
        let root_dir = &env::var("CARGO_MANIFEST_DIR").expect("$CARGO_MANIFEST_DIR");
        let mut source = PathBuf::from(root_dir);
        source.push(TEST_DATA_DIR);
        source.push(fixture_filename);

        // The "real" path of the file is going to be under a temporary directory:
        let tempdir = tempfile::tempdir().expect("Failed to generate temp directory");

        // Check if the given path is absolute... match only the filename if so. (allows recursion.)
        let mut path = PathBuf::from(&tempdir.path());
        let fixture_path = Path::new(&fixture_filename);
        match fixture_path.is_absolute() {
            true  => path.push(fixture_path.file_name().expect("Invalid filename")), // PathBuf overwrites buffer when path is absolute.
            false => path.push(fixture_filename), 
        };

        Fixture { _tempdir: tempdir, source, path }
    }
}

impl Fixture {
    pub fn copy(fixture_filename: &str) -> Self {
        let fixture = Fixture::blank(fixture_filename);
        if fixture.source.is_dir() {
            copy_dir_all(&fixture.source, &fixture.path).expect("Failed to copy directory");
        } else {
            fs::create_dir_all(fixture.path.parent().expect("No parent directory")).expect("Failed to create directory");
            fs::copy(&fixture.source, &fixture.path).expect("Failed to copy Fixture files.");
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

impl Deref for Fixture {
    type Target = Path;

    fn deref(&self) -> &Self::Target {
        self.path.deref()
    }
}

impl Display for Fixture {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.path.to_str().expect("Invalid path (non UTF8 characters ?)"))
    }
}