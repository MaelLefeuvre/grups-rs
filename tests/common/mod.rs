#[cfg(test)]
mod fixture;
pub use fixture::Fixture;

#[cfg(test)]
mod grups_runner;
pub use grups_runner::GrupsRunnerBuilder;

#[macro_export]
macro_rules! validate_file {
    ($ref_file:expr, $obtained_file:expr) => {
        let want = include_bytes!($ref_file);
        let got  = std::fs::read($obtained_file)
            .unwrap_or_else(|_| panic!("Failed to open {:?}", $obtained_file));
        assert_eq!(want.to_vec(), got)
    };
}

