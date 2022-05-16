use log::LevelFilter;
use env_logger::{Builder, Env};

pub struct Logger;

impl Logger {

    pub fn init(verbosity: u8) -> Self {
        let log_level = Self::u8_to_loglevel(verbosity);
        let env = Env::default()
            .filter("GRUPS_LOG");

        Builder::new().filter_level(log_level).parse_env(env).init();
        Self
    }

    fn u8_to_loglevel(verbosity: u8) -> LevelFilter {
        match verbosity {
            0            => LevelFilter::Error,
            1            => LevelFilter::Warn,
            2            => LevelFilter::Info,
            3            => LevelFilter::Debug,
            4..= u8::MAX => LevelFilter::Trace
        }
    }

    pub fn set_level(&self, verbosity: u8) {
        log::set_max_level(Self::u8_to_loglevel(verbosity))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_level(){
        let logger = Logger::init(0);
        for level in 0..u8::MAX {
            logger.set_level(level);

            let expected_level = match level {
                0           => LevelFilter::Error,
                1           => LevelFilter::Warn,
                2           => LevelFilter::Info,
                3           => LevelFilter::Debug,
                4..=u8::MAX => LevelFilter::Trace
            };

            assert_eq!(log::max_level(), expected_level);
        }
    }
}