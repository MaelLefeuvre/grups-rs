use indicatif::MultiProgress;
use indicatif_log_bridge::LogWrapper;
use log::LevelFilter;
use log::Level;
use env_logger::{Builder, Env, fmt::Color};
use std::io::Write;
use std::sync::OnceLock;
use thiserror::Error;

static INSTANCE: OnceLock<Logger> = OnceLock::new();

#[derive(Debug)]
pub struct Logger {
    multi_pg: MultiProgress,
}


#[derive(Debug, Error)]
pub enum LoggerError {
    #[error("Logger is already initialized. ({0})")]
    DoubleLogInitialization(#[from] log::SetLoggerError),

    #[error("Static Instance of Logger is already set")]
    DoubleCellInitialization
}

impl Logger {
    pub fn init(verbosity: u8) -> Result<(), LoggerError>{
        let log_level = Self::u8_to_loglevel(verbosity);
        let env = Env::default()
            .filter("GRUPS_LOG");

        let logger = Builder::new().filter_level(log_level)
            .format(|buf, record| {
                
                let traceback: String;
                let set_intensity: bool;
                if record.level() == LevelFilter::Error {
                    traceback = format!("(@ {}:{}) ", record.file().unwrap_or("unknown"), record.line().unwrap_or(0));
                    set_intensity = true;
                } else {
                    traceback = String::from("");
                    set_intensity = false;
                };

                let mut arg_style = buf.style();
                arg_style.set_intense(set_intensity);


                let mut level_style = buf.style();
                let color = match record.level() {
                    Level::Error => Color::Red,
                    Level::Warn  => Color::Yellow, 
                    Level::Info  => Color::Green,
                    Level::Debug => Color::Blue,
                    Level::Trace => Color::Cyan
                };
                level_style.set_color(color).set_bold(true);

                writeln!(
                    buf,
                    "[{} {: <5} {}] {traceback}{}",
                    chrono::Local::now().format("%Y-%m-%dT%H:%M:%S"),
                    level_style.value(record.level()),
                    record.target(),
                    arg_style.value(record.args())
                )
            })
            .parse_env(env)
            .build();
        // Progress bar support.
        let multi_pg = MultiProgress::new();
        LogWrapper::new(multi_pg.clone(), logger)
            .try_init()?;

        INSTANCE.set(Self{multi_pg}).map_err(|_| LoggerError::DoubleCellInitialization)
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

    pub fn set_level(verbosity: u8) {
        log::set_max_level(Self::u8_to_loglevel(verbosity));
    }

    pub fn multi() -> &'static MultiProgress {
        &INSTANCE.get().expect("Unitialized").multi_pg
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_level() -> Result<(), LoggerError> {
        Logger::init(0)?;
        for level in 0..u8::MAX {
            Logger::set_level(level);

            let expected_level = match level {
                0           => LevelFilter::Error,
                1           => LevelFilter::Warn,
                2           => LevelFilter::Info,
                3           => LevelFilter::Debug,
                4..=u8::MAX => LevelFilter::Trace
            };

            assert_eq!(log::max_level(), expected_level);
        }
        Ok(())
    }
}