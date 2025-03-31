use log::LevelFilter;
use log::Level;
use env_logger::{Builder, Env, fmt::Color};
use std::io::Write;
pub struct Logger;

impl Logger {

    pub fn init(verbosity: u8) {
        let log_level = Self::u8_to_loglevel(verbosity);
        let env = Env::default()
            .filter("GRUPS_LOG");

        Builder::new().filter_level(log_level)
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
            .init();
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_level(){
        Logger::init(0);
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
    }
}