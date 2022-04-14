use log::LevelFilter;
use env_logger::{Builder, Env};

pub fn init_logger(verbosity: &u8){
    let log_level = match verbosity {
        0           => LevelFilter::Error,
        1           => LevelFilter::Warn,
        2           => LevelFilter::Info,
        3           => LevelFilter::Debug,
        4..=u8::MAX => LevelFilter::Trace

    };
    let env = Env::default()
        .filter("GRUPS_LOG");

    Builder::new().filter_level(log_level).parse_env(env).init()
}