use std::{fmt::Display, panic::Location};

use anyhow::{Context, Result};

pub mod prelude {
    extern crate anyhow;
    pub use anyhow::{anyhow, bail, Context, Result};
    
    extern crate thiserror;
    pub use thiserror::Error;

    pub use super::{LocatedError, LocatedOption};
}

macro_rules! loc_caller {
    ($caller:expr) => {
        format!("[{}:{}:{}]", $caller.file(), $caller.line(), $caller.column())
    }
}

pub trait LocatedError<T, E> {
    /// Wrap the error value with additional context + the location at which it was called.
    fn loc<C>(self, context: C) -> Result<T, anyhow::Error>
    where
        C: Display + Send + Sync + 'static;

    /// Wrap the error value with additional context that is evaluated lazily
    /// only once an error does occur + the location at which it was called.
    fn with_loc<C, F>(self, f: F) -> Result<T, anyhow::Error>
    where
        C: Display + Send + Sync + 'static,
        F: FnOnce() -> C;
}

impl<T, E> LocatedError<T, E> for Result<T, E>
where
    E: Display + Send + Sync + 'static,
    Result<T, E>: Context<T, E>,
{
    #[track_caller]
    fn loc<C>(self, context: C) -> Result<T, anyhow::Error>
    where
        C: Display + Send + Sync + 'static
    {
        match self {
            Ok(ok) => return Ok(ok),
            Err(_) => {
                let loc = loc_caller!(Location::caller());
                self.context(format!("{loc} {}", context))
            }
        }
    }

    #[track_caller]
    fn with_loc<C, F>(self, f: F) -> Result<T, anyhow::Error>
    where
        C: Display + Send + Sync + 'static,
        F: FnOnce() -> C
    {
        match self {
            Ok(ok) => return Ok(ok),
            Err(_) => {
                let caller = std::panic::Location::caller();
                let loc = format!("[{}:{}:{}]", caller.file(), caller.line(), caller.column());
                self.with_context( || format!("{loc} {}", f()))
            }
        }
        
    }
}


pub trait LocatedOption<T> {
    /// Wrap the error value with additional context + the location at which it was called.
    fn loc<C>(self, context: C) -> Result<T, anyhow::Error>
    where
        C: Display + Send + Sync + 'static;

    /// Wrap the error value with additional context that is evaluated lazily
    /// only once an error does occur + the location at which it was called.
    fn with_loc<C, F>(self, f: F) -> Result<T, anyhow::Error>
    where
        C: Display + Send + Sync + 'static,
        F: FnOnce() -> C;
}


impl<T> LocatedOption<T> for Option<T>

{
    #[track_caller]
    fn loc<C>(self, context: C) -> Result<T, anyhow::Error>
    where
        C: Display + Send + Sync + 'static
    {
        match self {
            Some(ok) => return Ok(ok),
            None => {
                let loc = loc_caller!(Location::caller());
                self.context(format!("{loc} {}", context))
            }
        }
    }

    #[track_caller]
    fn with_loc<C, F>(self, f: F) -> Result<T, anyhow::Error>
    where
        C: Display + Send + Sync + 'static,
        F: FnOnce() -> C
    {
        match self {
            Some(ok) => return Ok(ok),
            None => {
                let caller = std::panic::Location::caller();
                let loc = format!("[{}:{}:{}]", caller.file(), caller.line(), caller.column());
                self.with_context( || format!("{loc} {}", f()))
            }
        }
        
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use thiserror::Error;
    use std::{fs::File};

    #[derive(Error, Debug)]
    pub enum BubbleError {
        #[error(transparent)]
        Bloup(#[from] anyhow::Error),

        #[error("Blop")]
        Blop,
    }

    fn error_source(file: &str) -> Result<(), BubbleError> {
        let _ = File::open(file).loc(format!("failed to open file: '{file}'"))?;
        Ok(())
    }

    fn error_bubble_0() -> Result<()> { 
        error_source("bloup/").with_loc(|| "Whoops!")
    }

    fn error_bubble_1() -> Result<()> {
        error_bubble_0().with_loc(|| BubbleError::Blop)
    }

    fn error_bubble_2() -> Result<()> {
        error_bubble_1().with_loc(|| "WTF?!")
    }

    #[test]
    fn print_with_loc_error() -> Result<()> {
        if let Err(err) = error_bubble_2() {
            // ---- Ensure file, line, and col matches.
            let mut chain = err.chain().into_iter();
            let results = [error_bubble_2(), error_bubble_1(), error_bubble_0()];
            for result in results { 
                assert_eq!(
                    format!("{}", chain.next().unwrap()),
                    format!("{}", result.err().unwrap())
                );
            }

            // ---- Display
            eprintln!("ERROR: {err:?}");
        }
        Ok(())
    }

    fn none_bubble() -> Option<()> {
        None
    }

    #[test]
    fn test_missing() -> Result<()> {
        let x = none_bubble().loc(BubbleError::Blop);
        eprintln!("{x:?}");
        Ok(())
    }
}