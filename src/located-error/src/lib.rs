use std::{fmt::Display, panic::Location};

use anyhow::{Context, Result};

/// Public prelude for LocatedError.
/// 
/// Note that this re-exports anyhow and thiserror
pub mod prelude {
    extern crate anyhow;
    pub use anyhow::{anyhow, bail, Context, Result};
    
    extern crate thiserror;
    pub use thiserror::Error;

    pub use super::{loc, LocatedError, LocatedOption};
}

macro_rules! loc_caller {
    ($caller:expr) => {
        format!("[{}:{}:{}]", $caller.file(), $caller.line(), $caller.column())
    }
}

#[macro_export]
macro_rules! loc {
    ($e: expr) => {
        Err(anyhow::anyhow!(format!("[{}:{}:{}] {}", file!(), line!(), column!(), $e)))
    }
}

/// Trait extending [`anyhow::Result<T>`] with additional information regarding the location of the error (e.g. file, line and column)
/// 
/// # Example 
/// ```should_panic
/// use anyhow::{anyhow, Result};
/// use crate::located_error::LocatedError;
/// 
/// // ---- Main runner
/// fn compute(path: &str) -> Result<()> {
///     let path: &str = "/invalid-file/";
///     let file = std::fs::File::open(path)
///         .with_loc(|| format!("Failed to open file {path}") )?;
///     /* ---- expensive computation ensues ---- */
///     Ok(())
/// }
/// 
/// // ---- Main
/// fn main() -> Result<()> {
///     let path: &str = "/invalid-file/";
///     let file = compute(path).loc("While running main function.")?;
///     Ok(())
/// }
/// ```
/// ## This should output the following lines 
/// ```Text
/// > Error: [src/lib.rs:14:26] While runnning main function.
/// > 
/// > Caused by:
/// >     0: [src/lib.rs:8:10] Failed to open file /invalid-file/
/// >     1: No such file or directory (os error 2)
/// ```
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
    /// Implement [`LocatedError`] for any [`Result<T, std::error::Error>`]
    /// 
    /// Note that this will will inevitably convert your error into  an `anyhow::Result<T>`  
    /// 
    /// Furthermore, note that [`LocatedError::loc()`] is eagerly evaluated. 
    /// For a lazy implementation, see [`LocatedError::with_loc()`]
    #[track_caller]
    fn loc<C>(self, context: C) -> Result<T, anyhow::Error>
    where
        C: Display + Send + Sync + 'static
    {
        match self {
            Ok(ok) => Ok(ok),
            Err(_) => {
                let loc = loc_caller!(Location::caller());
                self.context(format!("{loc} {context}"))
            }
        }
    }

    /// Implement [`LocatedError`] for any [`Result<T, std::error::Error>`]
    /// 
    /// Note that this will will inevitably convert your error into  an `anyhow::Result<T>`  
    /// 
    /// Furthermore, note that [`LocatedError::with_loc()`] is lazily evaluated. 
    /// For an eager implementation, see [`LocatedError::loc()`]
    #[track_caller]
    fn with_loc<C, F>(self, f: F) -> Result<T, anyhow::Error>
    where
        C: Display + Send + Sync + 'static,
        F: FnOnce() -> C
    {
        match self {
            Ok(ok) => Ok(ok),
            Err(_) => {
                let caller = std::panic::Location::caller();
                let loc = format!("[{}:{}:{}]", caller.file(), caller.line(), caller.column());
                self.with_context( || format!("{loc} {}", f()))
            }
        }
        
    }
}


/// Trait extending [`Option<T>`] with additional information regarding the location of the error (e.g. file, line and column)
/// 
/// # Example 
/// ```should_panic
/// use anyhow::{anyhow, Result};
/// use crate::located_error::{LocatedOption, LocatedError};
/// 
/// // ---- Main runner
/// fn maybe_bytes(n: usize) -> impl Iterator<Item=u8> + 'static {
///     vec![0; n].into_iter() // This is dumb and could fail at any time...
/// }
/// 
/// fn compute(n: usize) -> Result<()> {
///     let vec = maybe_bytes(n).next() // Oh, that's dangerous.
///         .loc("Unexpected empty vector while retrieving bytes")?; 
///     /* ---- expensive computation ensues ---- */
///     Ok(())
/// }
/// 
/// // ---- Main
/// fn main() -> Result<()> {
///     let n = 0;
///     let file = compute(n)
///         .with_loc(||format!("While attempting to run computations with n={n}"))?;
///     Ok(())
/// }
/// ```
/// ## This should output the following lines 
/// ```Text
/// > Error: [src/lib.rs:21:10] While attempting to run computations with n=0
///
/// > Caused by:
/// >     [src/lib.rs:12:17] Unexpected empty vector while retrieving bytes
/// ```
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
            Some(ok) => Ok(ok),
            None     => {
                let loc = loc_caller!(Location::caller());
                self.context(format!("{loc} {context}"))
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
            Some(ok) => Ok(ok),
            None     => {
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
            let mut chain = err.chain();
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