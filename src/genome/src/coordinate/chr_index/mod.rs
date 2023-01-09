use std::{str::FromStr, fmt::{Display, Formatter, self}, ops::{Deref, Add}, hash::{Hash, Hasher}, cmp::Ordering};

mod error;
use error::ChrIdxError;

use located_error::prelude::*;

#[derive(Debug, Clone, Copy)]
pub struct ChrIdx(pub u8);

impl FromStr for ChrIdx {
    
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let inner = s.replace("chr", "")
            .parse::<u8>()
            .loc(ChrIdxError(s.to_owned()))?;
        Ok(Self(inner))
    }
}

impl Display for ChrIdx {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl ChrIdx {
    pub fn into_inner(&self) -> u8 {
        self.0
    }
}

impl From<u8> for ChrIdx {
    fn from(value: u8) -> Self {
        Self(value)
    }
}

impl Into<u8> for ChrIdx {
    fn into(self) -> u8 {
        self.0
    }
}

impl AsRef<u8> for ChrIdx {
    fn as_ref(&self) -> &u8 {
        &self.0
    }
}

impl Deref for ChrIdx {
    type Target = u8;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Add for ChrIdx {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl Hash for ChrIdx {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl Ord for ChrIdx {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for ChrIdx {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq<Self> for ChrIdx {
    fn eq(&self, other: &Self) -> bool { 
        self.0 == other.0
    }
}

impl Eq for ChrIdx {}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display() {
        let pos = 22;
        let want = format!("{pos:_^12}");
        let got  = format!("{:_^12}", ChrIdx(pos));
        println!("want: {want}\ngot : {got}");
        assert_eq!(want, got);
    }
}