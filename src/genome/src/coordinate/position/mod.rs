use std::{hash::{Hash, Hasher}, fmt::{self, Display, Formatter}, str::FromStr, cmp::Ordering};
mod error;
use error::ParsePositionError;

#[derive(Debug, Clone, Copy)]
pub struct Position(pub u32);

impl Display for Position {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl FromStr for Position {
    type Err = ParsePositionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(s.parse::<u32>()?))
    }
}

impl From<u32> for Position {
    fn from(value: u32) -> Self {
        Self(value)
    }
}

impl Into<u32> for Position {
    fn into(self) -> u32 {
        self.0
    }
}

impl Hash for Position {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl Ord for Position {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for Position {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq<Self> for Position {
    fn eq(&self, other: &Self) -> bool { 
        self.0 == other.0
    }
}

impl Eq for Position {}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display() {
        let pos = 139543;
        let want = format!("{pos:_^12}");
        let got = format!("{:_^12}", Position(pos));
        println!("want: {want}\ngot : {got}");
        assert_eq!(want, got);
    }
}
