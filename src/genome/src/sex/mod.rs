use std::{fmt::{self, Display}, str::FromStr};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
    Unknown
}

impl Sex {
    pub fn random() -> Self {
        [Self::Female, Self::Male][fastrand::bool() as usize]
    }
}

impl FromStr for Sex {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        Ok(match s.to_lowercase().as_str() {
            "male" | "1"   => Self::Male,
            "female" | "2" => Self::Female,
            _              => Self::Unknown,
        })
    }
}

impl Display for Sex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", match self {
            Self::Female => "female",
            Self::Male   => "male",
            Self::Unknown => "unknown"
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display() {
        assert_eq!(format!("{}", Sex::Female), "female");
        assert_eq!(format!("{}", Sex::Male), "male");
        assert_eq!(format!("{}", Sex::Unknown), "unknown");
    }

    #[test]
    fn from_str() {
        assert_eq!(Sex::from_str("FEMALE"), Ok(Sex::Female));
        assert_eq!(Sex::from_str("2"), Ok(Sex::Female));
        assert_eq!(Sex::from_str("MALE"), Ok(Sex::Male));
        assert_eq!(Sex::from_str("1"), Ok(Sex::Male));
        assert_eq!(Sex::from_str(""), Ok(Sex::Unknown));
        assert_eq!(Sex::from_str("-"), Ok(Sex::Unknown));
        assert_eq!(Sex::from_str("?"), Ok(Sex::Unknown));
        assert_eq!(Sex::from_str("-9"), Ok(Sex::Unknown));
    }
}