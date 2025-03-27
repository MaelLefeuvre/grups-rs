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

