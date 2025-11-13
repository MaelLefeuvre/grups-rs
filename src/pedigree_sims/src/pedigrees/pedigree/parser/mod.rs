use std::{
    path::Path,
    io::{self, BufRead, BufReader},
    fs::File, fmt::{Display, self},
};

use located_error::prelude::*;
use log::{debug, warn};

mod error;
use error::PedigreeBuilderError;

use super::Pedigree;

use genome::Sex;

/// Simple enum differenciating between the available
/// pedigree definition file formats
///
/// # Variants
/// - [`PedigreeFormat::Ped`]: `grups-rs` standard pedigree
///   definition file format.
/// - [`PedigreeFormat::Legacy`]: legacy `GRUPS` pedigree definition
///   file format.
/// 
enum PedigreeFormat {
    Ped,
    Legacy
}

/// Simple enum representing the three 'parsing modes' or 'steps' required to fully define a pedigree, and found
/// within the pedigree definition file.
pub enum PedigreeSection{
    Individual,
    Relationship,
    Comparison,
    Content,
}

impl From<&str> for PedigreeSection {
    fn from(value: &str) -> Self {
        match value {
            "INDIVIDUALS"   => Self::Individual,
            "RELATIONSHIPS" => Self::Relationship,
            "COMPARISONS"   => Self::Comparison,
            _               => Self::Content
        }
    }
}


/// Simple pedigree definition file line. Contains the raw string
/// contents, and the corresponding line number within the file.
struct PedigreeLine{lineno: usize, contents: String}

impl PedigreeLine {
    /// Check whether or not the pedigree line descibres the start
    /// of a new section (legacy format).
    /// # Behaviour
    /// 
    /// Returns `true` if the PedigreeLine can be parsed into a
    /// [`PedigreeSection::Content`] variant.
    fn is_legacy_section(&self) -> bool {
        let is_content = matches!(
            PedigreeSection::from(self.contents.as_str()),
            PedigreeSection::Content
        );
        !is_content
    }

    /// Determine whether or not this line is pedigree comparison definition entry
    /// i.e.: does this line start with the `COMPARE` keyword
    fn is_comparison_definition(&self) -> bool {
        self.contents.starts_with("COMPARE")
    }

    /// Split PedigreeLine fields by tab and whitespace. Filter out empty fields.
    fn split(&self) -> impl Iterator<Item = &str> {
        self.contents.split(&['\t', ' ']).filter(|field| !field.is_empty())
    }

    /// Get the index of a specific [`PedFormatField`], given the provided [`PedFormatFieldOrder`]
    fn get_field(&self, field: PedFormatField, field_order: &PedFormatFieldOrder) -> Result<&str> {
        let index = field_order.find_field_index(field)
            .with_loc(|| PedigreeBuilderError::MissingField(field))?;
        
        self.split().nth(index).with_loc(||
            PedigreeBuilderError::RetrieveField(field, index)
        )
    }
}


/// Enum defining the expected fields commonly found within `.ped` of `.fam` files.
/// 
/// # Variants
/// - [`PedFormatField::FamId`]: Family Id
/// - [`PedFormatField::Iid`]  : Individual Id
/// - [`PedFormatField::Fid`]  : Parent 1 Id (Father Id)
/// - [`PedFormatField::Mid`]  : Parent 2 Id (Mother Id)
/// - [`PedFormatField::Sex`]  : Sex of the individual
/// - [`PedFormatField::Aff`]  : assigned Phenotype
/// 
/// 
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PedFormatField {
    FamId,
    Iid,
    Fid,
    Mid,
    Sex,
    Aff
}

#[derive(Debug, Error)]
#[error("Failed to parse field {0}")]
pub struct PedFormatFieldError(String);

impl TryFrom<&str> for PedFormatField {
    type Error = PedFormatFieldError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value.to_lowercase().as_ref() {
            "famid"      => Ok(Self::FamId),
            "id" | "iid" => Ok(Self::Iid),
            "fid"        => Ok(Self::Fid),
            "mid"        => Ok(Self::Mid),
            "sex"        => Ok(Self::Sex),
            "aff"        => Ok(Self::Aff),
            other        => Err(PedFormatFieldError(other.to_string()))
        }
    }
}

impl Display for PedFormatField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", match self {
            Self::FamId           => "Family Id",
            Self::Iid             => "Individual Id",
            Self::Fid | Self::Mid => "Parent Id",
            Self::Sex             => "Sex",
            Self::Aff             => "Phenotype"
        })
    }
}

/// Enum defining the expected field order for a standard grups-rs pedigree definition file. This generally
/// Mirrors the order commonly found in `.ped` and/or `.fam` files.
/// 
/// # Variants
/// - [`PedFormatFieldOrder::Three`]: Only three fields found, in that order:
///   1. [`PedFormatField::Iid`]: Individual Id
///   2. [`PedFormatField::Fid`]: Father Id
///   3. [`PedFormatField::Mid`]: Mother Id
///
/// - [`PedFormatFieldOrder::Four`]: Only four fields found, in that order (.ped file):
///   1. [`PedFormatField::Iid`]: Individual Id  
///   2. [`PedFormatField::Fid`]: Father Id 
///   3. [`PedFormatField::Mid`]: Mother Id
///   4. [`PedFormatField::Sex`]: Sex of the individual (ignored)
///
/// - [`PedFormatFieldOrder::Six`]: All six fields found, in that order (.fam file)
///   1. [`PedFormatField::FamId`]: Family Id (ignored)
///   2. [`PedFormatField::Iid`]: Individual Id
///   3. [`PedFormatField::Fid`]: Father Id
///   4. [`PedFormatField::Mid`]: Mother Id
///   5. [`PedFormatField::Sex`]: Sex of the individual (ignored)
///   6. [`PedFormatField::Aff`]: Individual's assigned phenotype (ignored)
///
/// - [`PedFormatFieldOrder::Custom`]: Any number of fields. Fields can be given in any order, provided the file starts with a header.
/// 
/// --- 
/// 
/// The expected header labels are given below (case-insensitive):  
/// 
///   1.[`PedFormatField::FamId`]: famid (always ignored)
///   2. [`PedFormatField::Iid`]: iid | id (required)
///   3. [`PedFormatField::Fid`]: fid (required)
///   4. [`PedFormatField::Mid`]: mid (required)
///   5. [`PedFormatField::Sex`]: sex (this field is always ignored)
///   6. [`PedFormatField::Aff`]: aff (this field is always ignored)
enum PedFormatFieldOrder {
    Three([PedFormatField; 3]),
    Four([PedFormatField; 4]),
    Six([PedFormatField; 6]),
    Custom(Vec<PedFormatField>)
}

impl PedFormatFieldOrder {
    pub fn default<const N: usize>() -> Self {
        use PedFormatField as Field;
        match N {
            3 => Self::Three([Field::Iid, Field::Fid, Field::Mid]),
            4 => Self::Four([Field::Iid, Field::Fid, Field::Mid, Field::Sex]),
            6 => Self::Six([Field::FamId, Field::Iid, Field::Fid, Field::Mid, Field::Sex, Field::Aff]),
            _ => panic!("Invalid PedFormat field count"),
        }
    }

    pub fn find_field_index(&self, field: PedFormatField) -> Option<usize> {
        let find_field = |vec: &[PedFormatField] | {
            vec.iter().position(|f| *f == field)
        };
        match self {
            Self::Three(v)  => find_field(v),
            Self::Four(v)   => find_field(v),
            Self::Six(v)    => find_field(v),
            Self::Custom(v) => find_field(v),
        }
    }
}

#[derive(Debug, Error)]
#[error("Invalid field in pedigree definition file header: '{0}'")]
struct PedFormatFieldOrderError(String);

impl TryFrom<&[&str]> for PedFormatFieldOrder {
    type Error = PedFormatFieldOrderError;

    fn try_from(value: &[&str]) -> Result<Self, Self::Error> {
        let format_fields = value.iter().flat_map(|field| PedFormatField::try_from(*field)).collect::<Vec<_>>();
        Ok(Self::Custom(format_fields))
    }
}
/// Parse a pedigree definition file and return a `Pedigree` struct.
/// # Fields:
/// - `format`: Specify the appropriate pedigree definition file format
///   for this file.
pub struct PedigreeBuilder {
    lines: Vec<PedigreeLine>,
    format: PedigreeFormat,
}

impl PedigreeBuilder {
    /// Read through the contents of a pedigree definition file, and return a
    /// `PedigreeBuilder`
    /// 
    /// # Arguments:
    /// - `path`: Path leading to the input pedigree definition file.
    ///
    /// # Behaviour
    /// This struct internaly parses, cleans-up and stores every valid line
    /// found within the file as a [`PedigreeLine`] struct.
    /// 
    /// The appropriate [`PedigreeFormat`] is determined by looking through every line
    /// and checking whether any of them matches a [`PedigreeSection`] variant other
    /// than [`PedigreeSection::Content`].
    /// 
    /// Empty lines, comment lines, and inline comments are ignored.
    /// 
    pub fn new(path: &Path) -> Result<Self> {
        let loc_msg = || format!("While attempting to parse {}", path.display());
        // ---- ressource acquisition
        let file = File::open(path)
            .map_err(PedigreeBuilderError::OpenFile)
            .with_loc(loc_msg)?;
        let reader    = BufReader::new(file);
        let mut lines = vec![];
        for (i, line) in reader.lines().enumerate() {
            let line = line
                .map_err(|source| PedigreeBuilderError::IoError{source, lineno: i+1})
                .with_loc(loc_msg)?;

            // ---- Remove comments and empty lines
            if line.starts_with('#') || line.is_empty() { continue }

            // ---- Remove inline comments 
            let line = match line.split_once('#') {
                Some((contents, _)) => contents.to_owned(),
                None                => line
            };
            lines.push(PedigreeLine{lineno: i+1, contents: line});
        }

        // ---- Decide on the appropriate parsing  by checking the existence
        // of legacy 'section' labels (i.e. 'INDIVIDUALS', 'RELATIONSHIPS', 'COMPARISONS')
        let format = match lines.iter().any(PedigreeLine::is_legacy_section) {
            true  => PedigreeFormat::Legacy,
            false => PedigreeFormat::Ped
        };

        Ok(Self{lines, format})
    }

    /// Public build method to generate a [`Pedigree`] from this builder.
    pub fn build(&self) -> Result<Pedigree> {
        match self.format {
            PedigreeFormat::Legacy => self.build_legacy_format(),
            PedigreeFormat::Ped    => self.build_ped_format(),
        }
    }
    
    /// Build a [`Pedigree`] from a template pedigree file following the 
    /// [`PedigreeFormat::Ped`] file format of `grups-rs`
    fn build_ped_format(&self) -> Result<Pedigree> {
        let first_fields: Vec<&str> = self.lines[0].split().collect();

        // ---- Count the number of fields and assign a default field order
        let mut field_order = match first_fields.len() {
            3 => Ok(PedFormatFieldOrder::default::<3>()),
            4 => Ok(PedFormatFieldOrder::default::<4>()),
            6 => Ok(PedFormatFieldOrder::default::<6>()),
            n => Err(PedigreeBuilderError::InvalidFieldNumber(n))
        }?;
        

        // First line is considered a header if every field fails to get parsed as numeric.
        // If that is the case, then we can allow for a more flexible, custom order
        let mut skip = 0;
        if first_fields.iter().all(|field| field.parse::<usize>().is_err()) {
            field_order = PedFormatFieldOrder::try_from(first_fields.as_ref())
            .loc(PedigreeBuilderError::InvalidHeader)?;
            skip = 1;
        }


        // ---- Start new pedigree
        let mut pedigree = Pedigree::new();

        // ---- Add individuals
        let loc_msg = |ctxt: &str, i: usize| format!("{ctxt} while parsing line n°{i} in the pedigree definition file");
        for line in self.lines.iter().skip(skip).filter(|line| !line.is_comparison_definition()) {
            let i = line.lineno;
            let iid = line.get_field(PedFormatField::Iid, &field_order).with_loc(||
                loc_msg("Failed to retrieve value of Individual id", i)
            )?;
            let sex: Option<Sex> = line.get_field(PedFormatField::Sex, &field_order).map_or(None, |s| s.parse().ok());
            // If None: No specified sex field => Continue on.
            // If Some(Sex::Unknown) emit a warning.
            if sex.is_some_and(|s| s == Sex::Unknown) {
                warn!("Unknown sex for individual {iid}: '{}'", line.get_field(PedFormatField::Sex, &field_order)?);
            }
            let _iid = pedigree.add_individual(iid, None, sex);//.with_loc(||
                //PedigreeBuilderError::AddIndividual(iid.to_string(), i)
            //)?;

            //println!("Adding individual: {iid}")
        }

        // ---- Add relationships
        for line in self.lines.iter().skip(skip).filter(|line| !line.is_comparison_definition()) {
            let i = line.lineno;
            let iid = line.get_field(PedFormatField::Iid, &field_order)
                .expect("Invalid Id");
            let fid = line.get_field(PedFormatField::Fid, &field_order).with_loc(||
                loc_msg("Failed to retrieve value of fid field", i)
            )?; 
            let mid = line.get_field(PedFormatField::Mid, &field_order).with_loc(||
                loc_msg("Failed to retrieve value of mid id", i)
            )?;

            // Skip founders
            if fid == "0" || mid == "0" {
                continue 
            }
            
            let [iid, fid, mid] =[iid, fid, mid].map(|label| {
                pedigree.individuals.get_ind_id(label).expect("Individual should be retrievable")
            });
            //let x = pedigree.get_ind_id(iid);
            
            pedigree.set_relationship(iid, [fid, mid]);
                //.with_loc(||loc_msg(&format!("Failed to set a valid relationship for {iid}"), i))?;
        }

        // ---- Add comparisons
        for line in self.lines.iter().skip(skip).filter(|line| line.is_comparison_definition()) {
            let i = line.lineno;
            let fields = line.split().collect::<Vec<&str>>();
            let (label, ind1, ind2) = (fields[1], fields[2], fields[3]);
            pedigree.add_comparison(label, [ind1, ind2])
                .with_loc(||loc_msg(&format!("Failed to set a valid comparison for {label}"), i))?;
        }
        
        Ok(pedigree)
    }

    /// Build a [`Pedigree`] from a template pedigree file following the 
    /// [`PedigreeFormat::Legacy`] file format of `GRUPS`. This is
    /// mainly kept to maximize backwards compatibility.
    fn build_legacy_format(&self) -> Result<Pedigree> {
        debug!("Legacy format detected...");
        let mut current_parse_mode = PedigreeSection::Content;
        let mut pedigree = Pedigree::new();
        let loc_msg = |ctxt: &str, i: usize| format!("{ctxt} while parsing line n°{i} in the pedigree definition file");
        for line in &self.lines {
            let (i, contents) = (line.lineno, &line.contents);
            // ---- Switch to the relevant ParseMode if the corresponding pattern has been found
            // and skip the current line if we've just switched to a different mode.
            if line.is_legacy_section() {
                current_parse_mode = PedigreeSection::from(line.contents.as_str());
                continue;
            }

            // ---- If this is neither a comment line, not a mode-switch,
            //      run the appropriate parsing method, according to the current ParseMode
            match current_parse_mode {
                PedigreeSection::Individual => {
                    let split = contents.split_ascii_whitespace().collect::<Vec<&str>>();
                    let iid = split.first().ok_or(anyhow!(io::ErrorKind::InvalidData))?;
                    let sex: Option<Sex> = split.get(1).and_then(|s| s.parse().ok());
                    // If None: No specified sex field => Continue on.
                    // If Some(Sex::Unknown) emit a warning.
                    if sex.is_some_and(|s| s == Sex::Unknown) {
                        warn!("Unknown sex for individual {iid}: '{}'",
                            split.get(1).ok_or(anyhow!(io::ErrorKind::InvalidData))?
                        );
                    }
                    pedigree.add_individual(iid, None, sex);
                        //.with_loc(|| PedigreeBuilderError::AddIndividual(contents.to_string(), i))?;
                },
                PedigreeSection::Relationship => {
                let (offspring, parent1, parent2) = Self::parse_legacy_pedline(contents, "=repro(")
                    .with_loc(||loc_msg("Failed to parse a valid relationship", i))?;
                let [offspring, parent1, parent2] = [offspring, parent1, parent2].map(|label| pedigree.individuals.get_ind_id(&label).expect("Individual should be retrievable"));
                pedigree.set_relationship(offspring, [parent1, parent2]);
                    //.with_loc(||loc_msg(&format!("Failed to set a valid relationship for {offspring}"), i))?;
                },
                PedigreeSection::Comparison => {
                    let (label, ind1, ind2) = Self::parse_legacy_pedline(contents, "=compare(")
                        .with_loc(||loc_msg("Failed to parse a valid comparison", i))?;
                    pedigree.add_comparison(&label, [&ind1, &ind2])
                        .with_loc(||loc_msg(&format!("Failed to set a valid comparison for {label}"), i))?;

                },
                PedigreeSection::Content => ()
            }
        }
        Ok(pedigree)
    }


    /// Parse a legacy pedigree definition line, using a given regular expression
    /// # Arguments:
    /// - `line`: buffer of the current pedigree definition line
    /// - `regex`: Optional regular expression used for splitting keys from value(s).
    /// 
    /// # Errors:
    ///  - return `std::io::Error::InvalidData` when parsing has failed.
    /// 
    /// # @ TODO:
    /// - return error if `values.len()` > 2
    fn parse_legacy_pedline (line: & str, regex: &str) -> io::Result<(String, String, String)> {
        use io::ErrorKind::InvalidData;

        // ---- Split line across the provided regex and create a <String> iterator.
        let mut temp=line.trim()
            .strip_suffix(')') 
            .ok_or(InvalidData)?
            .split(regex)
            .map(ToString::to_string);

        // ---- Extract the individual/comparison label.
        let key = temp.next().ok_or(InvalidData)?;

        // ---- Split values across ',' and extract the parents / compared individuals.
        let values: Vec<String> = temp.next()
            .ok_or(InvalidData)?
            .split(',')
            .map(ToString::to_string)
            .collect();
        Ok((key, values[0].clone(), values[1].clone()))
    }

}

