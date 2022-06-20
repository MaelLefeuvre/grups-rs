use std::{
    path::Path,
    io::{BufRead, BufReader},
    fs::File,
};

use crate::pedigrees::Pedigree;

/// Simple enum representing the three 'parsing modes' or 'steps' required to fully define a pedigree, and found
/// within the pedigree definition file.
enum ParseMode {Individuals, Relationships, Comparisons}

/// Parse a pedigree definition file and return a `Pedigree` struct.
/// # Arguments:
/// - `path`: Path leading to the input pedigree definition file.
pub fn pedigree_parser(path: &Path) -> std::io::Result<Pedigree> {
    // ---- Ressource acquisition
    let mut parse_mode = None;
    let mut pedigree = Pedigree::new();
    let reader = BufReader::new(File::open(path)?);

    // ---- Parse the pedigree definition file.
    for line in reader.lines() {
        let line= line?;
        let line: Vec<&str> = line.split('#').collect();
        match line[0].chars().next() { // Skip comments and empty lines.
            Some('#') | None => continue,
            Some(_)          => (),
        }

        // --- Switch to the relevant ParseMode if the corresponding pattern has been found.
        let mut change_parse_mode = |mode| {parse_mode = Some(mode); true};
        if let true = match line[0] {
            "INDIVIDUALS"   => change_parse_mode(ParseMode::Individuals),
            "RELATIONSHIPS" => change_parse_mode(ParseMode::Relationships),
            "COMPARISONS"   => change_parse_mode(ParseMode::Comparisons),
            _               => (false)
        }{continue} // ---- ...And skip the current line if we've just switched to a different mode.

        // ---- If this is neither a comment line, not a mode-switch,
        //      run the appropriate parsing method, according to the current ParseMode.
        match parse_mode {
            Some(ParseMode::Individuals)   => {
                let label = line[0].to_string();
                pedigree.add_individual(&label, None)?;
            },
            Some(ParseMode::Relationships) => {
                let (offspring, parent1, parent2) = parse_pedline(line, "=repro(")?;
                pedigree.set_relationship(&offspring, (&parent1, &parent2))?;

            },
            Some(ParseMode::Comparisons)   => {
                let (label, ind1, ind2) = parse_pedline(line, "=compare(")?;
                pedigree.add_comparison(&label, (&ind1, &ind2))?;
            },
            None                           => continue
        };
    }
    Ok(pedigree)
}

/// Parse a pedigree definition line, using a given regular expression
/// # Arguments:
/// - `line`: buffer of the current pedigree definition line
/// - `regex`: Optional regular expression used for splitting keys from value(s).
/// 
/// # Errors:
///  - return `std::io::Error::InvalidData` when parsing has failed.
/// 
/// # @ TODO:
/// - return error if `values.len()` > 2
fn parse_pedline (line: Vec<& str>, regex: &str) -> std::io::Result<(String, String, String)> {
    use std::io::ErrorKind::InvalidData;

    // ---- Split line across the provided regex and create a <String> iterator.
    let mut temp=line[0].trim()
        .strip_suffix(')') 
        .ok_or(InvalidData)?
        .split(regex)
        .map(|s| s.to_string());

    // ---- Extract the individual/comparison label.
    let key = temp.next().ok_or(InvalidData)?;

    // ---- Split values across ',' and extract the parents / compared individuals.
    let values: Vec<String> = temp.next()
        .ok_or(InvalidData)?
        .split(',')
        .map(|s| s.to_string())
        .collect();
    Ok((key, values[0].to_owned(), values[1].to_owned()))
}