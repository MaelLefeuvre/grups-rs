use std::{
    path::Path,
    io::{BufRead, BufReader},
    fs::File,
};

use genome::Genome;
use crate::pedigree::Pedigree;

pub fn pedigree_parser<'a> (path: &'a Path, genome: &'a Genome) -> std::io::Result<Pedigree> {
    #[derive(Debug)]
    enum ParseMode {Individuals, Relationships, Comparisons}
    let mut parse_mode = None;
    let mut pedigree = Pedigree::new();
    let reader = BufReader::new(File::open(path)?);
    for line in reader.lines() {
        let line= line?;
        let line: Vec<&str> = line.split('#').collect();
        match line[0].chars().next() { // Skip comments and empty lines.
            Some('#') | None => continue,
            Some(_)          => (),
        }

        let mut change_parse_mode = |mode| {parse_mode = Some(mode); true};
        if let true = match line[0] {
            "INDIVIDUALS"   => change_parse_mode(ParseMode::Individuals),
            "RELATIONSHIPS" => change_parse_mode(ParseMode::Relationships),
            "COMPARISONS"   => change_parse_mode(ParseMode::Comparisons),
            _               => (false)
        }{continue}

        match parse_mode {
            Some(ParseMode::Individuals)   => {
                let label = line[0].to_string();
                pedigree.add_individual(&label, None, genome.clone())?;
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

fn parse_pedline (line: Vec<& str>, regex: &str) -> std::io::Result<(String, String, String)> {
    use std::io::ErrorKind::InvalidData;
    let mut temp=line[0].trim()
        .strip_suffix(')')
        .ok_or(InvalidData)?
        .split(regex)
        .map(|s| s.to_string());

    let ind = temp.next().ok_or(InvalidData)?;

    let parents: Vec<String> = temp.next()
        .ok_or(InvalidData)?
        .split(',')
        .map(|s| s.to_string())
        .collect();
    Ok((ind, parents[0].to_owned(), parents[1].to_owned()))
}