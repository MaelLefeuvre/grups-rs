use criterion::{black_box, criterion_group, criterion_main, Criterion};
use genome::{coordinate::{Position, ChrIdx, Coordinate}, snp::Allele};
use pwd_from_stdin::pileup::{Pileup, Line};


use located_error::prelude::*;
use itertools::Itertools;

/// Returns an iterator Providing the start and stop index of each 3-field window within a pileup line
/// e.g.: 1|75681|G|5|..,..|JJEJJ|3|..C DDD$
///       + --------- + --------- + ------ +
///           step 0      step 1     step 2
fn iter_window(line: &str) -> impl Iterator<Item=(usize, usize)> + '_ {
    std::iter::once(0).chain(
        line.match_indices('\t').skip(2)
            .step_by(3)
            .map(|(index, _)| index + 1)
    )
    .chain(std::iter::once(line.len()))
    .tuple_windows()
}
/// Default, old implementation.
pub fn pileup_new_split(line: &str, ignore_dels: bool) -> Result<Line> {
    let split_line: Vec<&str>    = line.split('\t').collect();
    let chromosome: ChrIdx       = split_line[0].parse()?;
    let position  : Position     = split_line[1].parse()?;
    let reference : Allele       = split_line[2].parse()?;
    //Loop along individuals
    let mut individuals: Vec<Pileup> = Vec::new();
    for i in (3..split_line.len()).step_by(3) {
        let depth : u16  = split_line[i].parse()?;
        let bases : &str = split_line[i+1];
        let scores: &str = split_line[i+2];
        individuals.push(Pileup::new(reference, depth, bases, scores, ignore_dels)?);
    }
    Ok(Line {
        coordinate: Coordinate::new(chromosome, position),
        reference,
        individuals
    })
}

pub fn pileup_new_window(line: &str, ignore_dels: bool) -> Result<Line> {
    let mut window = iter_window(line);

    let header  = window.next().map(|(start, end)| line[start..end].split('\t')
        .collect::<Vec<&str>>())
        .expect("Window should be buildable from 'start' and 'end'");
    let chromosome: ChrIdx       = header[0].parse()?;
    let position  : Position     = header[1].parse()?;
    let reference : Allele       = header[2].parse()?;
    //Loop along individuals
    let mut individuals: Vec<Pileup> = Vec::new();
    for (start, end) in window {
        let chunk = line[start..end].split('\t').collect::<Vec<&str>>();
        let depth : u16  = chunk[0].parse()?;
        let bases : &str = chunk[1];
        let scores: &str = chunk[2];
        individuals.push(Pileup::new(reference, depth, bases, scores, ignore_dels)?);
    }
    Ok(Line {
        coordinate: Coordinate::new(chromosome, position),
        reference,
        individuals
    })
}

fn bench_pileup_split(c: &mut Criterion) {
    let mut group = c.benchmark_group("pileup");
    
    let mut line = String::from("1\t752566\tG");
    for _ in 0..5000 {
        line.push_str(&format!("\t{}\t{}\t{}", 9, "..,..,..,", "JJJJJJJJJJ"));
    }
    line.push('\n');

    group.sample_size(500);

    group.bench_function("Split_ignore_dels", |b| b.iter(|| {
        pileup_new_split(
            black_box(&line), 
            black_box(true)
        ).expect("'Split_ignore_dels' bench should start at this point")
    }));

    group.bench_function("Window_ignore_dels", |b| b.iter(|| {
        pileup_new_window(
            black_box(&line), 
            black_box(true)
        ).expect("Window_ignore_dels bench should start at this point")
    }));

    group.finish();
    
}

criterion_group!(benches, bench_pileup_split);
criterion_main!(benches);