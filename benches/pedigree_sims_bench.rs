use criterion::{black_box, criterion_group, criterion_main, Criterion};

use rand::{prelude::{ThreadRng, SliceRandom}, thread_rng};
use nanorand::WyRand;

use fastrand::{self};

fn simulate_observed_reads(n: u8, rng: &mut ThreadRng, contam_rate: f64, contam_pop_af: f64, seq_error_rate: f64, alleles: [u8; 2]) -> Result<Vec<u8>, String> {
    use rand::Rng;

    let mut reads = Vec::with_capacity(n as usize);
    //let mut rng = rand::thread_rng();

    // ---- Simulate n pileup observations.
    for _ in 0..n {
        // ---- Simulate modern human contamination. 
        let chosen_base: u8 = match rng.gen::<f64>() < contam_rate {
            true  => match rng.gen::<f64>() < contam_pop_af {
                true  => 1,  // Becomes the alternative reference allele, if contam_rate * contam_pop_af
                false => 0,  // otherwise, pick the reference allele.
            }
            false => *alleles.choose(rng)
                .ok_or_else(|| String::from("While simulating observed reads: failed to select a random allele"))?
        };

        // ---- Simulate sequencing error rate.
        let seqerror_choices: Vec<[u8; 3]> = vec![[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        if rng.gen::<f64>() < seq_error_rate {
            let wrong_base: u8 = *seqerror_choices[chosen_base as usize].choose(rng)
                .ok_or_else(|| {String::from(
                    "While simulating observed reads: failed to select a random \
                    erroneous base when simulating sequencing error."
                )
            })?;

            reads.push(wrong_base);
        }
        else {
            reads.push(chosen_base);
        }
    }
    Ok(reads)
}


fn simulate_observed_reads_wyrand(n: u8, rng: &mut WyRand, contam_rate: f64, contam_pop_af: f64, seq_error_rate: f64, alleles: [u8; 2]) -> Result<Vec<u8>, String> {
    use nanorand::Rng;
    
    let mut reads = Vec::with_capacity(n as usize);
    //let mut rng = rand::thread_rng();

    // ---- Simulate n pileup observations.
    for _ in 0..n {
        // ---- Simulate modern human contamination. 
        let chosen_base: u8 = match rng.generate::<f64>() < contam_rate {
            true  => match rng.generate::<f64>() < contam_pop_af {
                true  => 1,  // Becomes the alternative reference allele, if contam_rate * contam_pop_af
                false => 0,  // otherwise, pick the reference allele.
            }
            false => *alleles.get(rng.generate::<bool>() as usize)
                .ok_or_else(|| String::from("While simulating observed reads: failed to select a random allele"))?
        };

        // ---- Simulate sequencing error rate.
        let seqerror_choices: Vec<[u8; 3]> = vec![[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        if rng.generate::<f64>() < seq_error_rate {
            let wrong_base: u8 = *seqerror_choices[chosen_base as usize].get(rng.generate_range(0..3) as usize)
                .ok_or_else(|| {String::from(
                    "While simulating observed reads: failed to select a random \
                    erroneous base when simulating sequencing error."
                )
            })?;

            reads.push(wrong_base);
        }
        else {
            reads.push(chosen_base);
        }
    }
    Ok(reads)
}

fn simulate_observed_reads_fastrand(n: u8, rng: &mut fastrand::Rng, contam_rate: f64, contam_pop_af: f64, seq_error_rate: f64, alleles: [u8; 2]) -> Result<Vec<u8>, String> {
    
    let mut reads = Vec::with_capacity(n as usize);
    //let mut rng = rand::thread_rng();

    // ---- Simulate n pileup observations.
    for _ in 0..n {
        // ---- Simulate modern human contamination. 
        let chosen_base: u8 = match rng.f64() < contam_rate {
            true  => match rng.f64() < contam_pop_af {
                true  => 1,  // Becomes the alternative reference allele, if contam_rate * contam_pop_af
                false => 0,  // otherwise, pick the reference allele.
            }
            false => *alleles.get(rng.usize(0..1))
                .ok_or_else(|| String::from("While simulating observed reads: failed to select a random allele"))?
        };

        // ---- Simulate sequencing error rate.
        let seqerror_choices: Vec<[u8; 3]> = vec![[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        if rng.f64() < seq_error_rate {
            let wrong_base: u8 = *seqerror_choices[chosen_base as usize].get(rng.usize(0..seqerror_choices.len()))
                .ok_or_else(|| {String::from(
                    "While simulating observed reads: failed to select a random \
                    erroneous base when simulating sequencing error."
                )
            })?;

            reads.push(wrong_base);
        }
        else {
            reads.push(chosen_base);
        }
    }
    Ok(reads)
}

fn bench_simread(c: &mut Criterion) {
    let mut group = c.benchmark_group("readsim");
    
    let mut rand = thread_rng();
    let mut wyrand = nanorand::WyRand::new();
    let mut fastrand = fastrand::Rng::new();

    
    group.sample_size(500);

    group.bench_function("ThreadRNG", |b| b.iter(|| {
        simulate_observed_reads(
            black_box(1), 
            black_box(&mut rand),
            black_box(0.5),
            black_box(0.5),
            black_box(0.5),
            black_box([0,1])
        )
    }));


    group.bench_function("WyRand", |b| b.iter(|| {
        simulate_observed_reads_wyrand(
            black_box(1), 
            black_box(&mut wyrand),
            black_box(0.5),
            black_box(0.5),
            black_box(0.5),
            black_box([0,1])
        )
    }));

    group.bench_function("fastrand", |b| b.iter(|| {
        simulate_observed_reads_fastrand(
            black_box(1), 
            black_box(&mut fastrand),
            black_box(0.0),
            black_box(0.0),
            black_box(0.0),
            black_box([0,1])
        )
    }));
    group.finish();
    
}

criterion_group!(benches, bench_simread);
criterion_main!(benches);