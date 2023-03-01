use grups_io::read::genotype_reader::FSTReader;
use genome::coordinate::Coordinate;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_search_genotypes(c: &mut Criterion) {
    let mut group = c.benchmark_group("fst_RAM");
    let fst_path = "tests/test-data/fst/binary-2FIN-1ACB-virtual/ALL.wgs.phase3_v5b_20130502_n100-2FIN-1ACB-binary.fst";
    let mut reader = FSTReader::new(fst_path).unwrap();
    let coordinate = Coordinate::new(1, 10411523);
    group.bench_function("search_fst_genotypes_ram", |b| b.iter(|| {
            reader.search_coordinate_genotypes(&coordinate)
    }));

    group.bench_function("search_fst_frequencies_ram", |b| b.iter(|| {
        reader.search_coordinate_frequencies(&coordinate)
    }));
}

fn bytes_to_string(c: &mut Criterion) {
    let mut group = c.benchmark_group("test_string");
    
    let slice = black_box(b"HG0096".as_slice());
    let string  = black_box(String::from("HG0096"));
    group.bench_function("bytes2String", |b| b.iter(|| {
            let _my_string = black_box(String::from(unsafe { std::str::from_utf8_unchecked(slice)}));
        })
    );

    group.bench_function("String2bytes", |b| b.iter(|| {
            let _my_bytes: &[u8] = black_box(string.as_bytes());
        })
    );
}

criterion_group!(benches, bench_search_genotypes, bytes_to_string);
criterion_main!(benches);