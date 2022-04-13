# grups

![Build](https://github.com/MaelLefeuvre/grups/workflows/Build/badge.svg)

rust port & update of [grups](https://github.com/sameoldmike/grups).  

See:  
> Martin, M. D., Jay, F., Castellano, S., & Slatkin, M. (2017). Determination of genetic relatedness from low-coverage human genome sequences using pedigree simulations. Molecular ecology, 26(16), 4145–4157. https://doi.org/10.1111/mec.14188

## Dependencies
### 1. cargo  
This project is written in [Rust](https://www.rust-lang.org/), and thus requires [cargo](https://crates.io/) for source compilation.  

To install cargo:
```Bash
user@desktop:~$ curl --proto '=https' --tlsv1.2 https://sh.rustup.rs -sSf | sh
```
    
See `Cargo.toml` for a complete list of crate dependencies (these are automatically downloaded by cargo when setting up compilation)

## Installation from source

1. Clone this repository
   ```Bash
   user@desktop:~$ git clone git@github.com:MaelLefeuvre/grups.git
   ```

2. Run the test-suite from the project's root
    ```Bash
    user@desktop:~$ cd grups
    user@desktop:~$ cargo test --workspace
    ```

3. Compile 
   ```Bash
   user@desktop:~$ cargo build --release
   ```

4. Executables should be located within ./target/release/
    ```Bash
    user@desktop:~$ ./target/release/pwd_from_stdin --help
   pwd_from_stdin 
   Compute average PairWise Differences (PWD) within a pileup file
   
   USAGE:
      pwd_from_stdin [OPTIONS]

   OPTIONS:
    -b, --blocksize <BLOCKSIZE>
            Change the size of our jackknife blocks [default: 1000000]
    -c, --chr <CHR>...
            Restrict comparison to a given set of chromosomes
    -f, --filter-sites
            Do not perform comparison, but rather print out the pileup lines where a comparison
            could be made
    -g, --genome <GENOME>
            Fasta indexed reference genome
    -h, --help
            Print help information
    -i, --ignore-dels
            Ignore deletions when performing comparison
    -k, --known-variants
            Filter out tri-allelic sites when given a list of SNP-coordinates
    -m, --min-depth <MIN_DEPTH>...
            Provide with the minimal sequencing depth required to perform comparison [default: 1 1]
    -M, --min-qual <MIN_QUAL>
            Minimal required Base Quality (BQ) to perform comparison [default: 30]
    -o, --output-dir <OUTPUT_DIR>
            [default: grups-output]
    -o, --overwrite
            Enable output-file overwriting.
    -p, --print-blocks
            Print jackknife blocks for each individual
    -p, --pileup <PILEUP>
            Input pileup file
    -q, --quiet
            Disable warnings
    -s, --samples <SAMPLES>...
            0-based column index of the individuals that should be compared [default: 0 1]
    -s, --sample-names <SAMPLE_NAMES>...
            Provide with a list of sample names for printing
    -S, --self-comparison
            Enable self-comparison mode on individuals
    -t, --targets <TARGETS>
            Provide with a list of SNP coordinates to target within the pileup
    -v, --verbose
            Set the verbosity level (-v -vv -vvv -vvvv)
   ```
## Usage

- Basic example
    ```Bash
    ./pwd_from_stdin --pileup ./test-data/All_samples_RBq30Q30_1240k.mpileup 
    ```

    or 

    ```Bash
    cat ./test-data/All_samples_RBq30Q30_1240K.mpileup | pwd_from_stdin
    ```

- Support for multiple samples, with or without self-comparison.
    ```Bash
    ./pwd_from_stdin --pileup  ./test-data/All_samples_RBq30Q30_1240.mpileup \
                     --targets ./test-data/v50.0_1240K_public.txt            \
                     --samples 0 1 2                                         \
                     --sample-names ART16 ART21 ART24                        \
                     --min-depth 2                                           \
                     --self-comparison
    ```
