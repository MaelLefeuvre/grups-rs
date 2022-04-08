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

    USAGE:
        pwd_from_stdin [OPTIONS]

    OPTIONS:
        -b, --blocksize <BLOCKSIZE>             [default: 1000000]
        -c, --chr <CHR>...
        -f, --filter-sites
        -g, --genome <GENOME>
        -h, --help                              Print help information
        -i, --ignore-dels
        -k, --known-variants
        -m, --min-depth <MIN_DEPTH>...          [default: 1 1]
        -M, --min-qual <MIN_QUAL>               [default: 30]
        -p, --print-blocks
        -p, --pileup <PILEUP>
        -q, --quiet
        -s, --samples <SAMPLES>...              [default: 0 1]
        -s, --sample-names <SAMPLE_NAMES>...
        -S, --self-comparison
        -t, --targets <TARGETS>
        -v, --verbose
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
