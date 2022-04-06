![Build](https://github.com/MaelLefeuvre/grupsactions/workflow/rust.yml/badge.svg)
# grups
rust port + update of [grups](https://github.com/sameoldmike/grups).  

See:  
> Martin, M. D., Jay, F., Castellano, S., & Slatkin, M. (2017). Determination of genetic relatedness from low-coverage human genome sequences using pedigree simulations. Molecular ecology, 26(16), 4145–4157. https://doi.org/10.1111/mec.14188

## Installation from source

2. Run the test-suite
```Bash
cargo test --workspace
```

3. Compile 
   ```Bash
   cargo build --release
   ```
   executables should be in ./target/release/

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
