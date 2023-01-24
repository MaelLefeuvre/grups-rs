# GRUPS-rs : A high-performance rust port and update of [grups](https://github.com/sameoldmike/grups).  

![Build](https://github.com/MaelLefeuvre/grups/workflows/Build/badge.svg)

## Introduction

GRUPS (*Get Relatedness Using Pedigree Simulations*) is an ancient-DNA genetic relatedness kinship estimation software, revolving around pedigree simulations using existing genotype callsets (Typically, the [1000G-phase3 dataset](https://www.internationalgenome.org/category/phase-3/)), and based on the methods developed in [Martin D., et al. (2017)](https://doi.org/10.1111/mec.14188), 

GRUPS-rs currently supports the simulation of user-defined pedigree in thousand of replicates, to compute the expected genetic distances of two individuals under various relatedness scenarios (See section: [Defining custom pedigrees](#Defining-custom-pedigrees)).

Modern human contamination, sequencing errors and allele-fixation parameters can also be introduced within the simulations to account for the biases they may introduce in real-life datasets.

## Software Dependencies
### 1. cargo  
This project is written in [Rust](https://www.rust-lang.org/), and thus requires [cargo](https://crates.io/) for source compilation.  

To install cargo:
```Bash
user@desktop:~$ curl --proto '=https' --tlsv1.2 https://sh.rustup.rs -sSf | sh
```
    
See `Cargo.toml` for a complete list of crate dependencies (these are automatically downloaded by cargo when setting up compilation)

### 2. libsvm

Rust uses ordinally partitionned Sequence Vector Machines (SVMOPs) to find the best separation between expected distributions of relatedness and thus requires `libsvm-dev` to be installed on your workspace:
```Bash
user@desktop:~$ sudo apt-get install libsvm-dev
```

The original publication of libsvm can be found [here](https://doi.org/10.1145/1961189.1961199). More specifically, see [this paper](https://doi.org/10.1007/3-540-44795-4_13) for a more detailled explanation regarding SVMOPs.

---
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

3. Compile and install 
   ```Bash
   user@desktop:~$ RUSTFLAGS="-Ctarget-cpu=native" cargo install --path .
   ```

4. Grups should be located within `~/.cargo/bin/` and included in your PATH
    ```Bash
    user@desktop:~$ grups
    GRUPS: Get Relatedness Using Pedigree Simulations

    USAGE:
        grups [OPTIONS] <SUBCOMMAND>
    
    OPTIONS:
        -h, --help       Print help information
        -q, --quiet      Disable warnings
        -v, --verbose    Set the verbosity level (-v -vv -vvv -vvvv)
    
    SUBCOMMANDS:
        from-yaml         Run GRUPS using a previously generated .yaml config file
        fst               Convert VCF files into FST-indexes (Finite State Transcucer)
        help              Print this message or the help of the given subcommand(s)
        pedigree-sims     Run pwd-from-stdin and perform pedigree simulations in one go.
        pwd-from-stdin    Compute the average PairWise Difference (PWD) within a pileup file
   ```

## Data Dependencies:
### 1. SNP Callset

GRUPS requires an SNP-callset in the form of `.vcf` or `.vcf.gz` files to perform pedigree simulations. Any VCF file will work (see. [Caveats](#Caveats-(when-using-an-alternative-callset)), but GRUPS remains, as of now, mainly designed to work with the [1000G-phase3 dataset](https://www.internationalgenome.org/category/phase-3/).

The 1000G-phase3 dataset can be downloaded [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).

### 2. Input panel definition file

GRUPS will require an input panel definition file to distinguish (super-)populations and define samples within your SNP Callset.  

If you plan to use the `1000G-phase3` callset, a predefined panel can be downloaded [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_male_samples_v3.20130502.ALL.panel)

This file is in tab-separated format and must at least contain the following columns:
```Text
<SAMPLE-ID>    <POP-ID>    <SUPER-POP-ID>
```

### 3. Recombination Maps

GRUPS requires a genetic recombination map to simulate meiosis. For main intents and purposes, and when using the `1000G-phase3` callset, we recommend the HapMap-II-b37 map, which can be downloaded [here](https://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/)

### Caveats (when using an alternative SNP-callset)

If you plan to use an alternative SNP Callset, here are a few caveats you should keep in mind when preparing your input:

1. GRUPS will require an input panel definition file to your SNP Callset. See the [Input panel definition File](#2.-Input-panel-definition-file) section for the appropriate format, and/or check the 1000G-phase3 `.panel` file as a template [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_male_samples_v3.20130502.ALL.panel).

2. As of now, GRUPS does not calculate population allele frequencies by itself, and will rather look through the `INFO` field (column 7) of your callset for the appropriate tag. Thus, if you plan to simulate genomes using the `EUR` population, your SNP callset *must* carry a `EUR_AF` info tag, specifying the the allele frequency, at each SNP coordinate. See the [bcftools +fill-tag](https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html) plugin documentation to learn how to annotate population-specific allele frequencies.

3. As of now, GRUPS distinguishes relevant SNP coordinates within the callset using various `INFO` field annotations. If you plan to use a different SNP-callset, either ensure your `.vcf` files are correctly annotated with the following tags, or make sure all position that are not bi-allelic SNPs have been thoroughly filtered-out from your dataset before using GRUPS:
    - Poly-allelic sequence variations are distinguished (and ignored) by searching for the `MULTI_ALLELIC` tag.
    - SNPs are distinguished from other types of mutation by searching for the `VT=SNP` tag.  

4. By default, GRUPS will consider the provided SNP callset as being called on the `GRCh37` reference genome. If your callset has been generated using another reference genome, you'll have to provide the software with a fasta index file (`.fa.fai`) of your reference, using the `--fasta` argument (See the [pwd-from-stdin parameter list](#pwd-from-stdin) section for more information).

---
## Usage

### pwd-from-sdin

#### Basic example
```Bash
grups pwd-from-stdin --pileup ./tests/test-data/pileup/parents-offspring.pileup --output-dir ./test-grups
```
```
Name                 - Overlap - Sum PWD - Avg. Pwd - Avg. Phred
Ind0-Ind1            - 107     - 20      - 0.18692  - 37.57009
```

or 

```Bash
cat ./tests/test-data/pileup/parents-offspring.pileup | grups pwd-from-stdin --output-dir ./test-grups
```
```
Name                 - Overlap - Sum PWD - Avg. Pwd - Avg. Phred
Ind0-Ind1            - 107     - 20      - 0.18692  - 37.57009
```

#### Grups can also Support multiples samples, with or without self-comparison.
```Bash
grups pwd-from-stdin --pileup  ./tests/test-data/pileup/parents-offspring.pileup \
                     --samples 0-2                                               \
                     --sample-names MDH1 MDH2 MDH3                               \
                     --min-depth 2                                               \
                     --self-comparison          
```
```Text
Name                 - Overlap - Sum PWD - Avg. Pwd - Avg. Phred
MDH1-MDH1            - 16      - 0       - 0.00000  - 37.81250
MDH1-MDH2            - 1       - 0       - 0.00000  - 38.00000
MDH1-MDH3            - 4       - 1       - 0.25000  - 35.50000
MDH2-MDH2            - 5       - 0       - 0.00000  - 37.00000
MDH2-MDH3            - 2       - 2       - 1.00000  - 38.50000
MDH3-MDH3            - 31      - 0       - 0.00000  - 35.74194                             
```

### Pedigree-sims 

#### Basic example

Grups can then perform pedigree simulations to compute the most likely relationship of each individual in a single go.

```
grups pedigree-sims --pileup ./tests/test-data/pileup/parents-offspring.pileup \
                    --data-dir ./tests/test-data/vcf/binary-2FIN-1ACB-virtual/ \
                    --recomb-dir ./tests/test-data/recombination-map/          \
                    --pedigree ./tests/test-data/pedigree/tiny_pedigree.txt    \
                    --samples 0-2                                              \
                    --sample-names MDH1 MDH2 MDH3                              \
                    --reps 1000
                    --quiet
```
```Text
Name                 - Overlap - Sum PWD - Avg. Pwd - Avg. Phred
MDH1-MDH2            - 107     - 20      - 0.18692  - 37.49533
MDH1-MDH3            - 86      - 16      - 0.18605  - 36.87209
MDH2-MDH3            - 31      - 9       - 0.29032  - 37.96774
--------------------------------------------------------------
MDH1-MDH2            - First Degree         - 0.186916 - 0.180952 - 0.005964
MDH1-MDH3            - First Degree         - 0.186047 - 0.141176 - 0.044870
MDH2-MDH3            - Unrelated            - 0.290323 - 0.322581 - 0.032258
```

### Defining custom pedigrees

Defining pedigrees within grups is performed through simple definition files. See the example pedigree [here](resources/pedigrees/example_pedigree.txt)  

In essence, a pedigree in GRUPS is defined and parsed in three distinct steps, each one tied to a keyword within the definition file:

1. `INDIVIDUALS`: Define the individuals within the pedigree.
    - Individuals are then defined by a unique, line-separated id or name.
    - Ids must not contain any whitespace
    - Dashes, underscores, and other special characters are allowed, but we recommend that users stick to using alphanumeric characters.

    - Example:
      ```
      # Define individuals within the pedigree
      INDIVIDUALS
      father
      mother
      child
      ```

2. `RELATIONSHIPS`: Define the parents for each offspring.
    - Relationships are parsed by targeting the `=repro()` regular expression, through this nomenclature:
      ```
      <offspring-id>=repro(<parent1-id>, <parent2-id>)
      ```
    - Example:
      ```
      # Define offspring relationships
      RELATIONSHIPS
      child=repro(father,mother)
      ```

Note that comment lines are allowed within the definition file: any line starting with a `#` character is ignored during parsing.


3. `COMPARISONS` Define which pairwise comparisons should grups investigate to compute genetic distances.
    - Each comparison is defined by a unique, line-separated id or name (e.g. 'parents', 'siblings').
    - comparison ids can contain whitespaces, and various special characters (though we recommend sticking to alphanumeric characters and underscores).
    - Comparisons are then parsed by targeting the `=compare()` regular expression, through this nomenclature:
      ```
      <relationship-id>=compare(<individual1-id>, individual2-id>)
      ```
    - Example:
      ```
      # Define pairwise comparisons
      COMPARISONS
      unrelated=compare(father,mother)
      1st_degree=compare(child, mother)
      ```

### Fst indexation

#### Basic principles
The biggest performance bottleneck of GRUPS is I/O. Working with heavy vcf files such as the 1000G project can have its toll on performance. For extended uses of grups on a predefined dataset, FST-indexation of VCF files can in certain cases greatly increase performances, at the cost of a higher memory footprint.

```
grups fst --vcf-dir ./tests/test-data/vcf/binary-2FIN-1ACB-virtual/ --output-dir ./test-fst-index
```

#### Multithreading
FST-indexation can be quite long and resource intensive (altough it remains a one-time operation). Thus, the use of multithreading across vcf files is recommended, provided your computer is equipped with multiple cores:
```
grups fst --threads `nproc` --vcf-dir ./tests/test-data/vcf/binary-2FIN-1ACB-virtual/ --output-dir ./test-fst-index
```

If, your input VCF are BGZF-compressed, `.vcf.gz` files one can also leverage multiple decompression thread per vcf file.
```
grups fst --decompression-threads 8 --vcf-dir ./tests/test-data/vcf/binary-2FIN-1ACB-virtual/ --output-dir ./test-fst-index
```

**NB:** be aware that the combined use of `--threads` and `--decompression-threads` is ***multiplicative***. E.g.: calling `grups fst` with `--threads 3` and `--decompression-threads 4` will preempt 12 OS-threads.

#### FST pre-filtration
On top of this, FST-indexation has the added benefit of performing prefiltration of unwanted entries within the input vcf file. Most notably:
- Duplicate coordinate entries.
- Coordinates containing a `MULTI_ALLELIC` tag within the `INFO` field.
- Entries which are not of type `VT=SNP`.

#### FST population subsetting
More-over, if the user expects to use only a single pedigree and contaminating population, FST indexation can be used to filter-out unused samples from the original VCF file
```
grups fst --vcf-dir ./tests/test-data/vcf/binary-2FIN-1ACB-virtual/ --pop-subset EUR|AFR --output-dir ./test-fst-index
```

Once the indexation is completed, `.fst` and `.fst.frq` files can be used seamlessly when performing pedigree simulations. The user merely has to specify the input type using `--mode fst`. Specifying a target directory is performed in the same way, using `--data-dir`.

```Bash
grups pedigree-sims --pileup ./tests/test-data/pileup/parents-offspring.pileup \
                    --data-dir ./test-fst-index                                \
                    --recomb-dir ./tests/test-data/recombination-map/          \
                    --pedigree ./tests/test-data/pedigree/tiny_pedigree.txt    \
                    --samples 0-2                                              \
                    --sample-names MDH1 MDH2 MDH3                              \
                    --reps 1000                                                \
                    --quiet                                                    \
                    --mode fst                                                 
```
```
Name                 - Overlap - Sum PWD - Avg. Pwd - Avg. Phred
MDH1-MDH2            - 107     - 20      - 0.18692  - 37.55140
MDH1-MDH3            - 86      - 16      - 0.18605  - 36.87209
MDH2-MDH3            - 31      - 9       - 0.29032  - 37.87097
--------------------------------------------------------------
MDH1-MDH2            - First Degree         - 0.186916 - 0.153543 - 0.033373
MDH1-MDH3            - First Degree         - 0.186047 - 0.201541 - 0.015495
MDH2-MDH3            - Unrelated            - 0.290323 - 0.322581 - 0.032258
```

#### Tradeoffs.

- FST indexes are a quite compressed data form, but not as much as `.vcf.gz` files. Users should expect a ~two-fold input file size increase when going from `.vcf.gz` to `.fst`.

---
## Running GRUPS using yaml config files

When executing the `pedigree-sims` or `pwd-from-stdin` modules, GRUPS will automatically serialize your command line arguments and generate a timestamped `.yaml` configuration file containing every provided argument for the given run.

This file will be located at the root of your output directory (which can be specified using `--output-dir`).

To relaunch grups using the exact same configuration, simply run grups using the `from-yaml` module, and provide the path to the desired `.yaml` file

```Bash
grups from-yaml ./grups-output/2022-06-13T162822-pedigree-sims.yaml
```

---
## Building the GRUPS API docs.

The complete documentation for the GRUPS-rs public API can be built using `cargo`:
```Bash
cargo doc --workspace --open --no-deps --document-private-items
```
Using this command, the documentation's `html` entry-point will be located under `target/doc/pedigree_sims/index.html` and should automatically open in your default browser.

---

## Citing GRUPS-rs

If you plan to use GRUPS-rs, pleace include the original publication from [Martin D. et al (2017)](https://doi.org/10.1111/mec.14188) as a citation:
> Martin, M. D., Jay, F., Castellano, S., & Slatkin, M. (2017). Determination of genetic relatedness from low-coverage human genome sequences using pedigree simulations. Molecular ecology, 26(16), 4145–4157. https://doi.org/10.1111/mec.14188

--- 
## Parameter List 

### Common args
- `-v, --verbose`
    Set the verbosity level (-v -vv -vvv -vvvv)  

    Set the verbosity level of this program. With multiple levels -v : Info  |  -vv : Debug
    | -vvv : Trace By default, the program will still output Warnings. Use --quiet/-q to
    disable them

- `-q, --quiet`  
    Disable warnings.
            
    By default, warnings are emmited and redirected to the console, even without verbose
    mode on. Use this argument to disable this. Only errors will be displayed.


- `-h, --help`  
    Print help information
---
### pwd-from-stdin
#### Parameters
- `--pileup <PILEUP>`  
    Input pileup file.
    
    Note that in the absence of a '--pileup' argument, the program will except a data stream
    from the standard input. i.e 'samtools mpileup -B -q30 -Q30 ./samples/* | pwd_from_stdin
    [...]

- `-t, --targets <TARGETS>`  
    Provide with a list of SNP coordinates to target within the pileup.
        
    Pileup positions which are not found within the provided --targets file will be exluded
    from comparison. Accepted file formats: '.snp' : EIGENSTRAT format: space-separated
    '.vcf' : Variant Call Format: tab-separated '.tsv' :   tab-separated file with columns
    'CHR POS REF ALT' '.csv' : comma-separated file with columns 'CHR POS REF ALT' '.txt' :
    space-separated file with columns 'CHR POS REF ALT'

- `-o, --output-dir <OUTPUT_DIR>`  
    Output directory where results will be written.

    [default: grups-output]

- `-M, --min-qual <MIN_QUAL>`  
    Minimal required Base Quality (BQ) to perform comparison.
    
    Nucleotides whose BQ is lower than the provided treshold are filtered-out. Value should
    be expressed in scale PHRED-33.
    
    [default: 30]

- `-x, --min-depth <MIN_DEPTH>...`    
    Provide with the minimal sequencing depth required to perform comparison.
        
    Note that the length of the provided vector does not need to be the same as the one
    provided for individuals. When such is the case, values of --min-depth will "wrap"
    around, and start again from the beginning for the next individuals. Example:
    '--samples 0-4 --min-depth 2 3 5 ' will result in: - Ind  : 0 1 2 3 4 - Depth: 2 3 5 2 3
        
    [default: 1 1]

- `-b, --blocksize <BLOCKSIZE>`  
    Change the window size of each jackknife block

    [default: 1000000]

- `-s, --samples <SAMPLES>...`  
    0-based column index of the individuals that should be compared
    
    Argument may accept slices and/or discrete integers i.e. '--samples 1-4 7 8' will be
    parsed as: [1, 2, 3, 4, 7, 8] Note that slices are inclusive.
    
    [default: 0 1]

- `-n, --sample-names <SAMPLE_NAMES>...`  
    Provide with a list of sample names for printing.
    
    By default, individuals will be referred to using their pileup index. e.g. "Ind0, Ind1,
    Ind2, etc."
    
    The provided list must of course follow the sorted index order which was provided by
    '--samples'. Note that the length of --samples and --sample-names do not need to match.
    When such is the case, default values are used for the missing names.
    
    Example: '--sample 0-2 5 6 --sample-names ART04 ART16 ART20' - index: 0      1      2
    5     6 - name : ART04  ART16  ART20  Ind5  Ind6

- `-c, --chr <CHR> ...`  
    Restrict comparison to a given set of chromosomes.

    Argument may accept slices (inclusive) and/or discrete integers i.e. '--chr 9-11 13
    19-22 ' will be parsed as: [9,10,11,13,19,20,21,22]

- `-g, --genome <GENOME>`  
    Fasta indexed reference genome.

    Note that a '.fasta.fai' genome index file must be present at the same directory.

#### Flags 
- `-f, --filter-sites`  
    Do not perform comparison, but rather print out the pileup lines where a comparison
    could be made.
    
    Note that lines are printed as soon as there is a valid comparison between ANY pair of
    individuals. Thus, may not prove very effective when using --self-comparison.

- `-i, --ignore-dels`  
    Ignore deletions when performing comparison.
    
    By default, deletions ('*' character) will count as a selectible 'nucleotide' for
    comparison.

- `-J, --print-blocks`  
    Print jackknife blocks for each individual

- `-k, --known-variants`  
    Filter out tri-allelic sites when given a list of SNP-coordinates.
    
    Note that this arguement requires the use of --targets to provide the program with a
    list of coordinates. The given coordinate file must furthermore explicitely provide with
    the REF/ALT alleles.

- `-q, --quiet`  
    Disable warnings.
    
    By default, warnings are emmited and redirected to the console, even without verbose
    mode on. Use this argument to disable this. Only errors will be displayed.

- `-S, --self-comparison`  
    Enable self-comparison mode on individuals.
        
    Note that --min-depth>=2 is required when comparing individuals to themselves, since at
    least two alleles are required to effectively compute pairwise differences at a given
    site.

- `-w, --overwrite`  
    Overwrite existing output files.  

---
### Pedigree-sims
#### Parameters
- `-F, --data-dir <DATA_DIR>`  
    Path to input VCF genomes for founder individuals (1000G)

- `-G, --recomb-dir <RECOMB_DIR>`  
    Path to recombination genetic map data

- `-p, --panel <PANEL>`  

- `-P, --pedigree-pop <PEDIGREE_POP>`  
    Superpopulation with which to perform pedigree simulations

    [default: EUR]

- `-C, --contam-pop <CONTAM_POP>...`  
    Superpopulation with which to contaminate pedigree simulations.        
    Format: <POP> <POP> ...

    [default: AFR]

- `-N, --contam-num-ind <CONTAM_NUM_IND>...`  
    Number of random individual genomes with which to contaminate pedigree simulations.
    Format: <INT> <INT> ... 
    Default: Contaminate all pedigree individuals with a single contaminating individual.

    [default: 1]

- `-Q, --contamination-rate <CONTAMINATION_RATE>...`  
    Contamination rates (or rate ranges) for each pileup individual.
    Format: Format: <FLOAT> <FLOAT> ... or <FLOAT-FLOAT> <FLOAT-FLOAT> ...

    [default: 0 0]

- `-d, --af-downsampling-rate <AF_DOWNSAMPLING_RATE>`  
    Proportion of SNPs to keep at true frequency

    [default: 0.0]

- `-D, --snp-downsampling-rate <SNP_DOWNSAMPLING_RATE>`  
    Proportion of filtered SNP positions to include in the analysis 

    [default: 0.0]

- `-m, --maf <MAF>`  
    Minimum allele frequency within the pedigree superpopulation that is required to include
    a given SNP within the simulations    

    [default: 0.00]

- `-I, --mode <MODE>`  
    Define the expected data input type for pedigree simulations

    [default: vcf]  [possible values: vcf, fst]

- `-#, --decompression-threads <DECOMPRESSION_THREADS>`  
    Number of additional parallel decompression threads.

    Can increase performance when working with BGZF compressed vcf files.

    [default: 0]

---
### FST Index
#### Parameters 
- `-d, --vcf-dir <VCF_DIR>`  
    Path to input VCF genomes for founder individuals (1000G)


- `-o, --output-dir <OUTPUT_DIR>`  
    Output directory

- `-p, --panel <PANEL>`
    Path to input Panel Definition files

- `-P, --pop-subset <POP_SUBSET>...`  
    Population Subset    

    Subset the index by a given number of (super)-population. (e.g. EUR, AFR). Note that at
    least one pedigree population and one contaminating population are required to obtain
    valid index files.

- `-@, --threads <THREADS>` 
    Number of parallel CPU processes when performing FST-Indexation
        
    Parallelization is dispatched according to the number of separate vcf(.gz) files.  

    [default: 1]

- `-#, --decompression-threads <DECOMPRESSION_THREADS>`  
    Number of additional parallel decompression threads when decompressing (BGZF compressed
    files only).
        
    Parallelization is dispatched according to the number of replicates.   

    [default: 0]

## Contributing
---
## Obtaining code coverage metrics for GRUPS-rs
    
1. Install cargo tarpaulin  
    ```
    cargo install cargo-tarpaulin
    ```

 2. Run coverage tests  
    ```
    cargo tarpaulin --workspace --command test --out Lcov
    ```
