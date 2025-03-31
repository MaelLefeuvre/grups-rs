# GRUPS-rs : A high-performance rust port and update of [grups](https://github.com/sameoldmike/grups).  

[![Ubuntu](https://github.com/MaelLefeuvre/grups-rs/actions/workflows/Ubuntu.yml/badge.svg)](https://github.com/MaelLefeuvre/grups-rs/actions/workflows/Ubuntu.yml) [![MacOS](https://github.com/MaelLefeuvre/grups-rs/actions/workflows/MacOS.yml/badge.svg)](https://github.com/MaelLefeuvre/grups-rs/actions/workflows/MacOS.yml) [![DOI](https://zenodo.org/badge/477140542.svg)](https://zenodo.org/doi/10.5281/zenodo.10389506)
## Introduction

GRUPS-rs (_**G**et **R**elatedness **U**sing **P**edigree **S**imulations_) is an ancient-DNA genetic relatedness estimation software relying on pedigree simulations using existing genotype callsets - typically, the [1000G-phase3 dataset](https://www.internationalgenome.org/category/phase-3/). This software is a pure rust implementation and update of the "GRUPS" method, initially developed in [Martin, et al. (2017)](https://doi.org/10.1111/mec.14188).

GRUPS-rs currently supports the simulation of user-defined pedigree in thousand of replicates, to compute the expected genetic distances of two individuals under various relatedness scenarios (See section: [Defining custom pedigrees](#Defining-custom-pedigrees)), including complex inbreeding scenarios such as double first-cousins, three-quarters siblings, sesqui first-cousins, etc.

Modern human contamination, sequencing errors and allele-fixation rate parameters can also be introduced during simulations to account for the biases they may introduce in real-life datasets.

## Installation
### Software Dependencies

If you plan to install  from source, you'll need:
1. The cargo compiler [cargo](https://crates.io/). (version `>=1.70`).
2. The [libsvm](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) library (a version `>=3.24` is recommended)

If you plan to use the [`grups.plots`](https://github.com/MaelLefeuvre/grups.plots) companion shiny dashboard, you'll need to have [R](https://www.r-project.org/), along with the [`devtools`](https://github.com/r-lib/devtools) library. A version of `R` that is `>=4.1.2` is recommended. 

#### 1. cargo  
This project is written in [Rust](https://www.rust-lang.org/), and thus requires [cargo](https://crates.io/) for source compilation. The current minimum supported rust version is `1.66`.  

To install the latest version of cargo:
```Bash
curl --proto '=https' --tlsv1.2 https://sh.rustup.rs -sSf | sh
```

See `Cargo.toml` for a complete list of crate dependencies (these are automatically downloaded by cargo when setting up compilation)

#### 2. libsvm

Rust uses ordinally partitionned Sequence Vector Machines (SVMOPs) to find the best separation between expected distributions of relatedness and thus requires the [libsvm](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) library to be installed on your workspace. (a version `>=2.34` is recommended.)

- **On Ubuntu**:
   ```Bash
   sudo apt-get install libsvm-dev
   ```
- **On Mac OS**:
   ```bash
   brew install libsvm
   ```

The original publication of libsvm can be found [here](https://doi.org/10.1145/1961189.1961199). More specifically, see [this paper](https://doi.org/10.1007/3-540-44795-4_13) for a more detailled explanation regarding SVMOPs.

---

### Installating `grups-rs` from source

1. Clone this repository (Note the `--recursive` flag is required if you wish to download the appropriate version of `grups.plots`) 
   ```Bash
   git clone --recursive https://github.com/MaelLefeuvre/grups-rs.git
   cd grups-rs
   ```
2. (**MacOS only**) Manually target the `include` and `lib` directory of libsvm, by exporting the `LIBSVM_INCLUDE` and `LIBSVM_LIBRARY` environment variables
   ```bash
   export LIBSVM_INCLUDE="$(brew --prefix libsvm)/include"
   export LIBSVM_LIBRARY="$(brew --prefix libsvm)/lib"
   ```
   
3. (Optional) Run the test-suite from the project's root
    ```Bash
    cargo test --workspace
    ```

4. Compile and install 
   ```Bash
   RUSTFLAGS="-Ctarget-cpu=native" cargo install --path .
   ```

5. `grups-rs` should be located within `~/.cargo/bin/` and included in your PATH
    ```bash
    grups-rs
    ```
    <pre class="ansi2html-content">
    <span style="color: mediumseagreen">grups-rs</span> 0.2.0
    Maël Lefeuvre &lt;mael.lefeuvre@mnhn.fr&gt;
    GRUPS-rs: Get Relatedness Using Pedigree Simulations
    
    <span style="color: gold">USAGE:</span>
        grups-rs [OPTIONS] &lt;SUBCOMMAND&gt;
    
    <span style="color: gold">OPTIONS</span>:
        <span style="color: mediumseagreen">-h</span>, <span style="color: mediumseagreen">--help</span>       Print help information
        <span style="color: mediumseagreen">-q</span>, <span style="color: mediumseagreen">--quiet</span>      Disable warnings
        <span style="color: mediumseagreen">-v</span>, <span style="color: mediumseagreen">--verbose</span>    Set the verbosity level (-v -vv -vvv -vvvv)
        <span style="color: mediumseagreen">-V</span>, <span style="color: mediumseagreen">--version</span>    Print version information
    
    <span style="color: gold">SUBCOMMANDS</span>:
        <span style="color: mediumseagreen">cite</span>              
        <span style="color: mediumseagreen">from-yaml</span>         Run GRUPS using a previously generated .yaml config file
        <span style="color: mediumseagreen">fst</span>               Convert VCF files into FST-indexes (Finite State Transcucer)
        <span style="color: mediumseagreen">help</span>              Print this message or the help of the given subcommand(s)
        <span style="color: mediumseagreen">pedigree-sims</span>     Run pwd-from-stdin and perform pedigree simulations in one go
        <span style="color: mediumseagreen">pwd-from-stdin</span>    Compute the average `PairWise Difference` (PWD) within a pileup file
    
    </pre>

---

### Installing the `grups.plots` companion interface

Installing `grups.plots` can be installed using a single command, provided you already have `R` and the `devtools` library preinstalled.

1. (Optional) Install the devtools library, if it is not preinstalled
   ```bash
   R --slave -e 'if (!require("devtools")) install.packages("devtools", repos = c(CRAN = "https://cloud.r-project.org"))'
   ```

2. Install `grups.plots`
   ```bash
   R --slave -e 'devtools::install("./grups.plots")'
   ```

See the [documentation of `grups.plot`](./grups.plots/README.md) to learn how to launch this companion package

---

## Getting help

There are multiple ways to obtain information about the different submodules of `grups-rs`
- Typing `grups-rs [submodule] -h` will display a *short* description of every available command line argument for a given command 
- Typing `grups-rs [submodule] --help` or `grups-rs [submodule] help` will display a *verbose* description of every available command line argument.

See the section [Parameter List](#parameter-list), for a detailled description of every module's command line interface.

---

## Data Dependencies:
`grups-rs` requires 4 types of additional input files and datasets to function:
1. A reference dataset of phased, modern human genotypes, such as the [1000g-phase3 dataset](https://doi.org/10.1038/nature15393) (See subsection [SNP callset](#1-snp-callset)).
2. An input panel definition file, specifying the name, population, and super-population of each sample found within the reference dataset (See subsection [Input panel definition file](#2-input-panel-definition-file)). 
3. A set of per-chromosome recombination map files, such as the [phaseII HapMap dataset](https://doi.org/10.1038/nature06258) (See subsection [Recombination Maps](#3-recombination-maps)).
4. A user-defined pedigree definition file. A set of pre-defined files can be found in the `resources/pedigrees` directory of this repository. See section [Defining custom pedigrees](#defining-custom-pedigrees), for a detailled explanation on how to create custom template pedigrees.

### 1. SNP Callset
GRUPS-rs requires an SNP-callset in the form of `.vcf` or `.vcf.gz` files to perform pedigree simulations. For most intents and purposes, the [1000g-phase3 dataset](https://www.internationalgenome.org/category/phase-3/) may provide with a good start, but any dataset of input VCF files will work, provided they carry phased diploid genotypes, and contain the appropriate required tags within the `INFO` field  - see [Caveats](#Caveats-(when-using-an-alternative-callset)).

The 1000g-phase3 dataset can be downloaded [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).

### 2. Input panel definition file
GRUPS-rs will require an input panel definition file to distinguish (super-)populations and define samples within your SNP Callset. If you plan to use the `1000G-phase3` callset, a predefined panel can be previewed and downloaded [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_male_samples_v3.20130502.ALL.panel)

This file must be unheaded, tab-separated, and should at least contain the following columns: `<SAMPLE-ID>    <POP-ID>    <SUPER-POP-ID>`

### 3. Recombination Maps
GRUPS-rs requires a genetic recombination map to simulate meiosis. For main intents and purposes, and when using the `1000g-phase3` callset, we recommend the HapMap-II-b37 map, which can be downloaded [here](https://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/)

### Caveats (when using an alternative SNP-callset)
If you plan to use an alternative SNP Callset, here are a few caveats you should keep in mind when preparing your input:

1. GRUPS-rs will require an input panel definition file to your SNP Callset. See the [Input panel definition file](#2-input-panel-definition-file) section for the appropriate format, and/or check the 1000G-phase3 `.panel` file as a template [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_male_samples_v3.20130502.ALL.panel).

2. As of now, GRUPS-rs does not calculate population allele frequencies by default, and will rather look through the `INFO` field (column 7) of your callset for the appropriate tag. Thus, if you plan to simulate genomes using the `EUR` population, each entry within your VCF files should carry a `EUR_AF` info tag, specifying the alternative allele frequency for that population. The [bcftools +fill-tags](https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html) plugin documentation to may be of help, should you wish annotate population-specific allele frequencies on your dataset.

3. As of now, GRUPS-rs distinguishes relevant SNP coordinates within the callset using various `INFO` field annotations. If you plan to use a different SNP-callset, either ensure your `.vcf` files are correctly annotated with the following tags, or make sure to filter all position that are not bi-allelic SNPs have been thoroughly filtered-out from your dataset before using GRUPS-rs:
    - Poly-allelic sequence variations are distinguished (and ignored) by searching for the `MULTI_ALLELIC` tag.
    - SNPs are distinguished from other types of mutation by searching for the `VT=SNP` tag.

4. By default, GRUPS-rs will consider the provided SNP callset as being called on the `GRCh37` reference genome. If your callset has been generated using another reference genome, we recommend to provide the software with a fasta index file (`.fa.fai`) of your reference, using the [`--genome`](#g--genome) argument (See the [pwd-from-stdin parameter list](#pwd-from-stdin) section for more information).

---

## Usage

The following section provides an quick explanation of each module of `grups-rs`, skip to section [Quick Start](#quick-start) if you wish to jump in straight ahead with a more pragmatic example.

### The `pedigree-sims` module

`pedigree-sims` is the main module of `grups-rs`, and will both compute the observed pairwise mismatch rates for your samples, and perform pedigree simulations to estimate the most likely relationship between each pair. 

**A basic example, using provided dummy test files:**
  ```bash
  grups-rs pedigree-sims --pileup ./tests/test-data/pileup/parents-offspring.pileup \
  --data-dir ./tests/test-data/vcf/binary-2FIN-1ACB-virtual/ \
  --recomb-dir ./tests/test-data/recombination-map/          \
  --pedigree ./tests/test-data/pedigree/tiny_pedigree.txt    \
  --samples 0-2                                              \
  --sample-names MDH1 MDH2 MDH3                              \
  --output-dir pedigree-sims-output                          \
  --reps 1000                                                \
  --quiet
  ```
  This should create a single output directory called `grups-output`. More specifically, the `pedigree-sims` module will generate all of the output files generated by `pwd-from-stdin`, as well as:
  1. an additional [`.result`](#result-file) file, containing summary statistics and results for all pedigree simulations.
  2. a set of `.sims` files, one for each pairwise comparison. These file contain raw simulation results for each pairwise comparison, and are located in the `simulations` subdirectory.
  
  See the section [Output Files](#output-files) for a detailled explanation of each output files.

---

### The `pwd-from-sdin` module

This module is available if you simply wish to quickly examine the pairwise mismatch rate between samples *without* performing any pedigree simulations.

**A basic example, using provided dummy test files:**
```Bash
grups-rs pwd-from-stdin --pileup  ./tests/test-data/pileup/parents-offspring.pileup \
                     --samples 0 2                                               \
                     --sample-names MDH1 MDH3                                    \
                     --min-depth 2 2                                             \
                     --self-comparison
```

This should create a single output directory called `grups-output`. More specifically, `pwd-from-stin` is expected to output:
1. A [`.pwd`](#pwd-file) file, containing summary statistics and the raw observed pairwise mismatch rate for all pairwise comparisons.
2. A set of [`.blk`](#blk-files) files, one for each pairwise comparison. These file contain raw observed pairwise mismatch rates in non-overlapping blocks ([`--blocksize`](#b--blocksize) = 1Mb, by default).

See the section [Output Files](#output-files) for a detailled explanation of each input files.

---

### The `fst` module: Encoding VCF files as a set of deterministic Finite State Acceptors

#### Basic principles
As GRUPS-rs must repeatly fetch genotypes from randomly sampled individuals during pedigree simsulations, I/O operations can become a strong performance bottleneck when running this software. Working with heavy `.vcf.gz` files such as the 1000g project can thus have its toll on performance. For extended uses of `grups-rs` on a predefined dataset, indexing your `.vcf` files as a set of FSA can in many cases greatly increase performances.

This operation can be done with the `fst` module of `grups-rs`.

On top of this, FST-indexation has the added benefit of performing prefiltration of unwanted entries within the input `.vcf`` file. Most notably:
- Duplicate coordinate entries.
- Coordinates containing a `MULTI_ALLELIC` tag within the `INFO` field.
- Entries which are not of type `VT=SNP`.

Furthermore, the `fst` module can also be useful to filter out individuals from unwanted population entries, as well as (re-)computing population allele frequencies (see section [Performing population subsets with the `fst` module](#performing-population-subsets-with-the-fst-module)).

```
grups-rs fst --vcf-dir ./tests/test-data/vcf/binary-2FIN-1ACB-virtual/ --output-dir ./test-fst-index
```

In this example, `grups-rs` index any `.vcf[.gz]` file found within the provided `binary-2FIN-1ACB-virtual` input directory, and output its contents within the `test-fst-index`. The expected output is a set of two finite state automaton (`.fst` and `.fst.frq`), one for each discovered input `.vcf[.gz]` file:

- `.fst` files indexes the genotype information of all retained samples, at each valid genotype coordinate.
- `.fst.frq` files indexes population allele frequencies for each valid genotype coordinate.

#### Multithreading the `fst` module.
Altough it remains a one-time operation, FSA-indexation can be quite long and resource intensive (e.g.: around 40 minutes is required to encode the `ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` file of the 1000g-phase database).  Thus, the use of multithreading across `.vcf.gz` files is highly recommended, provided your computer is equipped with multiple cores.

```
grups-rs fst --threads 22 --vcf-dir ./tests/test-data/vcf/binary-2FIN-1ACB-virtual/ --output-dir ./test-fst-index
```

Note that multithreading is performed across the number of discovered input `.vcf[.gz]` files. Thus, if the directory specified by [`--vcf-dir`](#d--vcf-dir) contains 22 files, there is no point in recruiting more than 22 threads.

#### Performing population subsets with the `fst` module
If the user expects to use only a single pedigree and contaminating population, FSA indexation can be used to filter-out unused samples from the original VCF file. Furthermore, the use of the optional [`--compute-pop-afs`](#f--compute-pop-afs) flag can be useful to (re-)compute population allele frequencies.

```bash
grups-rs fst \
--vcf-dir ./tests/test-data/vcf/binary-2FIN-1ACB-virtual/  \
--output-dir ./test-fst-index                              \
--pop-subset FIN AFR                                       \
--compute-pop-afs
```

---

#### Using FSA-encoded files with `grups-rs`

Once the indexation is completed, `.fst` and `.fst.frq` files can be used seamlessly when performing pedigree simulations. The user merely has to specify the input type using the [`--mode`](#i--mode) argument. Specifying a target directory is performed in the same way, using [`--data-dir`](#f--data-dir).

```Bash
grups-rs pedigree-sims --pileup ./tests/test-data/pileup/parents-offspring.pileup \
                    --data-dir ./test-fst-index                                \
                    --recomb-dir ./tests/test-data/recombination-map/          \
                    --pedigree ./tests/test-data/pedigree/tiny_pedigree.txt    \
                    --samples 0-2                                              \
                    --sample-names MDH1 MDH2 MDH3                              \
                    --reps 1000                                                \
                    --quiet                                                    \
                    --mode fst-mmap                                                 
```

##### A note on the `--mode` argument:

FSA-encoded files can be used in one of two ways:

`--mode fst-mmap`: This mode is highly recommended if your files are located within an SSD drive, as it allows `grups-rs` to efficiently make use of memory-mapped files, to reduce its memory consumption.
`--mode fst`: This mode will load each `.fst` and `.fst.frq` file within `grups-rs` heap memory. This may make `grups-rs` faster, at the cost of a higher memory footprint. In general, this mode is only recommended if you are sure your computer has enough RAM, and your files are located within an HDD drive.

---

### The `from-yaml` module: Re-running `grups-rs` using `.yaml` configuration files

When executing the `pedigree-sims` or `pwd-from-stdin` modules, GRUPS-rs will automatically serialize your command line arguments and generate a timestamped [`.yaml`](#yaml-file) configuration file containing every provided argument for the given run.

This file will be located at the root of your output directory (which can be specified using [`--output-dir`](#o--output-dir)).

To relaunch `grups-rs` using the exact same configuration, simply run `grups-rs` using the `from-yaml` module, and provide the path to the desired [`.yaml`](#yaml-file) file

```bash
grups-rs from-yaml ./grups-output/2022-06-13T162822-pedigree-sims.yaml
```

Note that the [`--quiet`](#q--quiet) and [`--verbose`](#v--verbose) arguments found within the [`.yaml`](#yaml-file) file are ignored by the `from-yaml` module, and must be-respecified when running that command :
```bash
grups-rs from-yaml ./grups-output/2022-06-13T162822-pedigree-sims.yaml --verbose
```

---

### Quick start

In this example, our objective is to apply `grups-rs` on a set of three individuals found buried within a collective grave from Mentesh-Tepe (`MT23` `MT26` and `MT7`). Two of these samples - `MT23` and `MT26` - were previously determined to be siblings [(Garino-Vignon et al. 2023)](https://doi.org/10.1038/s42003-023-04681-w). 

In this example, we'll be starting from scratch, and assume that you don't have any of the data required to begin this analysis on your computer. Thus, our workflow can be subdivided in three steps:
1. Download and process all of the datasets and files required for this analysis (i.e: (i) [1000g-phase3](https://doi.org/10.1038/nature15393) dataset, (ii) [HapMap-phase2](https://doi.org/10.1038/nature06258) recombination map, (iii) [GRCh37](https://doi.org/10.1371/journal.pbio.1001091) reference genome, and (iv) [Allen Ancient DNA Resource (AADR)](https://doi.org/10.1101/2023.04.06.535797) *"1240k"* SNP variant callset).
2. Download and pileup our three input individuals `MT23` `MT26` and `MT7` using [`samtools mpileup`](http://www.htslib.org/doc/samtools-mpileup.html). Preprocessed alignment files for these samples are available on the European Nucleotide Archive, with accession number [PRJEB54894](https://www.ebi.ac.uk/ena/browser/view/PRJEB54894).
3. Index the 1000g-phase3 dataset into a set of FSA-encoded files, using the `fst` module of `grups-rs`.
4. Perform kinship analysis, by applying the `pedigree-sims` module on the input pileup file, and visualize results using `grups.plots`

---

#### Downloading the required datasets for `grups-rs`

1. Download the 1000g-phase3 dataset (chromosomes 1-22), along with its input panel definition file:
   ```bash
   mkdir -p data/1000g-phase3/
   URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
   wget -nd -r -P data/1000g-phase3/ ${URL}/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz*
   wget -nd -P data/1000g-phase3/ ${URL}/integrated_call_samples_v3.20130502.ALL.panel
   ```

2. Download the HapMapII recombination map (chromosomes 1-22).
   ```bash
   mkdir -p data/hapmap-phase2-b37/
   TARBALL="http://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz"
   wget -qO- ${TARBALL} | tar -xvzf- -C data/hapmap-phase2-b37/ --wildcards 'genetic_map_GRCh37_chr[0-9]*.txt'
   ```

3. Download a predefined pedigree definition file.
   ```bash
   mkdir -p data/pedigrees
   PED_URL="https://gist.githubusercontent.com/MaelLefeuvre/b6f3baff1964defea432ebd2c6dca6dc/raw/465ab9478724dba1575d20e45fb88f06055ac55f/simple-pedigree.txt"
   wget -P data/pedigrees ${PED_URL}
   ```

---

#### Generating an input pileup file

1. Download and decompress the `GRCh37` reference genome
   ```bash
   REFGEN_URL="http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
   mkdir -p data/reference/
   wget -P data/reference/ ${REFGEN_URL}
   bgzip -d data/reference/$(basename ${REFGEN_URL})
   ```

2. Download and format the [AADR](https://doi.org/10.1101/2023.04.06.535797) 1240K SNP callset
   ```bash
   mkdir -p data/targets
   AADR_SNP_URL="https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V52/V52.2/SHARE/public.dir/v52.2_1240K_public.snp"
   wget -O- ${AADR_SNP_URL} | awk 'BEGIN{OFS="\t"}{print $2, $4}' > data/targets/v52.2_1240K_public.bed
   ```

3. Download samples from ENA and pileup these using samtools
   1. Fetch the ENA URL of `MT23`, `MT26` and `MT7` ; Store these in `ena-bam-urls.tsv`
      ```bash
      ENA_TSV="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB54894&result=read_run&fields=submitted_ftp&format=tsv&download=true&limit=0"
      wget -qO- ${ENA_TSV} | grep -oP "ftp.*/MT[0-9]{1,2}.fastq.*.bam(?=;)" | sed 's/^/ftp:\/\//' > ena-bam-urls.tsv
      ```
   2. run `samtools mpileup` 
      ```bash
      samtools mpileup \
      -RBqQ25 \
      -l data/targets/v52.2_1240K_public.bed \
      -f data/reference/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
      -b ena-bam-urls.tsv \
      > mentesh-tepe.garino-vignon-commsbio-2023.RBqQ25.1240k.GRCh37.mpileup
      ```

---

#### Encoding the 1000g-phase3 dataset into the `.fst[.frq]` format

```bash
mkdir -p data/fst/
grups-rs fst --threads 22 --vcf-dir data/1000g-phase3/ --output-dir data/fst/EUR --pop-subset EUR --verbose
```

---

#### Running pedigree simulations with the `pedigree sims` modules

1. Run `pedigree-sims`
   ```bash
   grups-rs pedigree-sims \
   --pileup mentesh-tepe.garino-vignon-commsbio-2023.RBqQ25.1240k.GRCh37.mpileup \
   --data-dir data/fst/EUR \
   --recomb-dir data/hapmap-phase2-b37 \
   --output-dir grups-rs-output-mentesh-tepe \
   --pedigree data/pedigrees/simple-pedigree.txt \
   --samples 0-2 \
   --sample-names MT7 MT23 MT26 \
   --mode fst-mmap \
   --reps 1000 \
   --min-qual 25 \
   --pedigree-pop EUR \
   --contam-pop EUR \
   --verbose
   ```

2. Visualize results with `grups.plots`  
   This can be done in a single command, directly from bash...
   ```bash
   R --slave -e 'grups.plots::app(data_dir="grups-rs-output-mentesh-tepe")'
   ```
   ...or from an interactive R session.
   ```R
   library(grups.plots)
   grups.plots::app(data_dir="grups-rs-output-mentesh-tepe")
   ```

---

## Defining custom pedigrees

Defining pedigrees within GRUPS-rs is performed through simple definition files. See the main example template pedigree definition file [here](resources/pedigrees/example_pedigree.ped). Other examples may be found in the [resources/pedigrees](/resources/pedigrees) subdirectory of this repository. GRUPS-rs currently supports two alternative formats, which are described below.

### Standard format (`GRUPS-rs`)

The standard pedigree definition format of `grups-rs` can be subdivided into two sections:

1. The first section takes charge of defining the individuals found within the family tree, and its topology. Here, this particular section extensively, mirrors the commonly found [`.ped`](https://csg.sph.umich.edu/abecasis/Pedstats/tour/input.html) file format of the `PEDSTATS/QTDT` software, or PLINK's [`.fam`](https://www.cog-genomics.org/plink/1.9/formats#fam) files. 
   - Here, merely three columns are required from these previously mentionned file formats. Here, the `iid` column specifies the within-family id of the individual, while `fid` and `mid` both specify the ids of the individuals parents (See below):
     ```python
     # First section: define the pedigree's topology
     iid    fid   mid
     Ind1   0     0       # Ind1 and Ind2 are defined as founder individuals
     Ind2   0     0
     Ind3   Ind1  Ind2    # Ind3 is defined as the offspring of Ind1 and Ind2
     ```

2. The second section of the file takes charge of specifying which pairwise comparisons should `GRUPS-rs` specifically investigate within the template family tree. Here, every line beginning with the `COMPARE` keyword is considered as a "comparison definition" entry by `GRUPS-rs`, and is expected to adhere to the following scheme:
   ```
   COMPARE  <label> <iid-1> <iid-2>
   ```
   Where,
   - `<label>` is the user-defined name for the given comparison (e.g. 'first-degree', 'cousins', 'unrelated', etc.)
   - `<iid-1>` is the individual id of the first sample being compared
   - `<iid-2>` is the individual id of the second sample being compared

   Example:
   ```
   COMPARE Unrelated Ind1 Ind2
   COMPARE First     Ind1 Ind3
   COMPARE Self      Ind3 Ind3
   ```

Note that, while requiring only three columns, `GRUPS-rs` is able to directly parse the previously mentionned `.fam` and `.ped` formats, provided that users manually annotate the required second section at the bottom of these files. Hence, an quick and intuitive way to design custom pedigree files for GRUPS-rs is to:
1. *Visually* generate the first section of the file, using the interactive [`QuickPed`](https://magnusdv.github.io/pedsuite/articles/web_only/quickped.html) online software [(Vigeland M.D. 2022)](https://doi.org/10.1186/s12859-022-04759-y).
2. Export and save the output of QuickPed as a `.ped` file
3. Manually append the desired 'COMPARE' entries at the bottom of this file.
### Legacy format (`GRUPS`)

On top of the current standard format, `GRUPS-rs` remains entirely backwards compatible with the previous file format of `GRUPS``.

In essence, the legacy pedigree definition file of GRUPS is defined and parsed in three distinct steps, each one tied to a keyword within the definition file:

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

3. `COMPARISONS` Define which pairwise comparisons should GRUPS-rs investigate to compute genetic distances.
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

Note that comment lines are allowed within the definition file: any line starting with a `#` character is ignored during parsing.

### A more complete example

Below is a more complete definition file, where the intent would be to investigate the following six kinship ties:

| label           | description                                                                     |
| --------------- | ------------------------------------------------------------------------------- | 
| `inbred-self`   | compare an inbred individual to itself                                          | 
| `self`          | compare an outbred individual to itself (r=1)                                   | 
| `first`         | compare two first-degree relatives (r=0.5)                                      |
| `second`        | compare two outbred second-degree relatives (r=0.25)                            |
| `inbred-second` | compare two three-quarter siblings relatives due to inbreeding (r=0.25+0.125).  | 
| `unrelated`     | compare two unrelated individuals (r=0)                                         | 

Thus, an appropriate family tree topology could be as follows:  
  
<p align="center">
   <img src="https://github.com/MaelLefeuvre/grups-rs/assets/70585821/52616848-44e0-4005-a42f-03e2e7562438" />

</p>

Where founder and simulated individuals are colored in teal and lavander, respectively. Green arrows represents the comparisons that GRUPS-rs is requested to perform.  

Here, a template family tree such as this one can be defined as the following:

**Standard format**
```python
# standard format
Ind1  0     0
Ind2  0     0
Ind3  Ind1  Ind2                    # Ind3 and Ind4 are defined as the childreb of Ind1 and Ind2
Ind4  Ind1  Ind2
Ind5  Ind3  Ind4                    # Ind5 is defined as an inbred individual, since Ind3 and Ind4 are siblings.
Ind6  0     0
Ind7  Ind4  Ind6

COMPARE inbred-self    Ind5  Ind5
COMPARE self           Ind3  Ind3
COMPARE first          Ind2  Ind4
COMPARE inbred-second  Ind5  Ind7
COMPARE second         Ind2  Ind7
COMPARE Unrelated      Ind1  Ind2
```

<b><ins>Legacy format</b></ins>:

Alternatively, one could also define this family tree using the legacy format of `GRUPS` in such a manner: <b><inb><details><summary>Legacy format (click to unroll)</summary></ins></b>

```python
# legacy format
INDIVIDUALS
Ind1
Ind2
Ind3
Ind4
Ind5
Ind6
Ind7

RELATIONSHIPS
Ind3=repro(Ind1,Ind2)             # Ind3 is defined as the child of Ind1 and Ind2.
Ind4=repro(Ind1,Ind2)             # Ind4, as the child of Ind1 and Ind2.
Ind5=repro(Ind3,Ind4)             # Ind5, as the child of Ind3 and Ind4. (Ind5 is thus an inbred individual, since Ind3 and Ind4 are siblings)
Ind7=repro(Ind4,Ind6)             # Ind7, as the child of Ind4 and Ind6.

COMPARISONS
inbred-self=compare(Ind5,Ind5)    # Compare inbred individual Ind5 to itself.  label this relationship as 'inbred-self'
self=compare(Ind3,Ind3)           # Compare outbred individual Ind3 to itself. label this relationship as 'self'
first=compare(Ind2,Ind4)          # Compare Ind2 and Ind4. label this relationship as 'first'
inbred-second=compare(Ind5,Ind7)  # Compare Ind5 and Ind7. label this relationship as 'inbred-second'
second=compare(Ind2,Ind7)         # Compare Ind2 and Ind7. label this relationship as 'second'
unrelated=compare(Ind1,Ind2)      # Compare Ind1 and Ind2. label this relationship as 'unrelated'
```

Note that founder individuals (i.e. `Ind1`, `Ind2` and `Ind6`) are only defined in the `INDIVIDUALS` section and do *not* require to be defined in the `RELATIONSHIPS` section.  

</details>
<br>

Other example pedigree definition files may be found in the [resources/pedigrees](/resources/pedigrees) subdirectory of this repository.

### Specifying A genetic sex for simulated individuals 

Starting at version `v0.4.0`, Pedigree definition files may optionally contain an additional columns specifying information regarding the chromosomal sex of the individuals. This feature is only relevant when genetic relatedness analysis on the X-chromosoma, using [`--x-chromosome-mode`](#-x--x-chromosome-mode) (see section: [X-chromosomal analysis](#running-x-chromosomal-kinship-analysis-with-grups-rs) )or when using `--sex-specific-mode` 

A usable file is available for convenience here: [resources/pedigrees/chrX-pedigree.ped](/resources/pedigrees/chrX-pedigree.ped). 

An example of this optional column can be seen here :
<details><summary>Example X-chromosome aware pedigree (click to unroll) </summary>

  ```text
  # First section: define the pedigree's topology
  iid         fid     mid     sex
  father      0       0       1
  mother      0       0       2
  son1        father  mother  1
  son2        father  mother  1
  daughter1   father  mother  2
  daughter2   father  mother  2
  ```

</details>
<br>

---

## Applying genetic relatedness analysis on the X-chromosome with GRUPS-rs

Starting at version `v0.4.0`, GRUPS-rs now carries the ability to perform genetic relatedness analysis and pedigree simulations on the X-chromosome. Using `GRUPS-rs` in this mode will however require the use of specific pedigree definition files and panels

1. Your pedigree definition file `must` specify the sex of every simulated individual for this type of analysis to make sense. See the dedicated section here: [Specifying A genetic sex for simulated individuals](#specifying-a-genetic-sex-for-simulated-individuals). Example pedigree files may be consulted here in the [resources/pedigrees](/resources/pedigrees/) directory.

2. We *strongly* recommend that the provided [Input panel definition file](#2-input-panel-definition-file) contain an additional column specifying the sex of every reference individual, as in:
    ```txt
    HG00171 FIN     EUR     female
    HG00181 FIN     EUR     male
    HG01882 ACB     AFR     male
    ```
    ### The `fst` module: Encoding VCF files as a set of deterministic Finite State Acceptors

3. The set of [`.fst`](#the-fst-module-encoding-vcf-files-as-a-set-of-deterministic-finite-state-acceptors), or `vcf` files, provided trough the [`-F`|`--data-dir`](#-f--data-dir) parameter should of course contain genotype targets for the X-chromosome.

Once these conditions are met, one may use GRUPS-rs with the [`-X`|`--x-chromosome-mode`](#-x--x-chromosome-mode) argument, to specifically estimate genetic relatedness on X-chromosome:
```bash
grups-rs pedigree-sims --x-chromosome-mode --mode fst                         \
--data-dir tests/test-data/fst/binary-2FIN-1ACB-virtual-chrX/                 \
--recomb-dir tests/test-data/recombination-maps/chrX/                         \
--pileup tests/test-data/pileup/ash128-ash133.virtual.min-depth-5.chrX.pileup \
--pedigree resources/pedigrees/chrX-pedigree.ped                              \
--contam-pop AFR                                                              \
--reps 1000                                                                   \
--versbose
```

## Output Files

### `.pwd` file

`.pwd` output files contain summary statistics and results regarding all pairwise raw $\widehat{PWD}^{obs}$. This file is generated by the `pwd-from-stdin` module, are tab-separated and headed. 

| Column          | Type    | Description                                                                                                    |
| --------------- | ------- | -------------------------------------------------------------------------------------------------------------- |
| `Pair_name`     | string  | Descriptive label for a given pairwise comparison. Labels are in the form `<IND i>-<IND j>`                    |
| `Raw.Overlap`   | integer | Number of overlapping SNPs, which met the provided treshold of [`--min-depth`](#x--min-depth) on both samples  |
| `Raw.Sum.PWD`   | float   | Sum of the long-term average pairwise mismatch rates for all overlapping positions                             |
| `Raw.Avg.PWD`   | float   | Average Pairwise Mismatch Rate, i.e.: raw $\widehat{PWD}^{obs}$, or `Raw.Sum.PWD / Raw.Overlap`                |
| `Raw.CI.95`     | float   | Raw 95% Confidence interval for `Raw.Avg.PWD`                                                                  |
| `Raw.Avg.Phred` | float   | Average Phred score for all overlapping positions (Scale: PHRED-33)                                            |

### `.result` file

`.result` files contain summary statistics and results regarding pedigree simulations results. This most notably contains information regarding the most likely estimated relationship, given pedigree simulations results, as well as all pairwise corrected $\widehat{PWD}^{obs}. This file emanates from the `pedigree-sims` module, are tab-separated and headed.

| Column            | Type    | Description                                                                                                                                                 |
| ----------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Pair_name`       | string  | Descriptive label for a given pairwise comparison. Labels are in the form `<IND i>-<IND j>`                                                                 |
| `Most_Likely_rel` | string  | Estimated Most Likely Relationship, given pedigree simulations. Note that this column says nothing about significance                                       |
| `Corr.Overlap`    | integer | Corrected Number of overlapping SNPs, after filtering out positions not found within [`--data-dir`](#f--data-dir) or below the provided [`--maf`](#m--maf) threshold                   |
| `Corr.Sum.PWD`    | float   | Corrected um of long-term average pairwise mismatch rates, after filtering positions not found within [`--data-dir`](#f--data-dir), or below the provided [`--maf`](#m--maf) threshold |
| `Corr.Avg.PWD`    | float   | Corrected Average Pairwise Mismatch Rate, i.e.: $\widehat{PWD}^{obs}_{i,j}$, or `Corr.Sum.PWD / Corr.Overlap`                                               |
| `Corr.CI.95`      | float   | Corrected 95% confidence interval for `Corr.Avg.PWD`                                                                                                        |
| `Corr.Avg.Phred`  | float   | Corrected average  Phred score, after filtering positions not found within [`--data-dir`](#f--data-dir), or below the provided [`--maf`](#m--maf) threshold`                           |
| `Sim.Avg.PWD`     | float   | Average $\widehat{PWD}^{sim}$ for the most likely relationship (`Most_Likely_rel`) distribution.                                                            |
| `Min.Z_Score`     | float   | Z-score between `Corr.Avg.PWD` and the distribution of the most_likely relationship (`Most_Likely_rel`)                                                     |

### `.sims` files

`.sims` files contain pairwise-specific raw simulation results. In other terms, this file contains information, and most notably the $\widehat{PWD}^{sim}$ of each investigated kinship tie, within every pedigree replicate. These files are generated from the `pedigree-sims` module, are tab-separated and unheaded, and are located within the `simulations` subdirectory.

| Field number | Type    | Description                                                                                                                                             |
| ------------ | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 0            | integer | Pedigree simulation replicate index                                                                                                                     |
| 1            | string  | User-defined pedigree comparison label                                                                                                                  |
| 2            | string  | User-defined label of the first pedigree individual being compared                                                                                      |
| 3            | string  | User-defined label of the second pedigree individual being compared                                                                                     |
| 4            | string  | label of the randomly selected reference for the first individual being compared (Founder individuals only - `None`, when the individual is simulated)  |
| 5            | string  | label of the randomly selected refereence for the second individual being compared (Founder individuals only - `None` when the individual is simulated) |
| 6            | integer | Sum of simulated pairwise differences for this comparison                                                                                               |
| 7            | integer | Sum of simulated overlap for this comparison                                                                                                            |
| 8            | float   | Average simulated pairwise mismatch rate, or $\widehat{PWD}^{sim} for this comparison                                                                   |

### `.blk` files

`.blk` files contain pairwise-specific raw observed pairwise mismatch rates within non-overlapping windows. These files are generated from the `pwd-from-stdin` module, are tab-separated and unheaded, and located within the `blocks` subdirectory.

| Field number | Type    | Description                                                                                                           |
| ------------ | ------- | --------------------------------------------------------------------------------------------------------------------- |
| 0            | integer | Chromosome number                                                                                                     |
| 1            | integer | Block Start (bp)                                                                                                      |
| 2            | integer | Block End (bp)                                                                                                        |
| 3            | integer | Block-specific number of observed overlapping SNPs, which met the provided treshold of [`--min-depth`](#x--min-depth) |
| 4            | float   | Block-specific Sum of long-term average pairwise mismatch rates for all overlapping positions                         |

### `.probs` file

The `.probs` file is an optional output file of the `pedigree-sims` module, which is only emitted when [`--assign-method`](#--assign-method) is set to `svm`. This file will contain per-class SVM probabilities for every investigated relationship during simulations, and for every pairwise comparison.

### `.yaml` file

`.yaml` contain the serialized output of the `grups-rs` command line interface. These files may be generated by any submodule of `grups-rs`, strictly follow the [YAML format specification v1.2](https://yaml.org/spec/1.2.2/). The objective of these yaml files is mainly to ensure reproducibility, by keeping a record of all parameters used during a `grups-rs` command. These yaml files may also be used to re-run a `grups-rs` command using the exact same configuration, using the `from-yaml` submodule (See [the `from-yaml` module](#the-from-yaml-module-re-running-grups-rs-using-yaml-configuration-files) section, for a more detailled explanation as to how to use this submodule). 

## Citing GRUPS-rs

If you plan to use GRUPS-rs, please cite [(Lefeuvre M. et al 2024)](https://doi.org/10.47248/hpgg2404010001), as well as the the original publication of GRUPS [(Martin D. et al 2017)](https://doi.org/10.1111/mec.14188):

> Lefeuvre M, Martin M, Jay F, Marsolier M, Bon C. GRUPS-rs, a high-performance ancient DNA genetic relatedness estimation software relying on pedigree simulations. Hum Popul Genet Genom 2024; 4(1):0001. https://doi.org/10.47248/hpgg2404010001

> Martin MD, Jay F, Castellano S, Slatkin M. Determination of genetic relatedness from low-coverage human genome sequences using pedigree simulations. Mol Ecol. 2017;26(16):4145-4157. doi:10.1111/mec.14188 

Detailled citation instructions can be found by running the `grups-rs cite` module.

---

## Parameter List 
### Common flags
###### `-h`|`--help`
Print help information for a specific command.

###### `-v`|`--verbose`
Set the verbosity level of this program. With multiple levels:
- `-v`: Sets the INFO log level.  
- `-vv`: Sets the DEBUG log level.  
- `-vvv`: Sets the TRACE log level.  

###### `-V`|`--version`
Print version information.

Note that the program will still output warnings by default. Use [`--quiet`](#q--quiet) argument to
disable them

###### `-q`|`--quiet`
Disable warnings.

By default, warnings are emited and redirected to the console, even without verbose
mode on. Use this argument to disable this behavior and only display errors.


---

### `pwd-from-stdin` arguments and flags
#### Required arguments
###### `--pileup`
Specify an input pileup file containing samples to investigate.  

Note that in the absence of a `--pileup` argument, the program may accept a data stream
from the standard input. i.e:  
```bash
samtools mpileup -B -q30 -Q30 ./samples/* | pwd_from_stdin [args]
```

#### Optional arguments
###### `-t`|`--targets`
Provide with a list of SNP coordinates to target within the pileup.

Pileup positions which are not found within the provided `--targets` file will be exluded
from comparison.  

Accepted file formats:  
- '.snp' : EIGENSTRAT format. See the [convertf documentation](https://github.com/argriffing/eigensoft/blob/master/CONVERTF/README) for additional information regarding this type of file.
- '.vcf' : Variant Call Format. See the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf) for additional information regarding this type of file.
- '.tsv' : a tab-separated file with columns `<CHROMOSOME> <POSITION> <REFERENCE> <ALTERNATE>`
- '.csv' : a comma-separated file with columns `<CHROMOSOME> <POSITION> <REFERENCE> <ALTERNATE>`
- '.txt' : a space-separated file with columns `<CHROMOSOME> <POSITION> <REFERENCE> <ALTERNATE>`

###### `-b`|`--blocksize`
Change the window size of each jackknife block (found within the output [`.blk` files](#blk-files)).

Note that lower block values may drastically increase the memory footprint of `grups-rs`.

###### `-c`|`--chr`
Restrict comparison to a given set of chromosomes.

This argument may accept slices, such as `--chr 9-11` and/or discrete integers such as `--chr 1 4 13`.  

###### `-g`|`--genome`
Fasta indexed reference genome.

By default, `grups-rs` will use GRCh37 as a reference genome. Use this argument if you wish to specify an alternative reference. Note that a `.fasta.fai` genome index file must be present at the same directory.

**Example:** specifying `--chr 9-11 13 19-22` ...will be parsed as: `[9, 10, 11, 13, 19, 20, 21, 22]`

###### `-M`|`--min-qual`

Minimal required Base Quality (BQ) score to select an observation as a valid candidate for comparison.

Nucleotides whose base quality is lower than the provided threshold are simply filtered-out prior to comparison. Value should be expressed in the Phred+33 scale.

###### `-n`|`--sample-names`
Provide with a list of sample names for printing.

By default, individuals will be referred to using their pileup index. e.g. "Ind0, Ind1, Ind2, etc."

The provided list must of course follow the sorted index order which was provided by [`--samples`](#s--samples). Note that the length of [`--samples`](#s--samples) and `--sample-names` do not need to match. When such is the case, default values are used for the missing names.

###### `-o`|`--output-dir`
Specify the output directory where results will be written.

Note that grups-rs will create the specified leaf directory if it is not present. However, `grups-rs` does not allow itself from creating parent directories.

###### `-s`|`--samples`
0-based column index of the individuals that should be compared within the file specified by [`--pileup`](#pileup)

Argument may accept slices (inclusive), such as `--samples 0-3` and/or discrete integer values such as `--samples 7 8`.

**Example:** `--samples 0-3 7 8` will be parsed as: `[0, 1, 2, 3, 7, 8]`

###### `-x`|`--min-depth`
Provide with the minimal sequencing depth required to perform comparison.

Note that the length of the provided vector does not need to be the same as the one provided for individuals. When such is the case, values of --min-depth will "wrap" around, and start again from the beginning for the next individuals.

**Example:** Specifying `--samples 0-4 --min-depth 2 3 5` would result in:
| Sample | Depth |
| ------ | ----- |
| 0      | 2     | 
| 1      | 3     | 
| 2      | 5     | 
| 3      | 2     | 
| 4      | 3     | 

###### `-X`|`--X-chromosome-mode`
Run GRUPS-rs in X-chromosome comparison mode.

This mode will trigger specific rules of allele transmission and recombination during pedigree simulations, reflective of X-chromosomal inheritance rules.
Note that this mode requires the use of specific pedigree definition files and panels, containing information regarding the chromosomal sex of individuals. See the dedicated
section regarding the use of this mode here: [Applying genetic relatedness analysis on the X-chromosome with GRUPS-rs](#applying-genetic-relatedness-analysis-on-the-x-chromosome-with-grups-rs)

#### Optional flags 
###### `-f`|`--filter-sites`
Do not perform comparison, but rather print out the pileup lines where a comparison
could be made.

This is a legacy argument from the initial implementation of GRUPS. Note that lines are printed as soon as there is a valid comparison between ANY pair of individuals. Thus, applying [`--filter-sites`](#f--filter-sites) may not prove very effective when used in conjunction with [`--self-comparison`](#s--self-comparison)

###### `-i`|`--consider-dels`
Consider deletions when performing comparison.

By default, deletions ('*' character) are not counted as selectible positions for comparison. Using this flag will instead mark them as valid matching positions.

###### `--exclude-transitions`
Exclude transitions from the input targets file.

Note that this argument requires the use of [`--targets`](#t--targets) to provide the program with a list of coordinates. The given coordinate file must furthermore explicitely provide with the REF/ALT alleles. See the [`--targets`](#t--targets) argument for additiotnal information regarding valid file formats.

###### `-J|--no-print-blocks`
Do not output jackknife blocks for each pair of individuals (`.blk` files).

By default, `grups-rs` will keep track of the pairwise mismatch rate within windows of size [`--blocksize`](#b--blocksize). This information can prove useful to investigate PMR values in sliding windows (this can be visualized with the `grups.plots` companion shiny interface), or to compute jackknife estimates of variance.

Using --no-print-block will prevent `grups-rs` from computing average pairwise differences in non-overlapping blocks, and will prevent it from outputting any `.blk` files.

###### `-k`|`--known-variants`
Filter out tri-allelic sites when given a list of SNP-coordinates.

Note that this arguement requires the use of `--targets` to provide the program with a
list of coordinates. The given coordinate file must furthermore explicitely provide with
the REFERENCE and ALTERNATE alleles.

See the `--targets` argument for additional information regarding valid formats.

###### `-S`|`--self-comparison`
Enable self-comparison mode on individuals.

Note that a minimal depth of 2 is required when comparing individuals to themselves,
since at least two alleles are required to effectively compute pairwise differences at a
given site. This can be activated by specifying a library-wise value of `--min-depth`>=2
for all samples.

Note that in cases where the specified `--min-depth` was set to 1 for a given individual,
grups-rs will automatically *temporarily* rescale the given value of `--min-depth` to 2 when
performing self-comparisons.

###### `-w`|`--overwrite`
Overwrite existing output files.  

By default, grups-rs does not allow itself from overwriting existing results files. Use
this flag to force this behaviour.

---
### `pedigree-sims` module arguments and flags.
#### Required arguments
###### `-F`|`--data-dir`
Path to a directory containing a database of phased modern human genotypes, such as the 1000g-phase3 dataset.

This database may be in the form of a set of VCF files (default). In that case, grups-rs will look for, and load any file ending with the .vcf[.gz] file extension.

Alternatively, users may use a set of FSA-encoded files. To use FSA-encoded files, users must explicitly specify the use of this type of data input, with the --mode argument. When specified, grups-rs will search for any file
ending with the .fst[.frq] file extension. 

See the documentation of the [`fst`](#the-fst-module-encoding-vcf-files-as-a-set-of-deterministic-finite-state-acceptors) module, for a more detailled explanation regarding how to generate FSA-encoded databases.

###### `-G`|`--recomb-dir`
Path to a directory containing a set of chromosome-specific genetic recombination maps, such as the HapMap-phaseII dataset.
 
Note that GRUPS-rs will search for, and attempt to load any file ending with the '.txt' file extension within that directory. It is therefore highly recommended that this directory ONLY contains the required recombination map files, and nothing else.

The expected input is a headed, tab separated file with columns `<chr> <pos(bp)> <rate(cM/Mb)> <Map(cM)>`

<details><summary>Example recombination file (click to unroll)</summary>

    Chromosome  Position(bp)    Rate(cM/Mb)     Map(cM)
    chr22       16051347        8.096992        0.000000
    chr22       16052618        8.131520        0.010291
    ...         ...             ...             ...

</details>

###### `-T`|`--pedigree`
Path to input pedigree definition file. Examples of such definition files may be found within the `resources/pedigrees` subdirectory of this github repository. See section [Defining Custom pedigrees](#defining-custom-pedigrees) for a detailled explanation on how to write custom input pedigree definition files.

#### Optional arguments

###### `-d`|`--af-downsampling-rate`
Allele frequency downsampling rate, i.e. the proportion of SNPs to keep at true frequency (in percentage).

`grups-rs` will use the specified `--af-downsampling-rate`` as a probability to randomly simulate allele fixation during pedigree simulations.

Note that the defined rates should be specified as percentages.

###### `-D`|`--snp-downsampling-rate`
Proportion of filtered SNP positions to include in the analysis (in percentage).

`grups-rs` will use the specified `--snp-downsampling-rate` as a probability to randomly ignore positions during pedigree simulations. 

Note that the defined rates should be specified as percentages.

###### `-Q`|`--contam-rate`
Contamination rates (or rate ranges) for each pileup individual. 
 
The provided argument(s) may accept hard set values, such as `--contam-rate 3.0`, or ranges, such as `--contam-rate 0.0-5.0`.
When a contamination rate is provided in the form of a range, pedigree-specific values are picked from a uniform distribution, within that defined range.

Note that contamination rates are specified as percentage, and are tied to the individuals contained within the input pileup. Specifying a constant rate for all individuals is also possible, by defining a unique value or range.

**Example: `grups-rs pedigree-sims [...] --samples 0-5 --contam-rate 2.0` implies that all five pileup indidivuals will be assigned a set contamination rate of 2%, during pedigree simulations. 

In general, keep in mind that contamination rate values are recycled, if the number of specified values is lower than the number of examined pileup individuals.

###### `-U`|`--seq-error-rate`
Sequencing error rates (or rate ranges) for each pileup individual.
 
The provided argument(s) may accept hard set values, such as `--seq-error-rate 1.0`, or ranges, such as `--seq-error-rate 1.0-3.0`. When a sequencing error rate is provided in the form of a range, pedigree-specific values are picked from a uniform distribution, within that defined range.

Note that sequencing error rates are specified as percentage, and are tied to the individuals contained within the input pileup. Specifying a constant rate for all individuals is also possible, by defining a unique value or range.

Example: `grups-rs pedigree-sims [...] --samples 0-7 --seq-error-rate 1.0` implies that all seven pileup indidivuals will be assigned a set sequecing error rate of 1%, during pedigree simulations.
In general, keep in mind that sequencing error rate values are recycled if the number of specified values is lower than the number of examined pileup individuals.

###### `-I`|`--mode`
Define the expected data input type for pedigree simulations.

This argument is closely tied to the [`--data-dir`](#f--data-dir) argument, and will define which type of files should `grups-rs` look for, as well as how to load them into memory. 

(**tl;dr:** `--mode fst-mmap` is recommended for most applications. Use `--mode fst` when runtime performance is critical, but memory usage is not an issue.)

- When using `--mode fst-mmap` `grups-rs` will search for files ending with the `.fst[.frq]` file extention, and will make use of [memory-mapped files](https://en.wikipedia.org/wiki/Memory-mapped_file) to query this reference dataset. This mode is highly recommended if the directory targeted with [`--data-dir`](#f--data-dir) is located within an [SSD drive](https://en.wikipedia.org/wiki/Solid-state_drive). However, applying the *fst-mmap* when working with data found within an [HDD drive](https://en.wikipedia.org/wiki/Hard_disk_drive) may have a slight negative impact on performance.

- When using `--mode fst`, `grups-rs` will search for files ending with the `.fst[.frq]` file extention. These files are then sequentially loaded into RAM and queried. The *fst* mode may be faster if you have a lot of comparisons to apply, but has has a higher memory footprint compared to the *fst-mmap* mode. Thus, it is only recommended if your input `.fst[.frq]` files are located on an HDD drive.

- When using `--mode vcf` (default), `grups-rs` will search for files ending with the `.vcf[.gz]` extention, and will directly query these files. This mode is generally slower, but may yield higher performance if you have highly covered data, and your vcf files were previously pre-filtered to only contain bi-allelic SNPs.

###### `-R`|`--reps`
Number of pedigree simulation replicates to perform for each pairwise comparisons.

The default provided value of 100 should be considered a bare-minimum, for quick screening. Values in the range of 500 to 1000 replicates is recommended.

###### `-m`|`--maf`
Minimum required (super-)population allele frequency to include a given SNP during simulations.

Note that `grups-rs` will re-compute the observed average PWD after simulations, to excluding and correct for any position previously excluded through this threshold.

###### `-P`|`--pedigree-pop`
Source population from which founder individuals are selected during pedigree simulations.

Note that when using [`--mode vcf`](#i--mode)`--mode vcf`, grups-rs may only use populations for which a `<POP>_AF` annotation is present in each INFO field of the VCF files used. To generate finer grained population allele frequencies, we recommend either the use of the [`bcftools +fill-tags`](https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html) plugin, or the use of the ['--compute-pop-afs'](#f--compute-pop-afs) argument when generating FSA-encoded dataset with `grups-rs`' `fst` module.

###### `-C`|`--contam-pop`

Population from which contaminating individuals are selected during pedigree simulations.  

Note that when using [`--mode vcf`](#i--mode), grups-rs may only use populations for which a <POP>_AF annotation is present in each INFO field of the VCF files used. To generate finer grained population allele frequencies, we recommend
either the use of the [`bcftools +fill-tags`](https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html) plugin, or the use of the ['--compute-pop-afs'](#f--compute-pop-afs)'--compute-pop-afs' argument when generating FSA-encoded dataset with the `grups-rs`' `fst` module.

###### `-N`|`--contam-num-ind`
Number of random individual genomes with which to contaminate pedigree simulations.

Note that the number of contaminating individuals, are tied to the samples contained within the input pileup, and that the specified values will be recycled if their length is lower than the number of examined pileup samples.

###### `-p`|`--panel`
Path to an input panel definition file.

By default, `grups-rs` will automatically search for a file ending with the `.panel` extension within the directory targeted by [`--data-dir`](#f--data-dir). Use the `--panel` argument to override this behaviour, and specify a definition file located somewhere else.

###### `--decompression-threads`
Number of additional parallel decompression threads.
 
Can increase performance when working with BGZF compressed `.vcf.gz` files. Note that this parameter has no effect when working with uncompressed `.vcf` or `.fst[.frq]` files.

###### `--sex-specific-mode`
Run GRUPS-rs in sex-specific mode

By default, grups-rs will randomly pick reference samples as founder individuals, without consideration of their original chromosomal sex.
With the use of `--sex-specific-mode`, pedigree samples are instead randomly assigned a chromosomal sex. Reference samples are then selected
in accordance with the sex of the considered founder individual.

###### `--seed`
Provide the random number generator with a set seed.

###### `--assign-method`
Select the method for most likely relationship assignment

**tl;dr**: use the default `--assign-method svm` ; only use `--assign-method zscore` if you wish to replicate the classification method of the first version of GRUPS.  

- **svm**: Assign a most likely relationship using Ordinally Partitionned Support Vector Machines. Binary SVMs are instantiated sequentially, from the lowest relatedness order to the highest, in terms of average PWD. Each SVM is fitted against the hypothesis that the observed PWD belongs to a higher degree than the given distribution. These SVMs are then used to generate per-class probabilities of belonging to a given class of relationship. The most likely relationship is assigned by selected the relationship class with the highest per-class probability.

- **zscore**: Perform minimum zscore assignation, i.e. the distribution with the lowest z-score from the observed PWD is selected as the most likely candidate. This approach is computationally inexpensive, but may provide with spurious results, when the different distributions carry drastically different standard deviations.

---

### FST Index
#### Required parameters
###### `-d`|`--vcf-dir`
Path leading to the directory containing the set of input reference VCF files that you wish to encode.

###### `-o`|`--output-dir`
Output directory where the generated FSA-encoded files (ending with the `.frq[.fst]` extension) should be written to.


#### Optional parameters
###### `-p`|`--panel`
Path to an input panel definition file.

By default, `grups-rs` will automatically search for a file ending with the `.panel` extension within the directory targeted by [`--data-dir`](#f--data-dir). Use the `--panel` argument to override this behaviour, and specify a definition file located somewhere else.

###### `-P`|`--pop-subset`
Population Subset

Subset the index by a given number of (super)-population. (e.g. EUR, AFR). Note that at least one pedigree population and one contaminating population are required to obtain valid index files. (Although, the source and contaminating population may be the same.)

###### `-@`|`--threads`
Number of parallel CPU processes when performing FST-Indexation

Note that parallelization is dispatched according to the number of separate `.vcf[.gz] files`. Thus, there is no point in invoking more threads than there are reference VCF files to encode.

###### `--decompression-threads`
Number of additional parallel decompression threads when decompressing BGZF-compressed VCF files. (This parameter has no effect when working with uncompressed `.vcf` files.)

Also note that decompression threads have a multiplicative effect when combined with '--threads'. Thus, setting `--decompression-threads 2` and `--threads 22`, will in fact consume up to 44 worker threads.

#### Optional flags
###### `-F`|`--compute-pop-afs`
Recalculate population allele frequencies for each population and super-population tag that can be found within the provided input
panel definition file (specified with [`--panel`](#p--panel-1)). 

For some publicly available datasets, such as the 1000g-phase3 adding this flag can allow for the use of smaller populations as reference or to use customly defined populations.

When unspecified, the program will instead look for `<POP>_AF` tags within the VCF's INFO field. These tags can be generated using the [`bcftools +fill-tags`](https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html) plugin.

## Contributing

### Obtaining code coverage metrics for GRUPS-rs

1. Install cargo tarpaulin  
    ```
    cargo install cargo-tarpaulin
    ```

 2. Run coverage tests  
    ```
    cargo tarpaulin --workspace --command test --out Lcov
    ```

### Building the GRUPS API docs.

The complete documentation for the GRUPS-rs public API can be built using `cargo`:
```Bash
cargo doc --workspace --open --no-deps --document-private-items
```
Using this command, the documentation's `html` entry-point will be located under `target/doc/pedigree_sims/index.html` and should automatically open in your default browser.

