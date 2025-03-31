# 0.4.0
## Features
- <b><ins>New</ins></b>: `--x-chromosome-mode` and `--sex-specific-mode` features, now allows the use of GRUPS-rs to perform genetic relatedness analysis on the X-chromosome

## Dependencies
  - Update `fastrand` to version `v2.0.0`
  - Remove unneeded `rand` dependency in multiple crate
  - bump `grups.plots` version to `v0.3.2`
    - new function `grups.plots::plot_all`, allowing non interactive plotting of grups-rs results
    - Allow the flexible importation of headed `.blk` and `.sims` files
    - Added Zenodo DOI Badge
## Quality-of-life
- `.blk` and `.sims` files now contain a descriptive header
- Improved compile time, through the harmonization of workspace dependencies and versioning

# 0.3.3
## Documentation
- Update citing information in `README.md` and `grups-rs cite` module.

## Minor changes
- Update Cargo.lock
- Remove unneeded `atty` dependency (Not this change bumps the MSRV to `cargo>=1.70.0`)

## Bugfixes
- Prevent the program from displaying signed zeroes (i.e.: `-0.0`) on `x32`bit architectures
- Fix grups.plots submodule URL format (set it back to http)

# 0.3.2
## Features
- Implement a simpler pedigree definition file format, that is based upon the [`.ped`](https://csg.sph.umich.edu/abecasis/QTDT/docs/pedigree.html) and/or [`.fam`](https://www.cog-genomics.org/plink/1.9/formats#fam) file formats.
- This effectively makes GRUPS-rs *somewhat* compatible with these usual file format, as users still need to target specific pairwise comparisons within the constructed tree. This is done through the 'COMPARE' keyword (see the documentation for more information on the current standard)
- Note that GRUPS-rs is still backwards compatible with its *'legacy'* format, which mirrors the initial implement of GRUPS. Here the program is able to automatically detect and parse the appropriate format.


# 0.3.1
## Features
- Add Jemalloc memory allocator

## Bugfixes
- Fix libsvm-rust's use afterfree, causing random SIGSEGV (signal 11) crashes.

# 0.3.0
## Bugfixes
- Fix confusing double negation in '--consider-dels' flag
