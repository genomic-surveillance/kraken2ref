# Changelog

[2.1.0] 2024-10-03
---
### Changed

- [refactor] sort_reads now split into two separate processes
- [improvement] Unit tests updated and improved
- [improvement] Explicit handling of kraken report being an empty file
- [improvement] Updated README

### Added

- [feature] Default reference-selection method now "max"
- [feature] Default behaviour to select only a single reference from a given valid taxonomy tree
- [feature] `sort_reads` now produces a taxid-to-read_id mapping as a JSON file
- [feature] `dump_fastq` now a separate script, reads in `sort_reads` output JSON to produce per-taxon fastqs
- [feature] Unit test for `dump_fastq`

---
[2.0.0] 2024-05-17
---
### Changed
- [improvement] Major refactor; polling and taxonomy tree handling are now encapsulated in classes
- [improvement] Refactored `sort_reads` module
- [usage] kraken2ref executable changed from `kraken2r` to `kraken2ref`
- [format] Moved from `src/` structure to flat structure
- [format] Moved most setup code from `setup.cfg` to `pyproject.toml`

### Added
- [improvement] Added test coverage for all modules/functionalities
- [feature] Added kmeans- and quantile-based polling methods; kmeans now the default
- [feature] Modification to the control-flow to allow for rescue of `parent_selected` results
- [feature] Added unique mode in `sort_reads` module
- [feature] Allow user to specify what suffix to use for `parse_report` output JSON
- [resource] Added python notebook in `tests/resources` to visualise/compare polling method behaviours

---
[1.1.0] 2024-04-02
---
Minor version built to accomodate changes to the flu taxonomy structure in the kraken2 DB

### Added
- [improvement] Added condensed mode to reduce data duplilcation in `sort_reads`
- [improvement] Added `generate_report` as command that is installed with kraken2ref
- [improvement] `generate_report` now outputs segment numbers for Flu as well

---
[1.0.0] 2024-03-05
---
First stable release.  
- kraken2ref used two modes called from subparsers:
    - `parse_report` analyses kraken2 taxonomic report and outputs JSON file
    - `sort_reads` uses kraken2 output (readID to taxID mapping) and `parse_report` output JSON file to sort reads by taxID