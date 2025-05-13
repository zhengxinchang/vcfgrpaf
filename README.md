# vcfgrpaf
vcfgrpaf is a fast, Rust-based implementation of grpaf.py from Truvari, designed for calculating allele frequency statistics per group from VCF files.

# Features
Reads genotype data from standard VCF files.

Computes per-group statistics for user-defined sample groups:

AF, MAF, MAC, AC, AN, N_HEMI, N_MISS, N_HOMREF, N_HET, N_HOMALT

Writes these statistics to new INFO fields in the output VCF.

Automatically removes redundant or outdated tags related to these stats.

# Installation
Ensure you have the Rust toolchain installed:

```bash

cargo build --release
The compiled binary will be located at target/release/vcfgrpaf.
```

# Usage

```bash

vcfgrpaf <VCF> --labels <LABELS> --output <OUTPUT>
Arguments
<VCF>: Input VCF file.

--labels <LABELS>: A tab-delimited file with two columns <sample> <group>, mapping samples to groups.

--output <OUTPUT>: Path to the output VCF file (use - to write to stdout).
```
Example

```bash

vcfgrpaf input.vcf.gz --labels groups.tsv --output grouped_output.vcf
```

Contents of groups.tsv:

```
sample1    groupA
sample2    groupA
sample3    groupB
...
```

# Output INFO Tags
Each group will result in the following tags in the VCF header:

```mathematica
##INFO=<ID=AF_groupA,Number=1,Type=Float,Description="AF on N groupA samples">
##INFO=<ID=MAC_groupB,Number=A,Type=Integer,Description="MAC on N groupB samples">
...
```
