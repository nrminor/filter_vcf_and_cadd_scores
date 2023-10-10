# Filter VCF and CADD Scores together
[![Open Source Files](https://github.com/nrminor/filter_vcf_and_cadd_scores/actions/workflows/open-source-starter.yml/badge.svg)](https://github.com/nrminor/filter_vcf_and_cadd_scores/actions/workflows/open-source-starter.yml)

Filtering VCFs can be tedious. The Python script in this repo wraps `bcftools` and `vcftools` to make filtering to a desired subset of samples a little easier. It also brings large tables of [CADD Scores](https://cadd.gs.washington.edu/) with it, using [Polars LazyFrames](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) to do so. Because it coordinates many moving parts, the script also uses elegant, Rust-like error handling thanks to the Python [`result`](https://pypi.org/project/result/) module.
