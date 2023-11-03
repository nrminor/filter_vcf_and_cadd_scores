# Filter VCF and CADD Scores together
[![Open Source Files](https://github.com/nrminor/filter_vcf_and_cadd_scores/actions/workflows/open-source-starter.yml/badge.svg)](https://github.com/nrminor/filter_vcf_and_cadd_scores/actions/workflows/open-source-starter.yml)

Filtering VCFs can be tedious. The Python script in this repo wraps `bcftools` and `vcftools` to make filtering to a desired subset of samples a little easier. It also brings large tables of [CADD Scores](https://cadd.gs.washington.edu/) with it, using [Polars LazyFrames](https://pola-rs.github.io/polars/py-polars/html/reference/lazyframe/index.html) to do so. Because it coordinates many moving parts, the script also uses elegant, Rust-like error handling thanks to the Python [`result`](https://pypi.org/project/result/) module.

To make this module's results more reproducible, the module is placed within a Python environment managed with [Poetry](https://python-poetry.org/). To set up a Poetry environment for this project of your own, simply install Poetry at the link about run the following steps:
1. Clone this repository with `git clone https://github.com/nrminor/filter_vcf_and_cadd_scores.git .`
2. Run `poetry install` to make sure all the required Python libraries are present in an isolated environment.
3. Run the module with:
```
poetry run filter_vcf_and_cadd_scores/main.py \
--vcf /path/to/variants.vcf \
--animal_file /path/to/animals.txt \
--cadd_table /path/to/cadd_scores.tsv \
--label "filtered_snps" \
--verbose true
```

See the docs for more details. Of note, this module is not yet backward compatible and was build with Python 3.11. Feel free to raise an issue if support for an earlier versions is needed. For now, users can run the module as-is with a Docker image built from the provided Dockerfile.
