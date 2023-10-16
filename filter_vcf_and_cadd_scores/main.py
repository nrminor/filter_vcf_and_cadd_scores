#!/usr/bin/env python3

"""
DOWNSAMPLE A MULTI-SAMPLE VCF AND ASSOCIATED CADD SCORES
--------------------------------------------------------
"""

import os
import sys
import argparse
import subprocess
import glob
import gzip
from dataclasses import dataclass
from result import Result, Ok, Err
import polars as pl

@dataclass
class NewFile:
    """
    A data class for returning a new file path resulting from a successful
    subprocess run.
    """
    name: str
    path: str

def parse_command_line_args() -> Result[argparse.Namespace, str]:
    """
    Parses command line arguments.

    Returns:
        Result[argparse.Namespace, str]: Returns Ok(argparse.Namespace) if args could 
        be parsed, else returns Err(str).
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", "-v",
                        type=str,
                        required=True,
                        help="Multisample VCF to subset.")
    parser.add_argument("--cadd_table", "-c",
                        type=str,
                        required=True,
                        help="CADD Score spreadsheet in TSV format.")
    parser.add_argument("--animal_file", "-a",
                        type=str,
                        required=True,
                        help="Text file with one animal name per line.")
    parser.add_argument("--label", "-l",
                        type=str,
                        required=True,
                        help="A label for the data you are running on.")
    args = parser.parse_args()
    
    return Ok(args)

def unify_sample_names(animals: pl.DataFrame, vcf_path: str, label: str) -> Result[NewFile, str]:
    """
        This function double checks the VCF at the provided path to make
        sure the sample names in the VCF are the same as the sample names
        in the metadata. If they aren't, it looks for the metadata sample
        names as substrings of the VCF sample names. If all metadata 
        sample names are in fact substrings of the VCF sample names, it
        creates a tab-delimited text file to guide BCFTools, which it
        uses in a subprocess to rename VCF sample names.

    Args:
        animals: A Polars DataFrame of the sample IDs.
        vcf_path: a character string filepath leading to the input VCF.

    Returns:
        Result[NewFile, str]: A Result containing either an instance of the 
        NewFile dataclass, which contains the  file path to the new VCF output
        by VCFtools, or an informative error message.
    """

    if not vcf_path.endswith(".vcf.gz") and not vcf_path.endswith(".vcf"):
        return Err("VCF is not in a recognizable format. Please double \
                   check that your VCF file is either gzip-compressed and \
                   ends with '.vcf.gz', or is uncompressed and ends with \
                   '.vcf'.")

    # create a list of sample names from the VCF
    vcf_samples = []
    if vcf_path.endswith(".gz"):
        with gzip.open(vcf_path, 'rb') as vcf_file:
            for line in vcf_file:
                if line.startswith(b'#CHROM'):
                    # Split the line by tabs, and extract the sample IDs
                    # starting from the 9th column (0-indexed)
                    vcf_samples = line.decode('utf-8').strip().split("\t")[9:]
                    break
    if vcf_path.endswith(".vcf"):
        with open(vcf_path, 'r', encoding="utf-8") as vcf_file:
            for line in vcf_file:
                if line.startswith('#CHROM'):
                    # Split the line by tabs, and extract the sample IDs
                    # starting from the 9th column (0-indexed)
                    vcf_samples = line.strip().split('\t')[9:]
                    break
    if len(vcf_samples) == 0:
        return Err("Failed to define any sample names in the provided VCF.")

    # create a list of the sample names from the metadata
    meta_samples = animals.select(
        "Sample ID"
    ).drop_nulls().to_series().to_list()
    
    # make sure the animal IDs are the same types
    if not all(type(meta_sample) == type(vcf_sample) 
               for meta_sample, vcf_sample in zip(meta_samples, vcf_samples)):
        meta_samples = [str(sample) for sample in meta_samples]
        vcf_samples = [str(sample) for sample in vcf_samples]
        
    # match the two lists, permitting no samples to be unmatched
    matched_sample_ids = []
    for meta_sample in meta_samples:
        for vcf_sample in vcf_samples:
            if meta_sample in vcf_sample:
                matched_sample_ids.append((vcf_sample, meta_sample))
                break
    if len(matched_sample_ids) == 0:
        return Err("Failed to match sample IDs between VCF and provided animal names")
    
    # return the original VCF if the paths are the same
    if all(x == y for x, y in matched_sample_ids):
        original_vcf = NewFile(vcf_path, vcf_path)
        return Ok(original_vcf)

    # write a two column text file, where the first column is names
    # in the VCF (old names), and the second column is names in the
    # metadata (new names)
    with open('new_sample_names.txt', 'w', encoding="utf-8") as out_handle:
        for line in matched_sample_ids:
            print(f"{line[0]}\t{line[1]}", file=out_handle)

    # create the BCFTools command and run it in a subprocess
    new_vcf = f"{label}_renamed_samples.vcf.gz"
    with open(new_vcf, "w", encoding="utf-8") as vcf_handle:
        bcftools_process = subprocess.Popen(("bcftools", "reheader", "--samples",
                           "new_sample_names.txt", vcf_path), 
                           stdout=subprocess.PIPE)
        gzip_process = subprocess.Popen(["gzip", "-c"], stdin=bcftools_process.stdout,
                                        stdout=vcf_handle, stderr=subprocess.PIPE)

        # Allow vcftools to receive a SIGPIPE if gzip exits
        if bcftools_process.stdout:
            bcftools_process.stdout.close()

        # Wait for processes to complete and get their exit codes and stderr
        _, bcftools_stderr = bcftools_process.communicate()
        _, gzip_stderr = gzip_process.communicate()

    # Get the exit code and give an informative error message if it isn't 0
    bcftools_exit_code = bcftools_process.returncode
    if bcftools_exit_code != 0:
        return Err(f"VCF renaming encountered an error with exit code {bcftools_exit_code}:\n\
                       {bcftools_stderr}")
    gzip_exit_code = gzip_process.returncode
    if gzip_exit_code != 0:
        return Err(f"VCF compression encountered an error with exit code {gzip_exit_code}:\n\
                       {gzip_stderr}")

    # clear out temporary sample name file
    os.remove("new_sample_names.txt")

    new_file = NewFile(new_vcf, os.path.abspath(new_vcf))

    return Ok(new_file)

def filter_vcf_samples(vcf_path: str, animals: str, label: str) -> Result[NewFile, str]:
    """
        Run VCFTools to filter a VCF based on a text file of sample IDs to keep.

    Args:
        vcf_path (str): The path to the input VCF file.
        animals (str): The path to the text file containing sample IDs to keep.

    Returns:
        Result[NewFile, str]: A Result containing either the new file path to
        the filtered filtered VCF, or an informative error message.
    """

    # handle the possibility that the VCF is not named or formatted correctly
    if not vcf_path.endswith(".vcf.gz") and not vcf_path.endswith(".vcf"):
        return Err("VCF is not in a recognizable format. Please double \
                   check that your VCF file is either gzip-compressed and \
                   ends with '.vcf.gz', or is uncompressed and ends with \
                   '.vcf'.")

    # set the VCFtools input format flag based on the VCF file name is now
    # guaranteed to be present thanks to the above error-handling
    if vcf_path.endswith(".vcf.gz"):
        input_format = "--gzvcf"
    else:
        input_format = "--vcf"

    # construct the VCFTools subprocess command as a list
    command = [
        "vcftools",
        input_format, vcf_path,
        "--keep", animals,
        "--recode", 
        "--stdout"
    ]

    # run VCFtools
    new_vcf = f"{label}_filtered.vcf.gz"
    with open(new_vcf, "wb") as vcf_handle:
        vcftools_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        gzip_process = subprocess.Popen(["gzip", "-c"], stdin=vcftools_process.stdout,
                                        stdout=vcf_handle, stderr=subprocess.PIPE)

        # Allow vcftools to receive a SIGPIPE if gzip exits
        if vcftools_process.stdout:
            vcftools_process.stdout.close()

        # Wait for processes to complete and get their exit codes and stderr
        _, vcftools_stderr = vcftools_process.communicate()
        _, gzip_stderr = gzip_process.communicate()

    # Get the exit code and give an informative error message if it isn't 0
    vcftools_exit_code = vcftools_process.returncode
    if vcftools_exit_code != 0:
        return Err(f"VCF filtering encountered an error with exit code {vcftools_exit_code}:\n\
                       {vcftools_stderr}")
    gzip_exit_code = gzip_process.returncode
    if gzip_exit_code != 0:
        return Err(f"VCF filtering encountered an error with exit code {gzip_exit_code}:\n\
                       {gzip_stderr}")

    # remove VCFtools log files
    for file in glob.glob("*.log"):
        os.remove(file)

    new_file = NewFile(new_vcf, os.path.abspath(new_vcf))
    return Ok(new_file)

def collect_filtered_positions(vcf_path: str) -> Result[pl.LazyFrame, str]:
    """
        Use BCFtools to define the positions that are still present in the 
        filtered VCF. This positions will be cross-referenced downstream
        with a table of CADD scores, making it possible to filter a CADD
        table without requiring that table to be annotated with sample IDs.
    
    Args:
        vcf_path: the path to the VCF to be assessed.
    
    Returns:
        Result[pl.LazyFrame, str]: A result containing either a Polars
        LazyFrame of the positions, or an informative error message.
    """

    if not vcf_path.endswith(".vcf.gz") and not vcf_path.endswith(".vcf"):
        return Err("VCF is not in a recognizable format. Please double \
                   check that your VCF file is either gzip-compressed and \
                   ends with '.vcf.gz', or is uncompressed and ends with \
                   '.vcf'.")

    # define the command
    command = ["bcftools", "query", "-f", "%POS\n", vcf_path]

    # Run the command and collect the output
    result = subprocess.run(command,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, check=True)

    # Check if the command was unsuccessful
    exitcode = result.returncode
    if exitcode != 0:
        return Err(f"BCFtools failed to find positions in the provided VCF and\
                   terminated with exit code {exitcode}:\n{result.stderr}")

    # The stdout contains the positions, one per line
    positions = result.stdout.strip().split("\n")

    # Convert to integers and double check that the list isn't empty
    positions = list(map(int, positions))
    if len(positions) == 0:
        return Err("BCFtools encountered a silent error such that the constructed\
                   list of positions is empty")

    # convert the list to a polars dataframe and return it
    position_df = pl.LazyFrame({"Positions": positions})
    return Ok(position_df)

def main() -> None:
    """
    Main handles I/O, coordinates the flow of data through the above functions,
    and exposes any errors.
    """

    # collect command line argument result
    arg_result = parse_command_line_args()
    match arg_result:
        case Ok(result):
            vcf_path = result.vcf
            table_path = result.cadd_table
            animals = result.animal_file
            label = result.label
        case Err(result):
            sys.exit(f"{result}")

    # make sure animal names match the sample IDs in the VCF
    animals_df = pl.read_csv(animals, separator='\t', null_values = ["NA", ""],
                             has_header=False, new_columns=["Sample ID"]).select(
                                pl.col(["Sample ID"])).drop_nulls()
    unify_result = unify_sample_names(animals_df, vcf_path, label)
    match unify_result:
        case Ok(result):
            renamed_vcf_path = result.path
        case Err(result):
            sys.exit(f"{result}")

    # filter the VCF
    filtering_result = filter_vcf_samples(renamed_vcf_path, animals, label)
    match filtering_result:
        case Ok(result):
            filtered_vcf = result.path
        case Err(result):
            sys.exit(f"{result}")

    # define the positions
    position_result = collect_filtered_positions(filtered_vcf)
    match position_result:
        case Ok(result):
            positions = result
        case Err(result):
            sys.exit(f"{result}")

    # lazily read the CADD Score table and filter its rows to only those
    # that contain a mutation that is present in the newly filtered VCF
    cadd_score_df = pl.read_csv(table_path, separator = '\t', skip_rows=1, null_values=["NA"], ignore_errors=True)
    assert "Pos" in cadd_score_df.columns, "Position column named 'Pos' is missing in CADD score file."
    cadd_score_lazy = cadd_score_df.lazy()
    positions.join(
        cadd_score_lazy, left_on="Positions", right_on="Pos", how="left"
    ).sink_csv(f"{label}_filtered_cadd_scores.tsv", separator='\t')

if __name__ == "__main__":
    main()
