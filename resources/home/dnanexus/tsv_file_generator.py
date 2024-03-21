#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 09:30:56 2024

@author: arun
"""

import argparse
import pandas as pd


def parse_args():
    """
    Parse command line args

    Args: None

    Returns:
        - args (Namespace): object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=("Modify cgppindel vcf file to include VAF,")
    )
    parser.add_argument(
        "-v",
        "--vcfs",
        type=str,
        help="Panel filtered VCF(s) from which to generate excel workbook",
    )

    args = parser.parse_args()

    return args


def read_vcf_df(input_vcf):
    """
    Reads vcf into pandas df, returns header as a list for output (bsvi) vcf

    Args:
        - input_vcf (file): vcf file to read in

    Returns:
        - vcf_df (df): df of variants from vcf
    """
    cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
            "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOUR"]

    # read vcf records into df
    vcf_df = pd.read_csv(input_vcf, sep="\t", comment="#", compression="infer",
                         names=cols, header=None)
    return vcf_df


def get_field_index(column, field):
    """
    Returns the index of a given field in a column with values separated by ':'

    Args:
        - column (pd.Series): df column to get field index
        - field (str): field to get index of

    Returns: index of given field
    """
    return list(column.apply(lambda x: x.split(":").index(field)))[0]


def get_field_value(column, index):
    """
    Given a column with ':' separated values and index, return a
    series of values with the value split by the index

    Args:
        - column (pd.Series): column to split and return from
        - index (int): index to select

    Returns: series of values selected by index
    """
    return column.apply(lambda x: int(x.split(":")[index]))


def generate_annotation_df(vcf_df):
    """
    Creates an annotation df from the given vcf file, structure of the df
    should look like this:

   ----------------------------------------------------------------------------
   CHROM | POS | ID | REF | ALT | NORMAL AF | TUMOUR AF | NORMAL DP | TUMOUR DP
   ----------------------------------------------------------------------------
   chr1  |12345|450f| AGT | TCA |     0     |  0.9999   |     0     |    6
   ----------------------------------------------------------------------------
   chr1  |54321|219b| AC  | TCA |   0.53    |   0.64    |     1     |    90
   ----------------------------------------------------------------------------
   chr2  |12345|699m| TA  |  G  |   0.666   |  0.8432   |    243    |    614
   ----------------------------------------------------------------------------

    Parameters
    ----------
    vcf_df : TYPE
        dataframe of all variants from a vcf

    Returns
    -------
    annotation_df : pd.DataFrame
        Output dataframe with simmilar structure above

    """

    # extract first 5 columnbs for the annotation file
    annotation_df = vcf_df[["CHROM", "POS", "ID", "REF", "ALT"]].copy()

    # get indices of required fields
    pu_index = get_field_index(vcf_df["FORMAT"], "PU")
    nu_index = get_field_index(vcf_df["FORMAT"], "NU")
    pr_index = get_field_index(vcf_df["FORMAT"], "PR")
    nr_index = get_field_index(vcf_df["FORMAT"], "NR")

    # get values for given field from NORMAL and TUMOUR columns
    format_columns = ["NORMAL", "TUMOUR"]

    # Create columns for Allele frequency
    for column in format_columns:
        pu_values = get_field_value(vcf_df[column], pu_index)
        nu_values = get_field_value(vcf_df[column], nu_index)
        pr_values = get_field_value(vcf_df[column], pr_index)
        nr_values = get_field_value(vcf_df[column], nr_index)

        # calculate af & format as pct to 1dp
        af_values = (pu_values + nu_values) / (pr_values + nr_values)
        af_values = af_values.fillna(0)

        annotation_df[f"{column} AF"] = af_values

    # Create columns for Read depth
    for column in format_columns:
        pr_values = get_field_value(vcf_df[column], pr_index)
        nr_values = get_field_value(vcf_df[column], nr_index)

        annotation_df[f"{column} DP"] = pr_values + nr_values

    return annotation_df


def main():
    """
    Main entry points to run the script. Generates the annots.tsv file

    Returns
    -------
    None.
    """
    args = parse_args()
    vcf_df = read_vcf_df(args.vcfs)
    annotation_df = generate_annotation_df(vcf_df)
    annotation_df.to_csv("annots.tsv", sep="\t", header=None, index=False)


if __name__ == "__main__":
    main()
