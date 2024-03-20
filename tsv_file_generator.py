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
        description=(
            'Modify cgppindel vcf file to include VAF,'
        )
    )
    parser.add_argument(
        '-v', '--vcfs', type=str,
        help='Panel filtered VCF(s) from which to generate excel workbook'
    )
    parser.add_argument(
        '-o', '--output_filename',
        help='Output VCF from script'
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
    vcf_df = pd.read_csv(input_vcf, sep="\t", comment='#', compression='infer',
                         names=cols, header=None)

    return vcf_df

example_vcf = read_vcf_df('cgppindel_output/127995080-24037K0005-24NGSHO7-8128-F-96527893_vs_TA2_S59_L008_tumor.flagged.vcf.gz')


def get_field_index(column, field):
    """
    Returns the index of a given field in a column with values separated by ':'

    Args:
        - column (pd.Series): df column to get field index
        - field (str): field to get index of

    Returns: index of given field
    """
    return list(column.apply(lambda x: x.split(':').index(field)))[0]

def get_field_value(column, index):
    """
    Given a column with ':' separated values and index, return a
    series of values with the value split by the index

    Args:
        - column (pd.Series): column to split and return from
        - index (int): index to select

    Returns: series of values selected by index
    """
    return column.apply(lambda x: int(x.split(':')[index]))

NU_index = get_field_index(example_vcf['FORMAT'], 'NU')
NU_value = get_field_value((example_vcf['TUMOUR']), NU_index)

PU_index = get_field_index(example_vcf['FORMAT'], 'PU')
PU_value = get_field_value((example_vcf['TUMOUR']), PU_index)
test = NU_value + PU_value

def generate_annotation_df(vcf_df):
    
    # extract first 5 columnbs for the annotation file
    annotation_df = vcf_df[['CHROM','POS','ID','REF','ALT']].copy()
    
    # get indices of required fields
    pu_index = get_field_index(vcf_df['FORMAT'], 'PU')
    nu_index = get_field_index(vcf_df['FORMAT'], 'NU')
    pr_index = get_field_index(vcf_df['FORMAT'], 'PR')
    nr_index = get_field_index(vcf_df['FORMAT'], 'NR')
    
    # get values for given field from NORMAL and TUMOUR columns
    format_columns = ['NORMAL', 'TUMOUR']
    
    # Create columns for Allele frequency
    for column in format_columns:
        pu_values = get_field_value(vcf_df[column], pu_index)
        nu_values = get_field_value(vcf_df[column], nu_index)
        pr_values = get_field_value(vcf_df[column], pr_index)
        nr_values = get_field_value(vcf_df[column], nr_index)
    
        # calculate af & format as pct to 1dp
        af_values = (pu_values + nu_values) / (pr_values + nr_values)
        af_values = af_values.fillna(0)
        af_pcts = af_values.apply(lambda x: '{:.1f}'.format(float(x * 100)))
        
        annotation_df[f'{column} AF'] = af_pcts
    
    # Create columns for Read depth
    for column in format_columns:
        pr_values = get_field_value(vcf_df[column], pr_index)
        nr_values = get_field_value(vcf_df[column], nr_index)
        
        annotation_df[f'{column} DP'] = pr_values + nr_values

    return annotation_df

example_annotation = generate_annotation_df(example_vcf)

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
    
    
    
if __name__ == '__main__':
    main()
