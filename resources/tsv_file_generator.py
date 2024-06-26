#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 09:30:56 2024

@author: arun
"""
import argparse
import csv
from pysam import VariantFile


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
        "--vcf_file",
        type=str,
        help="cgppindel VCF",
    )

    args = parser.parse_args()

    return args


def generate_annotation_data(vcf_file):
    """
    Creates an list of lists from the given vcf file, each list has the
    following elements in the specified order
    
    [CHROM, POS, ID, REF, ALT, NORMAL AF, TUMOUR AF, NORMAL DP, TUMOUR DP]


    Parameters
    ----------
    vcf_file : VariantFile
        dataframe of all variants from a vcf

    Returns
    -------
    annotation_data : list
        Output list of lists with similar structure above

    """
    annotation_data = []

    for variant in vcf_file.fetch():
        # Setting up normal and tumour samples
        normal = variant.samples['NORMAL']
        tumour = variant.samples['TUMOUR']

        # Setting up Variant Allele Frequency
        try:
            normal_af = ((normal['PU'] + normal['NU']) /
                         (normal['PR'] + normal['NR']))
        except ZeroDivisionError:
            normal_af = 0

        try:
            tumour_af = ((tumour['PU'] + tumour['NU']) /
                         (tumour['PR'] + tumour['NR']))
        except ZeroDivisionError:
            tumour_af = 0

        # Setting up Read Depth
        normal_dp = normal['PR'] + normal['NR']
        tumour_dp = tumour['PR'] + tumour['NR']

        # Generate vcf dataframe
        annotation_data.append([variant.chrom, variant.pos, variant.id,
                     variant.ref, ", ".join(variant.alts),
                     normal_af, tumour_af, normal_dp, tumour_dp])

    return annotation_data


def main():
    """
    Main entry points to run the script. Generates the annots.tsv file

    Returns
    -------
    None.
    """
    args = parse_args()
    vcf_file = VariantFile(args.vcf_file)
    annotation_data = generate_annotation_data(vcf_file)
    with open('annots.tsv', 'w', newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(annotation_data)


if __name__ == "__main__":
    main()
