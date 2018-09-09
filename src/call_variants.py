from __future__ import absolute_import, division, print_function
#from builtins import map, zip 
import pandas as pd
import numpy as np
import argparse
import logging
import os
import sys 
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_het_tpm_ranges(het_range, hom_range, df_smoothed, chromosomes):
    #TODO: remove this maybe
    if het_range is None:
        # convert to dict
        het_range = {}
        if chromosomes:
            chrs = chromosomes
        else:
            chrs = range(1, 23) + ['X', 'Y']
        for chrom in chrs:
            #chrom_med = median_finder(str(chrom), 'salmon_all_sorted.bed')
            base_mean, base_std = tpm_mean_finder(str(chrom), df_smoothed)
            het_mean = base_mean / 2.0
            print("{}: base mean TPM: {}, base stddev: {}, het TPM {}".format(chrom, base_mean, base_std, het_mean))
            # TODO: is it the right way to go?
            scaling_factor = 1.5
            het_range_boundary = base_std / scaling_factor
            het_range[str(chrom)] = [max(het_mean - het_range_boundary, 0), het_mean + het_range_boundary]
    else:
        start, end = het_range[0], het_range[1]
        het_range = {}
        for chrom in list(range(1, 23)) + ['X', 'Y']:
            het_range[str(chrom)] = [start, end]

    # hom_range: [k, l], k<=l
    # het_range: {1: [k, l], 2: [x, y], ..., X: [a, b], Y: [c, d]}
    ranges = { 'ho': hom_range, 'he': het_range}
    print(ranges) # Log these ranges
    return ranges


def tpm_mean_finder(df, col_name):
    ch_non_zero = df[df[col_name] > 0]
    if ch_non_zero.empty:
        return 0.0, 0.0
    ch_non_zero_vals = ch_non_zero[col_name]
    base_mean = ch_non_zero_vals.mean() # baseline TPM
    base_std = ch_non_zero_vals.std()
    return base_mean, base_std


def get_het_tpm_range(df, base_mean, base_std):
    range_boundary = base_mean * 0.1
    het_tpm_range = [max(base_mean - range_boundary, 0), base_mean + range_boundary]
    return het_tpm_range


def get_dup_tpm_range(df, base_mean, base_std):
    dup_tpm_range = [base_mean + 2 * base_std, 10e6]
    return dup_tpm_range


def extract_deletions(df, del_het_range, del_hom_range, dup_range, col_name):
    if not del_hom_range:
        del_hom_range = [0, 0.01]
    
    base_mean, base_std = tpm_mean_finder(df, col_name) 
    logger.info("Mean TPM %s, stddev %s", base_mean, base_std)
    if not del_het_range:
        del_het_range = get_het_tpm_range(df, base_mean, base_std)
    if not dup_range:
        dup_range = get_dup_tpm_range(df, base_mean, base_std)

    logger.info("Using range %s for homozygous deletion detection", del_hom_range)
    df_hom = df.copy()
    df_hom['Event'] = (df_hom[col_name] >= del_hom_range[0]) & (df_hom[col_name] <= del_hom_range[1])
    df_hom['Event'] = df_hom['Event'].astype(int)

    logger.info("Using range %s for heterozygous deletion detection", del_het_range)
    df_het = df.copy()
    df_het['Event'] = (df_het[col_name] >= del_het_range[0]) & (df_het[col_name] <= del_het_range[1])
    df_het['Event'] = df_het['Event'].astype(int)

    logger.info("Using range %s for duplication detection", dup_range)
    df_dup = df.copy()
    df_dup['Event'] = (df_dup[col_name] >= dup_range[0]) & (df_dup[col_name] <= dup_range[1])
    df_dup['Event'] = df_dup['Event'].astype(int)

    return df_hom, df_het, df_dup


