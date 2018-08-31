from __future__ import absolute_import, division, print_function
#from builtins import map, zip 
import pandas as pd
import numpy as np
import argparse
import logging
import os
import sys 
import logging

import smoothing

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_het_ranges(het_range, hom_range, df_smoothed, chromosomes, scaling_factor):
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


def within_range(val, range):
    # range refers to a list of two elements [start, end]
    # returns true if val is within range, inclusively
    return val >= range[0] and val <= range[1]


def tpm_mean_finder(chrom, df):
    ch = df[df['Chromosome'] == chrom]
    ch_non_zero = ch[ch['Smoothed_TPM'] > 0]
    if ch_non_zero.empty:
        return 0.0, 0.0
    ch_non_zero_vals = ch_non_zero['Smoothed_TPM']
    base_mean = ch_non_zero_vals.mean() # baseline TPM
    base_std = ch_non_zero_vals.std()
    return base_mean, base_std


def extract_deletions(df, col_name, ranges):
    r = ranges['ho'] #TODO
    df_hom = df.copy()
    df_hom['Event'] = (df_hom[col_name] >= r[0]) & (df_hom[col_name] <= r[1])
    df_hom['Event'] = df_hom['Event'].astype(int)
    df_het = df_hom.copy()
    return df_hom, df_het


def get_column_name(strategy, window):
    return "Smoothed_{}_{}".format(strategy, window)


def main(args):

    ######################################################################
    # Smooth input data
    ######################################################################
    df = pd.read_csv(args.salmon_all_sorted, delimiter='\t', header = 0,
                     dtype = {'Chromosome': str, 'Start': np.int32, 'End': np.int32,
                              'NumReads': np.float32, 'TPM': np.float32},
                     names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM'])

    smooth_signal = True if ((args.smoothing_window >= 3) or args.smoothing_window and args.smoothing_strategy) else False
    if not smooth_signal:
        mesg = "Skipping smoothing. To smooth the data set smoothing window"
        mesg += " to be larger than 2 and select a smoothing strategy."
        logger.info(mesg)

    df_smoothed = df.copy()
    col_name = get_column_name(args.smoothing_strategy, args.smoothing_window)
    df_smoothed[col_name] = smoothing.smooth_salmon_output(df, smooth_signal, args.smoothing_window, args.smoothing_strategy)
    df_smoothed.to_csv('salmon_all_smoothed.bed', sep='\t', index=False, header=True)

    ######################################################################
    # Variant detection
    # - Figure out the het ranges for each chromosome and extract the deletions
    # - Extract the deletions (HOM and HET) and save them to separate files
    ######################################################################
    ranges = get_het_ranges(args.het_range, args.hom_range, df_smoothed, args.chromosomes, args.scaling_factor)
    hom_extr, het_extr = extract_deletions(df_smoothed, col_name, ranges)
    hom_extr.to_csv('salmon_all_events_ho.bed', sep='\t', index=None, header=True)
    het_extr.to_csv('salmon_all_events_he.bed', sep='\t', index=None, header=True)

    logger.info('Done')

def get_args():
    parser = argparse.ArgumentParser(description="Sashimi, a tool copy number analysis.")
    parser.add_argument('salmon_all_sorted',
                         type=str,
                         help='Path to the sorted output of salmon in BED format')
    parser.add_argument('--chromosomes',
                         help='List of chromosomes to limit the analysis to')
    parser.add_argument('--hom_range',
                         default=[0, 0.001],
                         help='Expected range for homozygous deletions')
    parser.add_argument('--het_range',
                         default=[0.01, 5],
                         help='Expected TPM range for heterozygous deletions')
    parser.add_argument('--scaling_factor',
                         default=1.5,
                         type=float,
                         help='Scaling factor for heterozygous deletion detection')
    parser.add_argument('--smoothing_strategy',
                         default='hanning',
                         type=str,
                         help='What smoothing strategy should be used')
    parser.add_argument('--smoothing_window',
                         default=5,
                         type=int,
                         help='Size of the smoothing window')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)

