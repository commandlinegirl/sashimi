from __future__ import absolute_import, division, print_function
#from builtins import map, zip 
import pandas as pd
import numpy as np
import argparse
import logging
import os
import sys 

import quantsf2bed
import smoothing
import detect_variants
import postprocess_variants

def limit_to_chromosomes(df, chromosomes):
    '''
    Returns only records that are in chromosomes (if chromosomes is not empty
    '''
    if chromosomes:
        df = df[df['Chromosome'].isin(chromosomes)]
    return df

def main(args):

    ######################################################################
    # Preprocessing
    ######################################################################

    # convert quant.sf to bed format
    inf, outf = args.salmon_output, 'quant.bed'
    quantsf2bed.convert_quantsf_to_bed(inf, outf)

    # sort by first and second columns (chromosome, start)
    #df_raw_sorted = df_raw.sort_values(['Chromosome', 'Start'], ascending=[True, True])
    quantsf2bed.sort_numerically('quant.bed', 'quant_sorted.bed', [0,1])
    df_raw = pd.read_csv('quant_sorted.bed', delimiter='\t', header = None,
                         dtype = {'Chromosome': str, 'Start': np.int32, 'End': np.int32,
                                  'NumReads': np.float32, 'TPM': np.float32},
                         names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM'])

    # filter out specific chromosomes for the analysis, if they were provided
    df_chr_filtered = limit_to_chromosomes(df_raw, args.chromosomes)

    # append a 'Smoothed_TPM' column and save
    df_smoothed = smoothing.smooth_salmon_output(df_chr_filtered, args.smooth_raw_output, args.smoothing_window_len, args.smoothing_strategy)
    df_smoothed.to_csv('salmon_all_sorted.bed', sep='\t', index=False, header=False)

    ######################################################################
    # Variant detection
    # - Figure out the het ranges for each chromosome and extract the deletions
    # - Extract the deletions
    ######################################################################

    ranges = detect_variants.get_het_ranges(args.het_range, args.hom_range, df_smoothed, args.chromosomes, args.scaling_factor)
    smoothed_tpm_column = 5 # must be 0-indexed
    detect_variants.extract_deletions('salmon_all_sorted.bed', ['salmon_dels_only_ho.bed', 'salmon_dels_only_he.bed'], ranges, args.neighbors, smoothed_tpm_column)

    ######################################################################
    # Post-processing:
    # - merges adjacent deletions
    # - adds length column
    # - removes short variants
    ######################################################################

    postprocess_variants.postprocess_variant_candidates('salmon_dels_only_ho.bed', 'salmon_merged_dels_ho.bed', args.min_variant_length, args.merge_distance)
    postprocess_variants.postprocess_variant_candidates('salmon_dels_only_he.bed', 'salmon_merged_dels_he.bed', args.min_variant_length, args.merge_distance)

    # add a column to "salmon_all_sorted.bed" that marks each region with the detection output:
    # 1 - for variant detected, 0 - no variant
    #ho_events = postprocess_variants.map_region_to_detection('salmon_dels_only_ho.bed')
    he_events = postprocess_variants.map_region_to_detection('salmon_dels_only_he.bed')
    #postprocess_variants.add_call_columns(ho_events, 'salmon_all_sorted.bed')
    postprocess_variants.add_call_columns(he_events, 'salmon_all_sorted.bed')

def get_args():
    parser = argparse.ArgumentParser(description="Sashimi, a tool copy number analysis.")
    parser.add_argument('salmon_output',
                         help='Path to the output of salmon (quant.sf file)')
    parser.add_argument('--blacklist',
                         help='Path to the bed file with regions that should not be included in the analysis')
    parser.add_argument('--chromosomes',
                         help='List of chromosomes to limit the analysis to')
    parser.add_argument('--smooth_raw_output',
                         action='store_true',
                         help='Apply smothing to the data')
    parser.add_argument('--smoothing_window_len',
                         default=5,
                         help='Size of the smoothing window')
    parser.add_argument('--smoothing_strategy',
                         default='hanning',
                         help='What smoothing strategy should be used')
    parser.add_argument('--neighbors',
                         default=4,
                         help='To better help make a decision over whether to call a deletion or insertion - look at what your friends are doing')
    parser.add_argument('--hom_range',
                         default=[0, 0.01],
                         help='Expected range for homozygous deletions')
    parser.add_argument('--het_range',
                         default=[0.01, 5],
                         help='Expected TPM range for heterozygous deletions')
    parser.add_argument('--min_variant_length',
                         default=500,
                         help='List of chromosomes to limit the analysis to')
    parser.add_argument('--merge_distance',
                         default=500,
                         help='List of chromosomes to limit the analysis to')
    parser.add_argument('--scaling_factor',
                         default=1.5,
                         help='List of chromosomes to limit the analysis to')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)

