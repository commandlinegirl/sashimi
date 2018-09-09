from __future__ import absolute_import, division, print_function
import pandas as pd
import numpy as np
import argparse
import logging
import os
import sys 
import copy 

import smooth
import quantsf2bed
import call_variants
import postprocess_outputs

def limit_to_chromosomes(df, chromosomes):
    '''
    Returns only records that are in chromosomes (if "chromosomes" is not empty)
    '''
    if chromosomes:
        df = df[df['Chromosome'].isin(chromosomes)]
    return df


def subtract_baseline_signal(sample_df, background_signal):
    bgsignal = pd.read_csv(background_signal, delimiter='\t',
                           dtype = {'Chromosome': str, 'Start': np.int32, 'End': np.int32,
                                    'NumReads': np.float64, 'TPM': np.float64, 'Smoothed_TPM': np.float64},
                           names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM', 'Smoothed_TPM'])

    bgsignal = bgsignal[~bgsignal['Chromosome'].isin(['X', 'Y'])]

    # Subtract "background signal" TPM from the HG002 data TPM and add the result as a new column in hg002_df_chrom_1
    subtracted = np.subtract(sample_df['TPM'].values, bgsignal['TPM'].values)
    sample_df['Delta'] = pd.DataFrame(subtracted)
     

def read_salmon_output(quantsf, chromosomes):
    # convert quant.sf to bed format
    inf, outf = quantsf, 'quant.bed'
    quantsf2bed.quantsf2bed(inf, outf)

    # sort by first and second columns (chromosome, start)
    #df_raw_sorted = df_raw.sort_values(['Chromosome', 'Start'], ascending=[True, True])
    quantsf2bed.sort_numerically(outf, 'quant_sorted.bed', [0, 1])
    df_raw = pd.read_csv('quant_sorted.bed', delimiter='\t', header = None,
                         dtype = {'Chromosome': str, 'Start': np.int32, 'End': np.int32,
                                  'NumReads': np.float32, 'TPM': np.float32},
                         names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM'])

    # filter out specific chromosomes for the analysis, if they were provided
    df_chr_filtered = limit_to_chromosomes(df_raw, chromosomes)
    return df_chr_filtered

def main(args):

    ######################################################################
    # Read in and preprocess
    ######################################################################
    df = read_salmon_output(args.quantsf, args.chromosomes)
    smooth_signal = True if ((args.smoothing_window >= 3) or args.smoothing_window and args.smoothing_strategy) else False
    if not smooth_signal:
        col_name = 'TPM'
        mesg = "Skipping smoothing. To smooth the data set smoothing window"
        mesg += " to be larger than 2 and select a smoothing strategy."
        logger.info(mesg)
    else:
        col_name = 'Smoothed_TPM'

    df_smoothed = df.copy()
    df_smoothed[col_name] = smooth.smooth_salmon_output(df, smooth_signal, args.smoothing_window, args.smoothing_strategy)

    ######################################################################
    # Detect variants
    # - Figure out the het ranges for each chromosome and extract the deletions
    # - Extract the deletions (HOM and HET) and save them to separate files
    ######################################################################
    
    del_hom_calls, del_het_calls, dup_calls = \
        call_variants.extract_deletions(df_smoothed, args.del_het_tpm_range, args.del_hom_tpm_range, args.dup_tpm_range, col_name)
    del_hom_calls.to_csv('marked_regions_del_ho.bed', sep='\t', index=None, header=True)
    del_het_calls.to_csv('marked_regions_del_he.bed', sep='\t', index=None, header=True)
    dup_calls.to_csv('marked_regions_dup.bed', sep='\t', index=None, header=True)

    ######################################################################
    # Postprocess data
    ######################################################################
    postprocess_outputs.postprocess('del_ho', del_hom_calls, args.merge_distance, args.min_variant_len, args.blacklist)
    postprocess_outputs.postprocess('del_he', del_het_calls, args.merge_distance, args.min_variant_len, args.blacklist)
    postprocess_outputs.postprocess('dup', dup_calls, args.merge_distance, args.min_variant_len, args.blacklist)


def get_args():
    parser = argparse.ArgumentParser(description="Sashimi, a tool copy number analysis.")
 
    parser.add_argument('quantsf',
                         type=str,
                         help='Path to the output of Salmon (the quant.sf file)')
    parser.add_argument('--chromosomes',
                         nargs='*',
                         type=str,
                         help='Limit the analysis to the list of provided chromosomes')

    # parameters for calling variants
    parser.add_argument('--del_hom_tpm_range',
                         default=[0, 0.01],
                         type=float,
                         nargs='*',
                         help='Expected TPM range for homozygous deletions')
    parser.add_argument('--del_het_tpm_range',
                         nargs='*',
                         type=float,
                         help='Expected TPM range for heterozygous deletions')
    parser.add_argument('--dup_tpm_range',
                         nargs='*',
                         type=float,
                         help='Expected TPM range for duplications')
    parser.add_argument('--smoothing_strategy',
                         default='medianfilter',
                         type=str,
                         choices=["medianfilter", "flat", "hanning", "hamming", "bartlett", "blackman"],
                         help='A smoothing strategy that should be used to normalize the TPM signal. Available strategies are: "medianfilter", "flat", "hanning", "hamming", "bartlett", "blackman".')
    parser.add_argument('--smoothing_window',
                         default=3,
                         type=int,
                         help='The size of the smoothing window. The value must be an odd number. If it is lower than 3 the smoothing is not applied.')

    # parameters for output postprocessing 
    parser.add_argument('--merge_distance',
                         type=int,
                         default=0,
                         help='Max merge distance for adjacent events')
    parser.add_argument('--min_variant_len',
                         type=int,
                         default=1,
                         help='Min length of an event that will be reported')
    parser.add_argument('--blacklist',
                         nargs='*',
                         type=str,
                         help='A list of paths to the bed files with regions that should be excluded from the result file')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)

