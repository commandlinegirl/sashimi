from __future__ import absolute_import, division, print_function
import pandas as pd
import numpy as np
import argparse
import logging
import os
import sys 

import smoothing

import quantsf2bed

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
     

def main(args):

    # convert quant.sf to bed format
    inf, outf = args.quantsf, 'quant.bed'
    quantsf2bed.convert_quantsf_to_bed(inf, outf)

    # sort by first and second columns (chromosome, start)
    #df_raw_sorted = df_raw.sort_values(['Chromosome', 'Start'], ascending=[True, True])
    quantsf2bed.sort_numerically(outf, 'quant_sorted.bed', [0, 1])
    df_raw = pd.read_csv('quant_sorted.bed', delimiter='\t', header = None,
                         dtype = {'Chromosome': str, 'Start': np.int32, 'End': np.int32,
                                  'NumReads': np.float32, 'TPM': np.float32},
                         names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM'])

    # filter out specific chromosomes for the analysis, if they were provided
    df_chr_filtered = limit_to_chromosomes(df_raw, args.chromosomes)
    df_chr_filtered.to_csv('salmon_all_sorted.bed', sep='\t', index=False, header=True)

def get_args():
    parser = argparse.ArgumentParser(description="Sashimi, a tool copy number analysis.")
    parser.add_argument('quantsf',
                         type=str,
                         help='Path to the output of salmon (quant.sf file)')
    parser.add_argument('--chromosomes',
                         nargs='*',
                         type=str,
                         help='List of chromosomes to limit the analysis to')
    parser.add_argument('--blacklist',
                         type=str,
                         help='Path to the bed file with regions that should not be included in the analysis')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)

