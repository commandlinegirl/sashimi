from __future__ import absolute_import, division, print_function
import pandas as pd
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import logging
import os
import sys 
import copy 
import smoothing
import quantsf2bed
import read_salmon_output
import call_variants
import integrate_outputs
import os

def main(args):

    # read in and preprocess
    read_salmon_output.main(args)

    # call variants
    caller_outputs_ho = []
    caller_outputs_he = []
    for strategy, window in zip(args.smoothing_strategies, args.smoothing_windows):
        callv_args = argparse.Namespace(**vars(args)) 
        callv_args.smoothing_strategy = strategy
        callv_args.smoothing_window = window
        callv_args.salmon_all_sorted = "salmon_all_sorted.bed"
        call_variants.main(callv_args)
        new_fn = "salmon_all_events_ho_{}_{}.bed".format(strategy, window)
        os.rename("salmon_all_events_ho.bed", new_fn)
        caller_outputs_ho.append(new_fn)
        new_fn = "salmon_all_events_he_{}_{}.bed".format(strategy, window)
        os.rename("salmon_all_events_he.bed", new_fn)
        caller_outputs_he.append(new_fn)

    # integrated outputs from different callers
    integr_args = argparse.Namespace(**vars(args))
    integr_args.call_outputs_ho = caller_outputs_ho
    integr_args.call_outputs_he = caller_outputs_he
    integr_args.salmon_all_sorted = "salmon_all_sorted.bed"
    integrate_outputs.main(integr_args)

def get_args():
    parser = argparse.ArgumentParser(description="Sashimi, a tool copy number analysis.", formatter_class=RawTextHelpFormatter)
    # 
    parser.add_argument('quantsf',
                         type=str,
                         help='Path to the output of salmon (quant.sf file)')
    parser.add_argument('--chromosomes',
                         nargs='*',
                         type=str,
                         help='List of chromosomes to limit the analysis to')

    # call variants
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
    parser.add_argument('--smoothing_strategies',
                         nargs='*',
                         default=['medianfilter'],
                         type=str,
                         choices=["medianfilter", "flat", "hanning", "hamming", "bartlett", "blackman"],
                         help='List of smoothing strategies that should be used. List size must be equal to the list of smoothing_windows argument. Each position in the list corresponds to each position in the smoothing windows list. Available strategies are: "medianfilter", "flat", "hanning", "hamming", "bartlett", "blackman". See numpy.<strategy> for more details on each strategy')
    parser.add_argument('--smoothing_windows',
                         nargs='*',
                         default=[3],
                         type=int,
                         help='List of the sizes of the smoothing window. Each position in the list corresponds to each position in the smoothing_strategies list. The value must be an odd integer. If the value is lower than 3 the smoothing is not applied.')

    # integrate outputs
    parser.add_argument('--merge_distance',
                         type=int,
                         default=500,
                         help='Max merge distance for adjacent events')
    parser.add_argument('--min_variant_len',
                         type=int,
                         default=1000,
                         help='Min length of a deletion that will be reported')
    parser.add_argument('--min_score',
                         type=float,
                         default=1.0,
                         help='Min score (confidence level) for classification as a deletion. Must be a float between 0 and 1')
    parser.add_argument('--blacklist',
                         nargs='*',
                         type=str,
                         help='Path to the list of bed files with regions that should be excluded from the analysis')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)

