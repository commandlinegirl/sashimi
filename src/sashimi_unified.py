from __future__ import absolute_import, division, print_function
import pandas as pd
import numpy as np
import argparse
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
    parser = argparse.ArgumentParser(description="Sashimi, a tool copy number analysis.")
    # 
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

    # call variants
    parser.add_argument('--hom_range',
                         default=[0, 0.01],
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
                         help='What smoothing strategy should be used')
    parser.add_argument('--smoothing_windows',
                         nargs='*',
                         default=[3],
                         type=int,
                         help='Size of the smoothing window')

    # integrate outputs
    parser.add_argument('--merge_distance',
                         type=int,
                         default=500,
                         help='Max merge distance for adjacent events')
    parser.add_argument('--min_variant_len',
                         type=int,
                         default=1000,
                         help='Min length of a deletion')
    parser.add_argument('--min_score',
                         type=float,
                         default=1.0,
                         help='Min score for classification as a deletion')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)

