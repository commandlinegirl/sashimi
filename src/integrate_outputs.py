#!/usr/bin/env python
from __future__ import division
import csv
import os, argparse

import subprocess
import pybedtools
import pandas as pd
import numpy as np

def add_output_columns(all_classifier_outputs, event_outputs):
    """
    """
    df_all = pd.read_csv(all_classifier_outputs, delimiter='\t', header = 0)
    for i, fn in enumerate(event_outputs):
        df = pd.read_csv(fn, delimiter='\t', header = 0,
                         dtype = {'Chromosome': str, 'Start': np.int32, 'End': np.int32,
                                  'NumReads': np.float32, 'TPM': np.float32,
                                  'Event': np.int32})
        df_all['Classifier_' + str(i)] = df['Event']
    return df_all 


def integrate_multiple_calls(merged_events_fn, sequence_of_ws, integrated_output_fn):
    summed_ws = sum([int(item) for item in sequence_of_ws])
    assert(summed_ws > 0)
    writable = open(integrated_output_fn, "w")
    writable.write("Chromosome\tStart\tEnd\tNumReads\tTPM\tConfidence\n")
    with open(merged_events_fn, "r") as f:
        contents = csv.reader(f, delimiter='\t')
        next(contents) # skip header
        for row in contents:
            assert(len(sequence_of_ws) == (len(row)-5))
            measure = sum([float(item) * sequence_of_ws[indx] for indx, item in enumerate(row[len(row)-len(sequence_of_ws):])])
            confidence = measure / summed_ws
            call = 1 if measure > 0 else 0
            writable.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(row[0], row[1], row[2], row[3], call, confidence))
    writable.close()


def merge_adjacent_variants(df, merge_distance):
    '''
    Input and output file format:
    Chromosome Start End NumReads TPM Smoothed_TPM
    '''
    bed = pybedtools.BedTool.from_dataframe(df)
    merged_df = bed.merge(d=merge_distance, c=[3, 4, 5], o=['sum', 'mean', 'mean']) \
                   .to_dataframe(names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM', 'Confidence'])
    return merged_df


def filter_by_score(df, min_score):
    df = df[df['Confidence'] >= min_score]
    return df


def filter_by_variant_length(df, min_variant_len):
    df['Length'] = df['End'] - df['Start']
    if min_variant_len > 500:
        df = df[df['Length'] >= min_variant_len]
    return df

def remove_blacklisted_events(salmon_df, blacklist_files):
    '''
    Remove specific regions from the output
    '''
    salmon_output = salmon_df.copy()
    for blacklist_fn in blacklist_files:
        print("Filtering out regions from {}".format(blacklist_fn))
        salmon = pybedtools.BedTool.from_dataframe(salmon_output)
        blacklist = pybedtools.BedTool(blacklist_fn)
        filtered_df = salmon.intersect(blacklist, v=True) \
            .to_dataframe(names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM', 'Confidence', 'Length'], sep='\t')
        salmon_output = filtered_df
    return filtered_df


def postprocess_called_events(zyg_type, call_outputs, args):
    all_calls_in_one_file = 'all_classifier_outputs_{}.bed'.format(zyg_type)
    score_file = "integrated_output_{}.bed".format(zyg_type)
    merged_calls = "result_dels_{}.bed".format(zyg_type)

    salmon_all_outputs = add_output_columns(args.salmon_all_sorted, call_outputs)
    salmon_all_outputs.to_csv(all_calls_in_one_file, sep='\t', index=False, header=True)

    integrate_multiple_calls(all_calls_in_one_file, args.smoothing_windows, score_file)

    integrated_df = pd.read_csv(score_file, delimiter='\t', header=0)
    filtered_df = filter_by_score(integrated_df, args.min_score)

    merged_df = merge_adjacent_variants(filtered_df, args.merge_distance)
    filtered_df = filter_by_variant_length(merged_df, args.min_variant_len)

    if args.blacklist:
        filtered_df = remove_blacklisted_events(filtered_df, args.blacklist)

    filtered_df.to_csv(merged_calls, sep='\t', index=False, header=False)


def main(args):
    if args.call_outputs_ho:
        postprocess_called_events('ho', args.call_outputs_ho, args)
    if args.call_outputs_he:
        postprocess_called_events('he', args.call_outputs_he, args)

def get_args():
    parser = argparse.ArgumentParser(description="Stitch together deletion events from different callers")

    parser.add_argument('--salmon_all_sorted',
                         type=str,
                         help='Path to the output of BED transformation of SALMON output')
    parser.add_argument('--call_outputs_ho',
                         nargs='*',
                         type=str,
                         help='List of HOM events classified; tangent to salmon_all_sorted')
    parser.add_argument('--call_outputs_he',
                         nargs='*',
                         type=str,
                         help='List of HET events classified; tangent to salmon_all_sorted')
    parser.add_argument('--smoothing_windows',
                         nargs='+',
                         type=int,
                         help='Smoothing window sizes')
    parser.add_argument('--merge_distance',
                         type=int,
                         default=0,
                         help='Max merge distance for adjacent events')
    parser.add_argument('--min_variant_len',
                         type=float,
                         default=500,
                         help='Min length of a deletion')
    parser.add_argument('--min_score',
                         type=float,
                         default=1.0,
                         help='Min score for classification as a deletion')
    parser.add_argument('--blacklist',
                         nargs='*',
                         type=str,
                         help='Path to the bed file(s) with regions that should not be included in the final output')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)
