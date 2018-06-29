#!/usr/bin/env python
import os, argparse
import pandas as pd
import numpy as np

def add_output_columns(salmon_all, event_outputs):
    """
    """
    for i, fn in enumerate(event_outputs):
        df = pd.read_csv(fn, delimiter='\t', header = 1,
                         dtype = {'Chromosome': str, 'Start': np.int32, 'End': np.int32,
                                  'NumReads': np.float32, 'TPM': np.float32, 'Smoothed_TPM': np.float32, 'Event': np.int32},
                         names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM', 'Smoothed_TPM', 'Event'])
        salmon_all['Classifier_' + str(i)] = df['Event']
    return salmon_all

def main(args):
    salmon_all_sorted = args.salmon_all_sorted
    event_outputs = args.event_outputs

    salmon_all = pd.read_csv(args.salmon_all_sorted, delimiter='\t', header = 1,
                     dtype = {'Chromosome': str, 'Start': np.int32, 'End': np.int32,
                              'NumReads': np.float32, 'TPM': np.float32, 'Smoothed_TPM': np.float32},
                     names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM', 'Smoothed_TPM'])

    salmon_all_outputs = add_output_columns(salmon_all, event_outputs)
    salmon_all_outputs.to_csv('merged_events.bed', sep='\t', index=False, header=True)


def get_args():
    parser = argparse.ArgumentParser(description="Stitch together deletion events from different callers")

    # salmon_all_sorted = args.salmon_all_sorted
    # event_outputs = args.event_outputs

    parser.add_argument('--salmon_all_sorted',
                         type=str,
                         help='Path to the output of BED transformation of SALMON output')
    parser.add_argument('--event_outputs',
                         nargs='*',
                         type=str,
                         help='List of events classified; tangent to salmon_all_sorted')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)
