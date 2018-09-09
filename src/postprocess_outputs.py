#!/usr/bin/env python
from __future__ import division
import csv
import os, argparse
import subprocess
import pybedtools
import pandas as pd
import numpy as np

import logging 
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

RESULT_COLUMNS = ['#Chromosome', 'Start', 'End', 'NumReads', 'TPM', 'Smoothed_TPM', 'Event', 'Length']

def merge_adjacent_variants(df, merge_distance):
    bed = pybedtools.BedTool.from_dataframe(df)
    merged_df = bed.merge(d=merge_distance, c=[4, 5, 6, 7, 8], o=['mean', 'mean', 'mean', 'mean', 'sum']) \
                   .to_dataframe(names = RESULT_COLUMNS)
    return merged_df


def filter_calls(df):
    df = df[df['Event'] == 1]
    return df


def filter_by_variant_length(df, min_variant_len):
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
            .to_dataframe(names = RESULT_COLUMNS, sep='\t')
        salmon_output = filtered_df
    return filtered_df


def postprocess(zyg_type, call_outputs, merge_distance, min_variant_len, blacklist):
    result_file = "all_events_{}_{}.bed".format("del", zyg_type)
    result_df = pd.read_csv(result_file, delimiter='\t', header=0)

    result_df = filter_calls(result_df)

    result_df['Length'] = result_df['End'] - result_df['Start']
    if merge_distance > -1:
        logger.info("About to merge %s variants with max distance %s", zyg_type, merge_distance)
        result_df = merge_adjacent_variants(result_df, merge_distance)

    if min_variant_len:
        logger.info("About to remove %s variants shorter than %s", zyg_type, min_variant_len)
        result_df = filter_by_variant_length(result_df, min_variant_len)

    if blacklist:
        result_df = remove_blacklisted_events(result_df, blacklist)

    merged_calls = "merged_events_{}_{}.bed".format("del", zyg_type)
    result_df.to_csv(merged_calls, sep='\t', index=False, header=True)

