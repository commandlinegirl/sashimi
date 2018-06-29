from __future__ import absolute_import, division, print_function
#from builtins import map, zip 
import pandas as pd
import numpy as np
import argparse
import logging
import os
import sys 

import Queue

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
        for chrom in range(1, 23) + ['X', 'Y']:
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


def median(lst):
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return sum(sorted(lst)[quotient - 1:quotient + 1]) / 2.


def panel_decision_deletion(neighbor, tpm_threshold, expected_size):
    if len(neighbor) < expected_size: return False
    med_value = median(neighbor)
    if (within_range(med_value, tpm_threshold)):
        return True
    return False

# Dirty hack to find the TPM median of a chromosome in salmon_all_sorted.bed
def median_finder(chrom, inf):
    cmd = "awk '{{ if($1 == {chrom}) {{ print $0 }} }}' {inf}".format(chrom=chrom, inf=inf)
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    stdout_value = proc.communicate()[0]
    return median([float(item.split('\t')[4].replace('\r', '')) for item in stdout_value.split('\n') if len(item.split('\t')) > 3])


def tpm_mean_finder(chrom, df):
    ch = df[df['Chromosome'] == chrom]
    ch_non_zero = ch[ch['Smoothed_TPM'] > 0]
    if ch_non_zero.empty:
        return 0.0, 0.0
    ch_non_zero_vals = ch_non_zero['Smoothed_TPM']
    base_mean = ch_non_zero_vals.mean() # baseline TPM
    base_std = ch_non_zero_vals.std()
    return base_mean, base_std

def extract_deletions(quant_bed, out_fnames, tpm_thresholds, neighbor_k, tpm_column_used):
    tpm_column = tpm_column_used
    chrom_column = 0
    dels_only_file = open(out_fnames[0], 'w')
    hemizygous_dels_only_file = open(out_fnames[1], 'w')

    homozygous_dels = tpm_thresholds['ho']
    hemizygous_dels = tpm_thresholds['he']

    # prev_neighbors is a queue containing all regions that have already been determined as a del
    # next_neighbors is a queue containing all regions that will still need to be processed
    # curr_object is a
    prev_neighbors = Queue.Queue(maxsize=neighbor_k)
    next_neighbors = Queue.Queue(maxsize=neighbor_k)
    placeholder = Queue.Queue(maxsize=neighbor_k)
    curr_chrom = 0

    with open(quant_bed) as f:
        def refill_next_neighbors():
            # refill next_neighbors with contents of placeholder
            # on transition to processing records in the next chromosome
            while not placeholder.empty():
                out = placeholder.get()
                next_neighbors.put(out)
            while next_neighbors.qsize() < neighbor_k:
                try:
                    line = f.readline()
                    spl = line.split()
                    next_neighbors.put(spl)
                except:
                    pass

        def dump_prev_neighbors():
            # Throw away items in this queue
            while not prev_neighbors.empty():
                out = prev_neighbors.get()

        def resolve_current_object(curr_record, prev_neighbors_tpm, next_neighbors_tpm):
            a_ho = within_range(float(curr_record[tpm_column]), homozygous_dels)
            prev_ho = panel_decision_deletion(prev_neighbors_tpm, homozygous_dels, neighbor_k)
            next_ho = panel_decision_deletion(next_neighbors_tpm, homozygous_dels, neighbor_k)
            a_he = within_range(float(curr_record[tpm_column]), hemizygous_dels[curr_record[chrom_column]])
            prev_he = panel_decision_deletion(prev_neighbors_tpm, hemizygous_dels[curr_record[chrom_column]], neighbor_k)
            next_he = panel_decision_deletion(next_neighbors_tpm, hemizygous_dels[curr_record[chrom_column]], neighbor_k)

            # Predicate: (a and (b or c)) or b or c; this should reduce the amount of reported outliers
            if (a_ho and (prev_ho or next_ho)) or prev_ho or next_ho:
                dels_only_file.write('\t'.join(curr_record) + '\n')
            elif (a_he and (prev_he or next_he)) or prev_he or next_he:
                hemizygous_dels_only_file.write('\t'.join(curr_record) + '\n')

#        line = f.readline() # skip header
        line = f.readline()
        spl = line.split()
        next_neighbors.put(spl)
        curr_chrom = spl[0]
        refill_next_neighbors()

        for line in f:
            spl = line.split()
            if next_neighbors.empty():
                refill_next_neighbors()
                dump_prev_neighbors()
                curr_chrom = spl[0]
            curr_record = next_neighbors.get()
            # For now, place new records in placeholder. Eventually, next_neighbors will run out
            if spl[0] != curr_chrom: placeholder.put(spl)
            else: next_neighbors.put(spl)

            prev_neighbors_tpm = list(map(lambda x: float(x[tpm_column]), list(prev_neighbors.queue)))
            next_neighbors_tpm = list(map(lambda x: float(x[tpm_column]), list(next_neighbors.queue)))
            resolve_current_object(curr_record, prev_neighbors_tpm, next_neighbors_tpm)

            if (prev_neighbors.qsize() >= neighbor_k): out = prev_neighbors.get() #pass thru 'out' to consume the first value
            prev_neighbors.put(curr_record)

        def finish_process_queue():
            while not next_neighbors.empty():
                line = next_neighbors.get()
                prev_neighbors_tpm = list(map(lambda x: float(x[tpm_column]), list(prev_neighbors.queue)))
                next_neighbors_tpm = list(map(lambda x: float(x[tpm_column]), list(next_neighbors.queue)))
                resolve_current_object(line, prev_neighbors_tpm, next_neighbors_tpm)
                if prev_neighbors.qsize() >= neighbor_k: out = prev_neighbors.get()
                prev_neighbors.put(line)

        finish_process_queue()
        dump_prev_neighbors()
        while not placeholder.empty(): next_neighbors.put(placeholder.get())
        finish_process_queue()

    dels_only_file.close()
    hemizygous_dels_only_file.close()


def map_region_to_detection(result_bed):
    region2event_map = {}
    with open(result_bed, 'r') as f:
        for line in f:
            splat = line.split()
            # Keyed by chrom, start
            region2event_map[(splat[0], int(splat[1]))] = float(splat[4])
    return region2event_map


def add_call_columns(events, df_smoothed):
    event_series = []
    for index, row in df_smoothed.iterrows():
        event = 0
        if (row['Chromosome'], row['Start']) in events:
            event = 1
        event_series.append(event)
    df_smoothed['Event'] = event_series
    return df_smoothed

def main(args):

    ######################################################################
    # Variant detection
    # - Figure out the het ranges for each chromosome and extract the deletions
    # - Extract the deletions
    ######################################################################
    df_smoothed = pd.read_csv(args.salmon_all_sorted, delimiter='\t', header = None,
                              dtype = {'Chromosome': str, 'Start': np.int32, 'End': np.int32,
                                       'NumReads': np.float32, 'TPM': np.float32, 'Smoothed_TPM': np.float32},
                              names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM', 'Smoothed_TPM'])

    ranges = get_het_ranges(args.het_range, args.hom_range, df_smoothed, args.chromosomes, args.scaling_factor)
    smoothed_tpm_column = 5 # must be 0-indexed
    extract_deletions(args.salmon_all_sorted, ['salmon_dels_only_ho.bed', 'salmon_dels_only_he.bed'], ranges, args.neighbor_size, smoothed_tpm_column)

    # add a column to "salmon_all_sorted.bed" that marks each region with the detection output:
    # 1 - for variant detected, 0 - no variant
    ho_events = map_region_to_detection('salmon_dels_only_ho.bed')
    df_events = add_call_columns(ho_events, df_smoothed)
    df_events.to_csv('salmon_all_events.bed', sep='\t', index=False, header=True)

def get_args():
    parser = argparse.ArgumentParser(description="Sashimi, a tool copy number analysis.")
    parser.add_argument('salmon_all_sorted',
                         type=str,
                         help='Path to the sorted output of salmon in BED format')
    parser.add_argument('--neighbor_size',
                         default=4,
                         type=int,
                         help='To better help make a decision over whether to call a deletion or insertion - look at what your friends are doing')
    parser.add_argument('--hom_range',
                         default=[0, 0.01],
                         help='Expected range for homozygous deletions')
    parser.add_argument('--het_range',
                         default=[0.01, 5],
                         help='Expected TPM range for heterozygous deletions')
    parser.add_argument('--chromosomes',
                         type=int,
                         help='List of chromosomes to limit the analysis to')
    parser.add_argument('--merge_distance',
                         default=500,
                         type=int,
                         help='List of chromosomes to limit the analysis to')
    parser.add_argument('--scaling_factor',
                         default=500,
                         type=float,
                         help='List of chromosomes to limit the analysis to')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)

