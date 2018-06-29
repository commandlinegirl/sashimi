from __future__ import division
import csv
import argparse

import subprocess
import pybedtools
import pandas as pd


def integrate_multiple_calls(sequence_of_ws, integrated_output_fn):
    summed_ws = sum([int(item) for item in sequence_of_ws])
    assert(summed_ws > 0)
    writable = open(integrated_output_fn, "w")
    with open(args.merged_events, "r") as f:
        contents = csv.reader(f, delimiter='\t')
        next(contents) # skip header
        for row in contents:
            assert(len(sequence_of_ws) == (len(row)-6))
            measure = sum([float(item) * sequence_of_ws[indx] for indx, item in enumerate(row[len(row)-len(sequence_of_ws):])])
            confidence = measure / summed_ws
            call = 1 if measure > 0 else 0
            writable.write("{}\t{}\t{}\t{}\t{}\n".format(row[0], row[1], row[2], call, confidence))
    writable.close()

#TODO: set this programmatically
#bedtools_location = '/Users/azalcman/Downloads/bedtools2/bin'
#pybedtools.helpers.set_bedtools_path(bedtools_location)


def run_cmd(cmd, returnOutput=False):
    print(cmd)
    if returnOutput:
        output = subprocess.check_output(cmd, shell=True).strip()
        return output
    else:
        subprocess.check_call(cmd, shell=True)

def merge_adjacent_variants(in_fname, out_fname, merge_distance, min_score):
    '''
    Input and output file format:
    Chromosome Start End NumReads TPM Smoothed_TPM
    '''
    df = pd.read_csv(in_fname, delimiter='\t', header=None,
                     names = ['Chromosome', 'Start', 'End', 'EventDetected', 'ConfidenceLevel'])

    df_filtered = df[df['ConfidenceLevel'] >= min_score]
    # df_filtered['Length'] = df_filtered['End'] - df_filtered['Start']
    df_filtered['Length'] = df_filtered['End'] - df_filtered['Start']

    print(df_filtered.head())
    
    bed = pybedtools.BedTool.from_dataframe(df_filtered)
    merged = bed.merge(d=merge_distance, c=[6], o=['sum']).saveas(out_fname)
    

# def filter_by_variant_length(in_fname, min_variant_length):
#     df = pd.read_csv(in_fname, delimiter='\t',
#                      names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM', 'Smoothed_TPM'])
#     df['Length'] = df['End'] - df['Start']
#     if min_variant_length > 500:
#         df = df[df['Length'] >= min_variant_length]
#     return df


# def postprocess_variant_candidates(in_fname, out_fname, min_variant_length, merge_distance):
#     '''
#     Input and output file format:
#     Chromosome Start End NumReads TPM Smoothed_TPM
#     Output file format:
#     Chromosome Start End NumReads TPM Smoothed_TPM Length
#     Note: NumReads, TPM, Smoothed_TPM store means of values for records in the merged region
#     '''
#     fname_tmp = 'tmp' + out_fname
#     merge_variants(in_fname, fname_tmp, merge_distance)
#     df = filter_by_variant_length(fname_tmp, min_variant_length)
#     df.to_csv(out_fname, sep='\t', index=False, header=False)


def main(args):
    integrate_multiple_calls(args.neighbor_sizes, "integrated_output.bed")
    merge_adjacent_variants("integrated_output.bed", "merged_calls_ho.bed", args.merge_distance, args.min_score)
    

def get_args():
    parser = argparse.ArgumentParser(description="Integrate multiple call outputs into one result")
    parser.add_argument('--merged_events',
                         type=str,
                         help='Path to the output of multiple classifier calls')
    parser.add_argument('--neighbor_sizes',
                         nargs='+',
                         type=int,
                         help='Neighbor sizes')
    parser.add_argument('--merge_distance',
                         type=int,
                         default=0,
                         help='Max merge distance for adjacent events')
    parser.add_argument('--min_score',
                         type=float,
                         default=1,
                         help='Min score for classification as a deletion')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)

