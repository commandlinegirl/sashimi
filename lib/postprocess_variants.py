##########################################
# Postprocessing of extracted deletions
##########################################

import subprocess
import pybedtools
import pandas as pd

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

def merge_variants(in_fname, out_fname, merge_distance):
    '''
    Input and output file format:
    Chromosome Start End NumReads TPM Smoothed_TPM
    '''

    # First check if file is empty, if so just create an empty file
    with open(in_fname) as f:
        first_line = f.readline()
        if len(first_line) == 0:
            run_cmd("touch {}".format(out_fname))
            return

    bed = pybedtools.BedTool(in_fname)
    bed.merge(d=merge_distance, c=[4, 5, 6], o=['mean', 'mean', 'mean']).saveas(out_fname)
    
def filter_by_variant_length(in_fname, min_variant_length):
    df = pd.read_csv(in_fname, delimiter='\t',
                     names = ['Chromosome', 'Start', 'End', 'NumReads', 'TPM', 'Smoothed_TPM'])
    df['Length'] = df['End'] - df['Start']
    if min_variant_length > 500:
        df = df[df['Length'] >= min_variant_length]
    return df

def postprocess_variant_candidates(in_fname, out_fname, min_variant_length, merge_distance):
    '''
    Input and output file format:
    Chromosome Start End NumReads TPM Smoothed_TPM
    Output file format:
    Chromosome Start End NumReads TPM Smoothed_TPM Length
    Note: NumReads, TPM, Smoothed_TPM store means of values for records in the merged region
    '''
    fname_tmp = 'tmp' + out_fname
    merge_variants(in_fname, fname_tmp, merge_distance)
    df = filter_by_variant_length(fname_tmp, min_variant_length)
    df.to_csv(out_fname, sep='\t', index=False, header=False)

def map_region_to_detection(result_bed):
    region2event_map = {}
    with open(result_bed, 'r') as f:
        for line in f:
            splat = line.split()
            # Keyed by chrom, start
            region2event_map[(splat[0], splat[1])] = splat[4]
    return region2event_map

def add_call_columns(events, raw_output_fn):
    with open ('tempcolumn', 'w') as f:
        res = []
        with open (raw_output_fn, 'r') as r:
            for line in r:
                spl = line.split()
                if (spl[0], spl[1]) in events:
                    res.append((spl[0], spl[1], str(int(spl[1]) + 500), 1))
                else:
                    res.append((spl[0], spl[1], str(int(spl[1]) + 500), 0))
        print >>f, "\n".join(map(lambda x: str(x[0]) + '\t' + str(x[1]) + '\t' + str(x[2]) + '\t' + str(x[3]), res))
    run_cmd("cut -f4 -d$'\t' tempcolumn | paste {} - | tee salmon_all_sorted_tmp.bed".format(raw_output_fn))
    run_cmd("mv salmon_all_sorted_tmp.bed {}".format(raw_output_fn))

