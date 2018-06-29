import os
import subprocess
import argparse

def run_cmd(cmd, returnOutput=False):
    print cmd
    if returnOutput:
        try:
            output = subprocess.check_output(cmd, shell=True, executable='/bin/bash').strip()
            return output
        except subprocess.CalledProcessError as e:
            print "error code", e.returncode, e.output
            return
    else:
        try:
            subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        except subprocess.CalledProcessError as e:
            print "error code", e.returncode, e.output


def get_overlaps(in_a, in_b, out_fname, F_val=None, f_val=None, r=False):
    cmd = "bedtools intersect -wao -a {a} -b {b}".format(a=in_a, b=in_b)
    if F_val:
        cmd += " -F {}".format(F_val)
    if f_val:
        cmd += " -f {}".format(f_val)
    if r:
        cmd += " -r"
    cmd += " > {outf}".format(outf=out_fname)
    run_cmd(cmd)


def get_tp(overlap_fname, out_fname):
    cmd = 'cat {overlap} | grep -v \"\-1\" '.format(overlap=overlap_fname)
    cmd += ' > {outf}'.format(outf=out_fname)
    run_cmd(cmd)


def get_delta(overlap_fname, out_fname):
    cmd = 'cat {overlap} | grep \"\-1\" '.format(overlap=overlap_fname)
    cmd += ' > {outf}'.format(outf=out_fname)
    run_cmd(cmd)


def count_lines(f_name):
    return sum(1 for line in open(f_name))


def main(args):
    #TODO: for now dont fully expose these options
    F_val = args.overlap
    f_val = args.overlap

    # note: using both -F and -f options is equivalent to using both -r and -F
    # using all -f, -F, -r will throw an error
    r = False 

    sample = args.sample_events
    truth = args.truth_events

    get_overlaps(sample, truth, 'overlaps_sample_side.bed', F_val, f_val, r)
    get_overlaps(truth, sample, 'overlaps_truth_side.bed', F_val, f_val, r)

    get_tp('overlaps_sample_side.bed', 'true_positives.bed')
    get_delta('overlaps_sample_side.bed', 'false_positives.bed')
    get_delta('overlaps_truth_side.bed', 'false_negatives.bed')

    tp = count_lines('true_positives.bed')
    fp = count_lines('false_positives.bed')
    fn = count_lines('false_negatives.bed')

    precision = 1.0 * tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = 1.0 * tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f_score = 2.0 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    # save the outputs, one per file (we save them in files, an output
    # per line so that it can be easily read by WDL)
    with open("score_precision.txt", "w") as pscore:
        pscore.write("{}\n".format(precision))
    with open("score_recall.txt", "w") as rscore:
        rscore.write("{}\n".format(recall))
    with open("score_f_score.txt", "w") as fscore:
        fscore.write("{}\n".format(f_score))


def get_args():
    parser = argparse.ArgumentParser(description="Integrate multiple call outputs into one result")
    parser.add_argument('--sample_events',
                         type=str,
                         help='Path to the sample bed file')
    parser.add_argument('--truth_events',
                         type=str,
                         help='Path to the truth bed file')
    parser.add_argument('--overlap',
                         type=float,
                         default=0.5,
                         help='Minimum reciprocal overlap required')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)
