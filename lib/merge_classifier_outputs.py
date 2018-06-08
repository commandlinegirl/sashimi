#!/usr/bin/env python
import os

def main(args):
    # The following line(s) initialize your data object inputs on the platform
    # into dxpy.DXDataObject instances that you can start using immediately.
    # salmon_all_sorted = dxpy.DXFile(salmon_all_sorted)
    # event_outputs = [dxpy.DXFile(item) for item in event_outputs]
    # The following line(s) download your file inputs to the local file system
    # using variable names for the filenames.
    salmon_all_sorted = args.salmon_all_sorted
    event_outputs = args.event_outputs
    # dxpy.download_dxfile(salmon_all_sorted.get_id(), "salmon_all_sorted")
    event_output_hash = {}
    for i, f in enumerate(event_outputs):
        # dxpy.download_dxfile(f.get_id(), "event_outputs-" + str(i))
        fname = "event_outputs-" + str(i)
        with open(fname, 'r') as event_output:
            for event in event_output:
                spl = event.split()
                if not fname in event_output_hash:
                    event_output_hash[fname] = set()
                event_output_hash[fname].add((str(spl[0]), str(spl[1])))

    # Piecewise merge must happen here...
    # Use salmon_all_sorted to track what coordinate system we're on
    writer = open("merged_deletion", "w")
    # Expects chr, start, end,
    with open("salmon_all_sorted", 'r') as f:
        for line in f:
            spl = line.split()
            a = []
            for i, f in enumerate(event_outputs):
                fname = "event_outputs-" + str(i)
                if (str(spl[0]), str(spl[1])) in event_output_hash[fname]:
                    a.append('1')
                else:
                    a.append('0')
            generateRow = "{c}\t{s}\t{e}\t{numReads}\t{TPM}\t{Smoothed_TPM}\t".format(
                                c=spl[0], s=spl[1], e=spl[2], numReads=spl[3], TPM=spl[4], Smoothed_TPM=spl[5]) + "\t".join(a) + "\n"

            writer.write(generateRow)
        writer.close()
    # The following line(s) use the Python bindings to upload your file outputs
    # after you have created them on the local file system.  It assumes that you
    # have used the output field name for the filename for each output, but you
    # can change that behavior to suit your needs.

    # merged_deletion = dxpy.upload_local_file("merged_deletion")

    # The following line fills in some basic dummy output and assumes
    # that you have created variables to represent your output with
    # the same name as your output fields.

    # output = {}
    # output["merged_deletion"] = dxpy.dxlink(merged_deletion)

    # return output

def get_args():
    parser = argparse.ArgumentParser(description="Stitch together deletion events from different callers")

    salmon_all_sorted = args.salmon_all_sorted
    event_outputs = args.event_outputs

    parser.add_argument('--salmon_all_sorted',
                         help='Path to the output of BED transformation of SALMON output')
    parser.add_argument('--event_outputs',
                         help='List of events classified; tangent to salmon_all_sorted')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)
