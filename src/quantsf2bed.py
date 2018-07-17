import sys
import csv

def convert_quantsf_to_bed(in_fname, out_fname, append_chr=False):
    in_file = open(in_fname, 'r')
    out_file = open(out_fname, 'w')
    header = '{}\t{}\t{}\t{}\t{}\n'.format('Chromosome', 'Start', 'End','NumReads', 'TPM')
    reader = csv.DictReader(in_file, delimiter='\t')
    for row in reader:
        s = row['Name']
        fields = s.split(':')
        chrom = fields[0]
        (start,end) = fields[1].split('-')
        out_file.write("{c}\t{s}\t{e}\t{numReads}\t{TPM}\n".format(
            c= "chr" + chrom if append_chr else chrom,
            s=int(start) - 1,
            e=end,
            numReads=row['NumReads'],
            TPM=row['TPM'])
        )
    out_file.close()

def sort_numerically(in_fname, out_fname, col_nums):
    reader = csv.reader(open(in_fname, 'r'), delimiter="\t")
    writer = csv.writer(open(out_fname, 'w'), delimiter='\t')
    for line in sorted(reader, key=lambda x: (x[col_nums[0]], float(x[col_nums[1]]))):
        writer.writerow(line)
