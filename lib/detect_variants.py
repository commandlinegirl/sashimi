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


def detect_variants():
    # Figure out the het ranges for each chromosome and extract the deletions
    ranges = get_het_ranges(het_range, hom_range, df_smoothed, chromosomes, scaling_factor)

    # Now extract the deletions
    smoothed_tpm_column = 5 # must be 0-indexed
    extract_deletions('salmon_all_sorted.bed', ['salmon_dels_only_ho.bed', 'salmon_dels_only_he.bed'], ranges, neighbors, smoothed_tpm_column)

