#!/usr/bin/env python
"""Convert Salmon '.sf' output to CNVkit's '.cnr' format."""
from __future__ import absolute_import, division, print_function
import os
import sys

import numpy as np
import pandas as pd
from cnvlib.cnary import CopyNumArray as CNA
from cnvlib.rna import safe_log2
from cnvlib.params import NULL_LOG2_COVERAGE
from skgenome import tabio

READ_LENGTH = 150  # Not super important


def parse_coords(coords):
    chrom, rest = coords.split(':', 1)
    start, end = rest.split('-')
    return chrom, int(start) - 1, int(end)


table = pd.read_table(sys.argv[1])
chroms, starts, ends = zip(*table['Name'].apply(parse_coords))
depths = READ_LENGTH * table['NumReads'] / table['Length']
norm_depth = table['TPM'] / table['TPM'][depths > 0].median()
log2_ratios = safe_log2(norm_depth, NULL_LOG2_COVERAGE)
weights = table['EffectiveLength'] / table['EffectiveLength'].max()

cnarr = CNA.from_columns({
    'chromosome': chroms,
    'start': starts, # np.array(starts) - 1,
    'end': ends,
    'gene': '-',
    'log2': log2_ratios,
    'depth': depths,
    'weight': weights,
    })
cnarr.sort()
tabio.write(cnarr, sys.stdout)
