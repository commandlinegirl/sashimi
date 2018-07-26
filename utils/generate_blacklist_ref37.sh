#!/usr/bin/env bash

# Generate a set of readpairs from the reference genome (37)
# at 30x coverage (num_reads corresponds to the -n arg
# of mason, which sets the number of pairs)
# https://github.com/seqan/seqan/tree/master/apps/mason2

REFSEQUENCE=file-FG91QGQ0zq64P4Xv4jbZ7V5Y
MASON_APPLET=project-F2Yj92j0B91505Kq3Jk5v2fp:applet-F7k7g0Q0B91K4zbXBg035vvk\

dx run $MASON_APPLET\
   -isequence=$REFSEQUENCE\
   -inum_reads=340000000 \
   -iseed=103 \
   -iread_name_prefix=hom_ref37 \
   --project project-F2Yj92j0B91505Kq3Jk5v2fp \
   -y --brief
