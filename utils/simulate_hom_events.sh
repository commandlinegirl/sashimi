#!/usr/bin/env bash

# Generate a sequence with inserted / deleted regions 
# at 30x coverage (num_reads corresponds to the -n arg
# of mason, # which sets the number of pairs)
# https://github.com/seqan/seqan/tree/master/apps/mason2

dx run project-F2Yj92j0B91505Kq3Jk5v2fp:workflow-F9GYB3Q0B91FbJ76JYZZf7Jy \
   -istage-F7k7gg80B91JKvpY13ygKg67.reference_sequence=file-FG91QGQ0zq64P4Xv4jbZ7V5Y \
   -istage-F7k7gg80B91JKvpY13ygKg67.prefix=wgsimulated_0.5k_to_30k \
   -istage-F7k7gg80B91JKvpY13ygKg67.DUPLICATION_minimum_length=500 \
   -istage-F7k7gg80B91JKvpY13ygKg67.DUPLICATION_maximum_length=30000 \
   -istage-F7k7gg80B91JKvpY13ygKg67.DUPLICATION_number=1000 \
   -istage-F7k7gg80B91JKvpY13ygKg67.INDEL_minimum_length=500 \
   -istage-F7k7gg80B91JKvpY13ygKg67.INDEL_maximum_length=30000 \
   -istage-F7k7gg80B91JKvpY13ygKg67.INDEL_number=1000 \
   -istage-F7k7gp00B9180ZV7146KKvjf.num_reads=340000000 \
   -istage-F7k7gp00B9180ZV7146KKvjf.seed=100 \
   -istage-F7k7gp00B9180ZV7146KKvjf.read_name_prefix=hom \
   --project project-F2Yj92j0B91505Kq3Jk5v2fp \
   -y --brief
