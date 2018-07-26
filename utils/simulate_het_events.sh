#!/usr/bin/env bash

# Generate a sequence with inserted / deleted regions 
# at 30x coverage (num_reads corresponds to the -n arg
# of mason, which sets the number of pairs. In total we 
# want 340e6 reads 150kb long, 170e6 before concatenation)
# https://github.com/seqan/seqan/tree/master/apps/mason2
#
# file-FG9Jfq80gvzQ1xkG1xjqQgQG is the output of 
# analysis-FG9JZP80B9101Z92G74J7gpB.stage-F7k7gg80B91JKvpY13ygKg67.simulated_sequence

dx run project-F2Yj92j0B91505Kq3Jk5v2fp:workflow-F9pQf680B91BXKYz11xgByGg \
   -imason_refgen.sequence=file-FG91QGQ0zq64P4Xv4jbZ7V5Y \
   -imason_refgen.num_reads=170000000 \
   -istage-mason.sequence=file-FG9Jfq80gvzQ1xkG1xjqQgQG \
   -istage-mason.num_reads=170000000 \
   --project project-F2Yj92j0B91505Kq3Jk5v2fp \
   -y --brief
