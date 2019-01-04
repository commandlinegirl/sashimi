import argparse
import logging
import pandas as pd

#chromosome	start	end	gene	log2	depth	probes	weight
#1	0	10000	-	-20	0	20	11.9972
#1	10000	54500	-	-0.674517	3.5489	89	67.0479
#1	54500	71000	-	-5.82593	1.38414	33	24.5069

def main(args):

    df = pd.read_csv(args.cns, delimiter='\t', header=0)
    df = df[['chromosome', 'start', 'end', 'log2']]
    df = df[df['log2'] < 0]
    df.to_csv('cns_result.bed', sep='\t', index=False, header=True) 

def get_args():
    parser = argparse.ArgumentParser(description="")
 
    parser.add_argument('cns',
                         type=str,
                         help='Path to the cns file')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    main(args)

