#!/usr/bin/env python3
import numpy as np
from itertools import islice
from collections import defaultdict
from operator import itemgetter
from itertools import groupby
import argparse

def split_line (line):
	ary = line.rstrip().split()
	current_ref, current_pos, kmer, current_rd, current, dwell = ary[0], ary[1], ary[2], ary[3], float (ary[6]), float(ary[8])
	return current_ref, current_pos, kmer,current_rd, current, dwell 
def collapsing_reads(inputfile):	 
	out = open (inputfile+'.collapse','w')
	rd_events = defaultdict(list)
	with open (inputfile,'r') as fh:
		for l in  fh:
			if l.startswith ('contig'):
				continue 
			ref, pos, kmer,rd, c, d = split_line(l) 
			k = f'{ref}\t{pos}\t{kmer}\t{rd}'
			rd_events[k].append (c)			
	for k in sorted(rd_events.keys(), key= lambda x:(x.split()[0], int(x.split()[1]))):
		out.write (f"{k}\t{np.mean(np.mean(rd_events[k]))}\n")
	return inputfile+'.collapse'

def window (s, n=5):
    it = iter (s)
    result = tuple (islice(it,n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def check_consecutive(lst):
    return sorted (lst) == lsit (range(min(lst), max(lst)+1 ))

def get_consecutive(seq):
    '''
    in case, where corrdinates are interrupted
    '''
    return [list (map (itemgetter (1), g)) for _ , g in groupby (enumerate (seq), lambda ix: int (ix[0]) - int(ix[1]) )]

def parse_args():
    """Parse command line arguments passed to script invocation."""
    parser = argparse.ArgumentParser(
        description='\n Generated 5mers and orgnaize all associated features in one line; Input has to be sorted!!!!')
    parser.add_argument('-f','--file', help='input file to be slided; 1st column is chromsome; 2ns colomun is position')
    parser.add_argument('-w','--window_size', type=int, default=5, help='sliding window size / stride length; default is 5')
    return parser.parse_args()

'''
New Input
contig  pos     kmer    readidx Current_mean
18s     0       TATCT   4533    86.71
18s     0       TATCT   4534    109.27
18s     0       TATCT   4537    91.285
contig pos	refKmer	readIndex	Curren
18s     0       TATCT   4967    91.78
'''

def combine_multiple_reads (inf_lst):
	def join_str (t):
		return ":".join ([ str(x) for x in list (t)])
	for x in inf_lst:
		rd, current = zip(*x) #for x in inf_lst)]
		print (rd)
		print (current)
		exit()
	oneliner = "{pos}\t{kmer}\t{read}\t{current_mean}\t{current_median}\t{current_std}".format(
	pos=join_str (set(pos)),
	kmer= join_str (base).replace(':',''),
	read=join_str (rd),
	current_mean=np.mean(current), 
	current_median = np.median (current), 
	current_std = np.std(current,ddof=1)
	)
	return oneliner

def main():
#18s     0       TATCT   4574    92.76666666666667
	args = parse_args()
	header = ''
	base = dict()
	ref_pos = defaultdict(list) # sotre positions on ref 
	ref_pos_rd = defaultdict(list) # sotre reads at ref pos 
	rd_rf_current = dict() # store currnet intensity at ref pos 
	collapsed_file = collapsing_reads (args.file)
	with open (collapsed_file) as fh:
		print ("ref\twindown\tkmer\tread\tmean_current")
		for l in fh:
			if l.startswith ('contig'):
				header = l.strip()
				continue
			ary = l.strip().split()
			ref_pos[ary[0]].append(ary[1])
			k = ary[0] + '\t' + ary[1]
			base[k]=ary[2][0]
			ref_pos_rd[k].append (ary[3])
			rd_rf_current[ary[0]+'\t'+ary[1]+'\t'+ary[3]] = ary[4] 
			
	for ref, v in ref_pos.items():
		v = sorted (map (int, set(v)))
		windows = window (v, args.window_size)
		for w in windows:
			consec = get_consecutive (list(w))
			for win in consec:
				kmer = "".join ([base[ref+'\t'+str(p)] for p in win])
				reads = []
				for p in win:
					reads.extend (ref_pos_rd[ref+'\t'+str(p)])
				for r in reads:
					intensities = []
					for p in win:
						kk = "{}\t{}\t{}".format(ref,p,r)
						intensities.append (rd_rf_current.get (kk,'NA'))
					intensities = ":".join (map(str,intensities))
					print ("{}\t{}\t{}\t{}\t{}".format(ref,":".join (map(str,win)),kmer,r,intensities) )	
if __name__ == '__main__':
    main ()
