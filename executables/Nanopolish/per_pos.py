import sys
import pandas as pd
infile = sys.argv[1]
inp=pd.read_csv(infile,sep='\t') #this step will take ages if the file is huge!
#For the simple figure of (per median/position) 
grouped_multiple_mean_inp = inp.groupby(['contig', 'position','reference_kmer']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_inp = grouped_multiple_mean_inp.reset_index()
grouped_multiple_mean_inp.columns =  grouped_multiple_mean_inp.columns.droplevel(-1)
grouped_multiple_mean_inp=grouped_multiple_mean_inp.assign(readidx='1')
grouped_multiple_mean_inp = grouped_multiple_mean_inp.reset_index()
grouped_multiple_mean_inp = grouped_multiple_mean_inp.iloc[:,[0,1,2,4,3]] 
grouped_multiple_mean_inp.to_csv('{}_processed_perpos_mean.csv'.format(infile), sep='\t', index = False) 