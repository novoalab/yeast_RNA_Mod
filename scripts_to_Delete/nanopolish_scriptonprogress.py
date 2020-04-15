# New python script


# Non-executable
import pandas as pd
bc1=pd.read_csv("bc1.reads-ref.eventalign.txt",sep='\t')
#For the simple figure of (per median/position) 
grouped_multiple_mean_bc1 = bc1.groupby(['contig', 'position','reference_kmer']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_bc1 = grouped_multiple_mean_bc1.reset_index()
grouped_multiple_mean_bc1.columns =  grouped_multiple_mean_bc1.columns.droplevel(-1)
grouped_multiple_mean_bc1=grouped_multiple_mean_bc1.assign(readidx='1')
grouped_multiple_mean_bc1 = grouped_multiple_mean_bc1.reset_index()
grouped_multiple_mean_bc1 = grouped_multiple_mean_bc1.iloc[:,[0,1,2,4,3]] # Remember that Python does not slice inclusive of the ending index.
grouped_multiple_mean_bc1.to_csv(r'bc1_processed_perpos_mean.tsv', sep='\t', index = False)


# executable

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



# Non-executable
grouped_multiple_mean_bc1 = bc1.groupby(['contig', 'position','reference_kmer','read_index']).agg({'event_level_mean':['mean']})
grouped_multiple_mean_bc1 = grouped_multiple_mean_bc1.reset_index()
grouped_multiple_mean_bc1.columns =  grouped_multiple_mean_bc1.columns.droplevel(-1)
grouped_multiple_mean_bc1.to_csv(r'bc1_processed_perread_mean.tsv', sep='\t', index = False)









module load Python/3.7.2-GCCcore-8.2.0
python /users/enovoa/hliu/SHARE/huanle_useful_snippets/src/oz_sliding_win/5.Sliding_event_aln_tabl_on_ref_with_reads_info.py -f bc1_processed_perpos_modified.tsv -w 15 > bc1_perpos_sliding15.tsv
python /users/enovoa/hliu/SHARE/huanle_useful_snippets/src/oz_sliding_win/5.Sliding_event_aln_tabl_on_ref_with_reads_info.py -f bc2_processed_perpos_modified.tsv -w 15 > bc2_perpos_sliding15.tsv
python /users/enovoa/hliu/SHARE/huanle_useful_snippets/src/oz_sliding_win/5.Sliding_event_aln_tabl_on_ref_with_reads_info.py -f bc3_processed_perpos_modified.tsv -w 15 > bc3_perpos_sliding15.tsv
python /users/enovoa/hliu/SHARE/huanle_useful_snippets/src/oz_sliding_win/5.Sliding_event_aln_tabl_on_ref_with_reads_info.py -f bc4_processed_perpos_modified.tsv -w 15 > bc4_perpos_sliding15.tsv