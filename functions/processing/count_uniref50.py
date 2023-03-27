""" 
Script which filters mmseqs output to find the best annotation for each read 
"""

import pandas as pd
import glob
from collections import Counter 
import re 
#samples = ["/home/grig0076/CF_hackathon/mini_tab.csv", 
#           "/home/grig0076/CF_hackathon/mini_tab_2.csv"] 

#samples = ['/home/lidd0026/cfhackathon/788707_20171213_S.functions.gz',
#          '/home/lidd0026/cfhackathon/788707_20180301_S.functions.gz', 
#          '/home/lidd0026/cfhackathon/788707_20180313_S.functions.gz', 
#          '/home/lidd0026/cfhackathon/788707_20181126_S.functions.gz',] 

colnames = ["Read","Uniref_fxn","SeqID","Aln_length","Mismatches","Openings","q_start","q_end","t_start","t_end","Eval","Bit_score"]

samples = [i for i in glob.glob('/home/nala0006/scratch/atavide-dev/atavide/atavide.out/Illumina_read_based_annotations/*') if 'function' in i][70:]

#ignore the index files
samples = [s for s in samples if 'index' not in s] 

sample_counts = pd.DataFrame() 
counter = 0 

for s in samples: 
    print('processing: ' + s, flush = True)
    samp1 = pd.read_csv(s, sep="\t", compression='gzip',header = None)
    samp1.columns = colnames
    
    unique_reads_samp1 = list(set(samp1["Read"]))
    samp1f = samp1[samp1["Eval"]<1e-5] #filter by evalue 
    
    best_evalue = samp1f.groupby("Read", as_index=False).Eval.min().merge(samp1f)[['Read', 'Uniref_fxn']]
    
    #counts = Counter(tt['Uniref_fxn'])
    counts = Counter(best_evalue['Uniref_fxn'])

    #gives us the counts for this specific item 
    this_count = pd.DataFrame.from_dict(counts, orient = 'index')
    this_count.columns = [s]

    #save the current dataframe becuase of all of the memory issues 
    this_count.to_csv(re.split('/', s)[-1] + '.illuminaknowncountsonly.csv')
    
