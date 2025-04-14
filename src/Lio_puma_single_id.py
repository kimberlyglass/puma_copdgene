from utils_puma import get_processed_copdgene,run_lio_puma
import pandas as pd

from tqdm import tqdm
import numpy as np
import sys,os,json
sys.path.insert(1, '/udd/remge/netZooPy/')
import netZooPy.puma.puma as ntz


path_file_output = snakemake.output.path_output
path_input_pickle = snakemake.input.path_input_pickle
subject_id = snakemake.wildcards.sbj

print(subject_id)
print('Computing: ',subject_id)
dict_input = get_processed_copdgene(path_input_pickle)

df_exp,df_corr,df_motif,set_mirna = dict_input['df_exp'],dict_input['df_corr'],dict_input['df_motif'],dict_input['set_mirna']
del dict_input
df_corr = df_corr.loc[df_exp.index,df_exp.index]
list_mirna = list(set_mirna)

run_lio_puma(subject_id, df_exp,df_corr.values,df_motif, list_mirna,path_file_output)


