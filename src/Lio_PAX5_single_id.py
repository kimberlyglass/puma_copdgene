from utils_puma import get_processed_copdgene,get_dict_ens_gene,to_iterator
import pandas as pd

from tqdm import tqdm
import numpy as np
import sys,os,json
from config import Config
import ray, psutil,pickle
from scipy.spatial.distance import cdist

import gseapy as gp


config = Config('bcell')
path_file_output = os.path.join(config.path_output,'PAX5_lioness_mRNA.csv') #snakemake.output.path_output
#path_input_pickle = snakemake.input.path_input_pickle
#subject_id = snakemake.wildcards.sbj
dict_input = get_processed_copdgene(config.path_dict_processed_input)

df_exp,df_corr,df_motif,set_mirna,df_pheno = dict_input['df_exp'],dict_input['df_corr'],dict_input['df_motif'],dict_input['set_mirna'], dict_input['df_pheno']
del dict_input
dict_ens_gene = get_dict_ens_gene()




### Ray Initialisation
num_cpus = psutil.cpu_count(logical=False) - 1
print('Using ',num_cpus,'cpus')
ray.init(log_to_driver=False, num_cpus=num_cpus)

dict_gene_ens = {y:x for x,y in dict_ens_gene.items()} 
ens_pax5 = dict_gene_ens['PAX5']

df_corr = pd.DataFrame(df_exp.drop(ens_pax5).T.corrwith(df_exp.loc[ens_pax5]))
list_results = []

@ray.remote
def run_lio_subj(subject_id):

    # Compute Lioness Correlation
    df_exp_subset = df_exp.drop(subject_id,axis=1)
    df_corr_subset = pd.DataFrame(df_exp_subset.T.corrwith(df_exp_subset.loc[ens_pax5]))

    n_conditions = df_exp.shape[1]
    df_lioness = n_conditions * (df_corr - df_corr_subset) + df_corr_subset
    df_lioness.columns = [subject_id]
    return df_lioness

    

res_id = [run_lio_subj.remote(subject_id) for subject_id in df_exp.columns]
list_results = [x for x in  tqdm(to_iterator(res_id), total=len(res_id))]

df_pax5_lio = pd.concat(list_results,axis=1)
df_pax5_lio.to_csv(path_file_output)


# run GSEA
#df_pax5_lio = pd.read_csv(path_file_output,index_col=0)
df_pax5_lio = df_pax5_lio.rename(index=dict_ens_gene)
tmp_common_sub = df_pax5_lio.columns.intersection(df_pheno.index)
#df_corr_PAX5_lio_pheno_bcell = pd.DataFrame(1 - cdist(df_pax5_lio[tmp_common_sub],df_pheno.loc[tmp_common_sub,config.target_variable].values.reshape(1,-1), metric='correlation'),index = df_pax5_lio.index).T

rank = pd.Series(1 - cdist(df_pax5_lio[tmp_common_sub],df_pheno.loc[tmp_common_sub,config.target_variable].values.reshape(1,-1), metric='correlation').reshape(-1),index = df_pax5_lio.index)
gene_sets=['GO_Biological_Process_2023']#'KEGG_2021_Human','GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023']

res_gsea = gp.prerank(rnk=rank,
                        gene_sets=gene_sets,threads=10,
                        min_size=10,max_size=50,permutation_num=10000, # reduce number to speed up testing
                        outdir=None, seed=6,verbose=True,no_plot=False)
res_gsea.res2d.Term = res_gsea.res2d.Term.str.replace('GO_Molecular_Function_2023__|GO_Cellular_Component_2023__|KEGG_2021_Human__|GO_Biological_Process_2023__','').str.split(' \(GO').str[0]

res_gsea.res2d.to_csv(config.path_output + 'GSEA_lio_PAX5.csv', index=None)