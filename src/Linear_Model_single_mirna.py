import numpy as np
import os
import pandas as pd
from utils_puma import get_processed_copdgene,get_stats_all_edges
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
from config import Config




mirna = snakemake.wildcards.mirna
path_output = snakemake.output.path_output
lm_exp = snakemake.wildcards.lm_exp
exp = snakemake.params.exp
config = Config(exp)
config.update_lm(lm_exp)

dict_input = get_processed_copdgene(config.path_dict_processed_input)
df_exp, df_pheno = dict_input['df_exp'], dict_input['df_pheno']

if config.lm_randomize:
    #shuffle the labels according to a random order defined in the preprocessing
    df_pheno = df_pheno.rename(index=dict_input['dict_rnd_pheno_index'])

del dict_input


df_pheno = df_pheno.loc[df_exp.columns]
df_pheno = df_pheno[[config.target_variable] + config.lm_confounder_columns]
df_pheno = df_pheno.dropna() 

print(path_output)
print(config.lm_confounder_columns)
print()


df_puma_by_mirna = pd.read_csv(config.path_puma_mirna + mirna+'.csv',index_col=0)
list_genes = df_puma_by_mirna.columns
list_pz = df_puma_by_mirna.index.intersection(df_pheno.index)
y = df_puma_by_mirna.loc[list_pz].values
y = StandardScaler().fit_transform(y)

if len(y)==0:
    pd.DataFrame().to_csv(path_output)
else:
    X = df_pheno.loc[list_pz,[config.target_variable] + config.lm_confounder_columns].copy()
    X = StandardScaler().fit_transform(X)
    lm_stats=get_stats_all_edges(X,y,range(1))
    lm_stats = np.swapaxes(lm_stats, 0, 2)
    m,n,r = lm_stats.shape
    out_arr = np.column_stack((np.repeat(np.arange(m),n),lm_stats.reshape(m*n,-1)))
    df_stats = pd.DataFrame(out_arr,index = list(list_genes)*m,columns = ['Variable','Coefficients','Standard_Errors','t','pvalues'])
    df_stats.Variable = df_stats.Variable.replace(dict(enumerate([config.target_variable])))
    
    df_stats['adj_pval'] = multipletests(df_stats.pvalues,method='fdr_bh')[1]
    print('Done statistics')
    os.makedirs('/'.join(path_output.split('/')[:-1]),exist_ok=True)
    df_stats.to_csv(path_output)




