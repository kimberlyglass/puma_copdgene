import pandas as pd
import os 
import sys,os
sys.path.append('./')
from utils_puma import *
import pyreadr
from utils_puma import *
from config import Config


gene = snakemake.wildcards.gene
exp = snakemake.params.exp
path_output = snakemake.output.path_output
path_data_bcr_file = snakemake.input.path_data_bcr_file
path_pickel_input = snakemake.input.path_pickel_input
'''
gene='ENSG00000160957'
exp = 'bcell'
path_pickel_input = '/proj/regeps/regep00/studies/COPDGene/analyses/remge/PUMA_COPDGene/output_bcell/dict_input_COPDGene.pickle'
path_data_bcr_file = '/proj/regeps/regep00/studies/COPDGene/analyses/remge/PUMA_COPDGene/data/BCR/pheno_airr2.Rdata'
path_output = '/d/tmp/remge/bcr.csv'
'''

config = Config(exp)



dict_input = get_processed_copdgene(path_pickel_input)
df_exp_bcell, df_pheno_bcell = dict_input['df_exp'], dict_input['df_pheno']
df_exp_bcell = df_exp_bcell[df_pheno_bcell.index]


df_bcr = pyreadr.read_r(path_data_bcr_file) # also works for Rds
df_bcr = df_bcr['pheno'].set_index('sid')

list_common_sbj = df_bcr.index.intersection(df_exp_bcell.columns)
df_exp_bcell=df_exp_bcell[list_common_sbj]
df_bcr = df_bcr.loc[list_common_sbj]
df_bcr.columns = df_bcr.columns.str.replace('.','_')
df_bcr = df_bcr.dropna(subset=['smoking_status'])
df_bcr.smoking_status = df_bcr.smoking_status.astype(int)
df_bcr['prop_classSwitch'] = (df_bcr['prop_unmut_classSwitch']+df_bcr['prop_mut_classSwitch'])/2
# add gene expression
df_bcr[gene] = df_exp_bcell.loc[gene]


target_columns = ['prop_classSwitch','propCounts_IGHA1','propCounts_IGHA2','propCounts_IGHD','propCounts_IGHG1','propCounts_IGHG2','propCounts_IGHG3','propCounts_IGHM',
               'totMut_', 'totMut_IGHA1','totCount','totMut_IGHA2',
                'totMut_IGHD','totMut_IGHG1','totMut_IGHG2','totMut_IGHG3','totMut_IGHM',]


formula_covariates = 'standardize(%s)+standardize(FEV1_FVC_utah) +standardize(age_visit) + sex + race2 + smoking_status + standardize(BMI)'%gene


df_stat = get_stats_all_edges_bcr(df_bcr,formula_covariates=formula_covariates,y_columns= target_columns)
df_stat.to_csv(path_output)

'''try:
    df_stat = get_stats_all_edges_bcr(df_bcr,formula_covariates=formula_covariates,y_columns= target_columns)
    df_stat.to_csv(path_output)
except:
    # if fail save empty file
    pd.DataFrame().to_csv(path_output)'''