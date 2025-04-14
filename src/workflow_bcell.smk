import pandas as pd
from config import Config
from utils_puma import get_lm_coef_file_name, get_processed_copdgene

exp = 'bcell'
config = Config(exp)

# cc: cell composition, py:pack years
list_lm_exp = ['base','cc','py','cc_py','random'] 

list_sbj = pd.read_csv(config.path_sbj_list, header=None)[0].values
list_mirna = pd.read_csv(config.path_mirna_list, header=None)[0].values
list_genes = pd.read_csv(config.path_gene_list, header=None)[0].values


# Force snakemake not to run puma again
#import os 
#run = os.listdir('/d/tmp/remge/PUMA_final/Lioness_Puma/by_subj/')
#run = [x.split('.')[0] for x in run]
#list_sbj = [x for x in list_sbj if x not in run]
#print(list_sbj)

rule all: 
    input: 
        #expand(config.path_puma_sbj+'{sbj}.csv',sbj = list_sbj),
        #expand(config.path_puma_mirna+'{mirna}.csv',mirna=list_mirna),
        expand(config.get_lm_coef_file_name('{lm_exp}',config.target_variable,'{mirna}'), mirna=list_mirna,lm_exp=list_lm_exp),
        expand(config.get_lm_bcr_gene('{gene}'),gene = list_genes)


rule run_lio_puma:
  input:
    path_input_pickle = config.path_dict_processed_input,

  output:
    path_output = config.path_puma_sbj+'{sbj}.csv',

  conda:
      'pypandaenv1'
  script :
      'Lio_puma_single_id.py'



rule run_save_by_mirna:
  input:
    #path_input_pickle = config.path_dict_processed_input,
    path_input_lio_file = expand(config.path_puma_sbj+'{sbj}.csv',sbj = list_sbj),

  output:
    config.path_puma_mirna+'{mirna}.csv'
  params:
    path_puma_sbj = config.path_puma_sbj,
    path_puma_mirna = config.path_puma_mirna,
    list_mirna = list_mirna 
  resources:
    mem_mb=8000,  # Memory in MB
    time="01:30:00"  # Time in HH:MM:SS format
  script:
    'Create_files_by_mirna.py'
    


rule run_linear_model:
  input:
    path_input_pickle = config.path_dict_processed_input,
    path_mirna_input = config.path_puma_mirna+'{mirna}.csv'
  output:
    path_output = config.get_lm_coef_file_name('{lm_exp}',config.target_variable,'{mirna}')
  params:
    exp = exp 
  conda:
      'pypandaenv1'
  script :
      'Linear_Model_single_mirna.py'


rule run_bcr_lm:
  input:
    path_dict_processed_input = config.path_output+'dict_input_COPDGene.pickle',
    path_data_bcr_file = config.path_data_bcr_file,
    path_pickel_input = config.path_dict_processed_input,
  output:
      path_output = config.get_lm_bcr_gene('{gene}'),
  params:
    exp = exp
  conda:
      'pypandaenv1'
  script :
      'BCR_mRNA_groups_lm_single_gene.py'

