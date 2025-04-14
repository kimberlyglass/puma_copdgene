import pandas as pd
from config import Config
from utils_puma import get_lm_coef_file_name, get_processed_copdgene

exp = 'bcell'
config = Config(exp)
dict_input = get_processed_copdgene(config.path_dict_processed_input)
list_genes = list(dict_input['df_exp'].index)
del dict_input



rule all: 
    input: 
        expand(config.get_lm_bcr_gene('{gene}'),gene = list_genes)


# all pz call -> ray parallelization
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

