import pickle
from config import Config
from utils_puma import get_processed_copdgene,run_deseq_mrna_mirna


def run_deseq(config):
    

    df_deseq,df_deseq_cc,df_deseq_mirna,df_deseq_mirna_cc = run_deseq_mrna_mirna(config,save_file=True,filter_genes=True)

    dict_input = get_processed_copdgene(config.path_dict_processed_input)
    dict_input['df_deseq']=df_deseq
    dict_input['df_deseq_mirna']=df_deseq_mirna
    dict_input['df_deseq_cc']=df_deseq_cc
    dict_input['df_deseq_mirna_cc']=df_deseq_mirna_cc

    with open(config.path_dict_processed_input,'wb') as f:
        pickle.dump(dict_input,f)