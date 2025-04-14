from utils_puma import get_processed_copdgene,run_puma_from_exp,to_iterator
import pandas as pd
import ray,psutil,os
from tqdm import tqdm
import numpy as np
import sys,os,json,warnings
sys.path.insert(1, '/udd/remge/netZooPy/')
import netZooPy.puma.puma as ntz


def run_lioness_puma(config):

    path_puma_sbj = config.path_puma_sbj

    ### Ray Initialisation
    num_cpus = (psutil.cpu_count(logical=False) - 1)
    print('Using ',num_cpus,'cpus')
    ray.init(log_to_driver=False, num_cpus=num_cpus)



    @ray.remote
    def run_puma_ray(subject_id,df_exp,np_corr_full,df_motif,miR):
        if os.path.isfile(path_puma_sbj+subject_id+'_chandl02.csv'):
            return
        # Compute Lioness Correlation
        np_corr_subset = np.corrcoef(df_exp.drop(subject_id, axis=1))
        n_conditions = df_exp.shape[1]
        np_lioness = n_conditions * (np_corr_full - np_corr_subset) + np_corr_subset
        df_corr = pd.DataFrame(np_lioness,index=df_exp.index,columns=df_exp.index)#*np.sign(np_corr_full)

        # the order of the genes is given by the expression file on which we computed the correlation before
        df_puma = run_puma_from_exp(None, df_motif, miR,df_corr=df_corr)
        try:
            df_puma.to_csv(path_puma_sbj+subject_id+'_chandl02.csv')
        except Exception as e:
            warnings.warn(f"An error occurred: {e}")
            warnings.warn(path_puma_sbj+subject_id+'_chandl02.csv')


    dict_input = get_processed_copdgene(config.path_dict_processed_input)
    df_exp,df_corr,df_motif = dict_input['df_exp'],dict_input['df_corr'],dict_input['df_motif']
    del dict_input


    id_df_exp = ray.put(df_exp)
    np_corr_full = df_corr.loc[df_exp.index,df_exp.index].values
    id_np_corr_full = ray.put(np_corr_full)

    list_mirna = list(df_motif.source.unique())

    print('Motif and Genes: ',df_motif[['source', 'target']].nunique())
    print('miRNA: ',len(list_mirna))



    id_df_motif = ray.put(df_motif)
    id_miR = ray.put(list_mirna)

    print('df_motif',df_motif.head(), df_motif.shape)

    res_id = [run_puma_ray.remote(subject_id, id_df_exp,id_np_corr_full,id_df_motif, id_miR) for subject_id in ['15340D']]#df_exp.columns[::-1]]

    # Get the results.
    for x in tqdm(to_iterator(res_id), total=len(res_id)):
        pass
