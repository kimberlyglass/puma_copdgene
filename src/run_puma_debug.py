import pandas as pd
import netZooPy.puma.puma as ntz

from utils_puma import get_processed_copdgene,run_puma_from_exp
import numpy as np


dict_old = get_processed_copdgene('/d/tmp/remge/Full_COPDGene_TS8_mirnafam_ncons_cont_new_mirna_2fam_merg/dict_input_COPDGene.pickle')


subject_id = '15340D'

np_corr_subset = np.corrcoef(dict_old['df_exp'].drop(subject_id, axis=1))
n_conditions = dict_old['df_exp'].shape[1]
np_lioness = n_conditions * (dict_old['df_corr'] - np_corr_subset) + np_corr_subset
df_corr = pd.DataFrame(np_lioness,index=dict_old['df_exp'].index,columns=dict_old['df_exp'].index)

# the order of the genes is given by the expression file on which we computed the correlation before
df_puma = run_puma_from_exp(None, dict_old['df_motif'], list(dict_old['set_mirna']),df_corr=df_corr)


df_puma.to_csv('/udd/remge/Projects/PUMA_final/tmp_Puma_debug.csv')