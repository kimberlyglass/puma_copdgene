import numpy as np
import sys,os
import pandas as pd
from tqdm import tqdm 
from utils_puma import get_processed_copdgene,to_iterator,get_lm_coef_file_name,get_stats_all_edges
import psutil,ray

from sklearn.preprocessing import StandardScaler
from scipy import stats
from statsmodels.stats.multitest import multipletests
#import statsmodels.api as sm
from tqdm import tqdm


def linear_model(config):

    num_cpus = psutil.cpu_count(logical=False) - 1
    print('Using ',num_cpus,'cpus')
    ray.init(log_to_driver=False, num_cpus=num_cpus)
    

    ### Linear Model

    @ray.remote
    def run_linear_model_ray(mirna,df_pheno,list_conf_col,lm_regressions,path_puma_mirna,path_dir_results,skip_done):
        
        
        if skip_done and all([os.path.exists(get_lm_coef_file_name(path_dir_results,list_reg_var,mirna)) for list_reg_var in lm_regressions]):
            print('Done',mirna)
            return
        
        df_puma_by_mirna = pd.read_csv(path_puma_mirna + mirna+'.csv',index_col=0)
        list_genes = df_puma_by_mirna.columns
        list_pz = df_puma_by_mirna.index.intersection(df_pheno.index)
        y = df_puma_by_mirna.loc[list_pz].values
        try:
            y = StandardScaler().fit_transform(y)
        except:
            print(mirna)
            return

        for list_reg_var in lm_regressions:
            print(list_reg_var)
            try:
                file_output = get_lm_coef_file_name(path_dir_results,list_reg_var,mirna)
                if len(y)==0:
                    pd.DataFrame().to_csv(file_output)
                    continue
                X = df_pheno.loc[list_pz,list_reg_var + list_conf_col].copy()
                X = StandardScaler().fit_transform(X)
                lm_stats=get_stats_all_edges(X,y,range(len(list_reg_var)))
                lm_stats = np.swapaxes(lm_stats, 0, 2)
                m,n,r = lm_stats.shape
                out_arr = np.column_stack((np.repeat(np.arange(m),n),lm_stats.reshape(m*n,-1)))
                df_stats = pd.DataFrame(out_arr,index = list(list_genes)*m,columns = ['Variable','Coefficients','Standard_Errors','t','pvalues'])
                df_stats.Variable = df_stats.Variable.replace(dict(enumerate(list_reg_var)))
                df_stats['adj_pval'] = multipletests(df_stats.pvalues,method='fdr_bh')[1]
                print('DOne statistics')
                os.makedirs('/'.join(file_output.split('/')[:-1]),exist_ok=True)
                df_stats.to_csv(file_output)
            except Exception as e: 
                print(e)
                #pd.DataFrame().to_csv(file_output)
                continue
            
        return 



    def run_linear_model(df_pheno,path_puma_mirna,set_mirna,list_conf_col,randomize,path_output,skip_done=False):
        all_target_columns = list(set(sum(lm_regressions, []))) # all the variables used
        tmp_df_pheno = df_pheno[list(set(all_target_columns + list_conf_col))]
        tmp_df_pheno = tmp_df_pheno.dropna()

        if randomize:
            np.random.shuffle(tmp_df_pheno.index.values)  
            tmp_df_pheno.to_csv(path_output+'Pheno_random.csv')     
        id_tmp_df_pheno = ray.put(tmp_df_pheno)
        print(path_output)
        print(lm_regressions)
        print(list_conf_col)
        print()
        
        res_id = [run_linear_model_ray.remote(mirna,id_tmp_df_pheno,list_conf_col,lm_regressions,path_puma_mirna,path_output,skip_done) for mirna in set_mirna]

        
        for res_edge in tqdm(to_iterator(res_id), total=len(res_id)):
            pass

            
    path_lm_coeff = config.path_lm_coeff
    path_puma_mirna = config.path_puma_mirna 
    path_lm_coeff_randomize = config.path_lm_coeff_randomize
    
    list_conf_col = config.lm_confounder_columns
    lm_regressions = config.lm_regressions
    

    dict_input = get_processed_copdgene(config.path_dict_processed_input)
    df_exp, df_pheno,set_mirna = dict_input['df_exp'], dict_input['df_pheno'], dict_input['set_mirna']
    del dict_input

    list_all_pz = list(df_exp.columns)

    df_pheno = df_pheno.loc[list_all_pz]

    # Run Linear model
    print('Run Linear Model')
    run_linear_model(df_pheno,path_puma_mirna,set_mirna,list_conf_col,randomize=False,path_output=path_lm_coeff, skip_done=True)

    #list_conf_col_no_bmi = [x for x in list_conf_col if 'bmi' not in x.lower()]
    #if len(list_conf_col_no_bmi)!=list_conf_col:
    #    run_linear_model(df_pheno,path_puma_mirna,set_mirna,list_conf_col_no_bmi,randomize=False,path_output=path_lm_coeff[:-1]+'_nobmi/')
    run_linear_model(df_pheno,path_puma_mirna,set_mirna,list_conf_col + config.cell_count_cols,randomize=False,path_output=path_lm_coeff[:-1]+'_cellcount/')
    run_linear_model(df_pheno,path_puma_mirna,set_mirna,list_conf_col +["ATS_PackYears_P2"],randomize=False,path_output=path_lm_coeff[:-1]+'_packyears/')
    run_linear_model(df_pheno,path_puma_mirna,set_mirna,list_conf_col + config.cell_count_cols + +["ATS_PackYears_P2"],randomize=False,path_output=path_lm_coeff[:-1]+'_packyears_cellcount/')
    if path_lm_coeff_randomize is not None:
        run_linear_model(df_pheno,path_puma_mirna,set_mirna,list_conf_col,randomize=True,path_output=path_lm_coeff_randomize)
