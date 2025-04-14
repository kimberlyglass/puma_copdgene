import pandas as pd
from utils_puma import *
from scipy.spatial.distance import cdist
from config import Config
import sys
import shutil


def run_preprocess(config):
    #save the prints of the numbers for future usage
    path_log_preproc = config.path_output + 'Preprocessing_log.txt'
    sys.stderr = sys.stdout = open(path_log_preproc, "w")


    print('Creating copdgene expression files')

    df_exp_mirna, df_pheno_mirna = preprocess_mirna_expression(config.path_mirna_pheno,config.path_mirna_exp)
    df_exp,df_pheno = preprocess_mrna_expression(config.path_mrna_pheno,config.path_mrna_master,config.path_exp,config.deco_bcell, config.path_file_annot)
    df_motif,set_mirna,dict_mirna_merging = create_motif_prior_fam(config.path_targetscan,set_genes=set(df_exp.index),set_mirna = set(df_exp_mirna.index),
                                                                threshold_score=config.motif_threshold_score,path_data = config.path_data,path_output = config.path_output)


    # save only intersection of gene expression and mirna expression
    df_motif = df_motif.query('target in @df_exp.index and source in @df_exp_mirna.index')
    list_common_genes = df_exp.index.intersection(df_motif.target)
    list_common_mirna = df_exp_mirna.index.intersection(df_motif.source)
    df_exp = df_exp.loc[list_common_genes]
    df_motif = df_motif.query('target in @list_common_genes and source in @list_common_mirna')
    # keep all the mirna because the motif has merged mirnas, but i want to analyze all of them
    # df_exp_mirna = df_exp_mirna.loc[list_common_mirna]
    set_mirna = set(df_motif.source)

    print('Intersection genes targeted in motif')
    print('Gene expression:',df_exp.shape)
    print('Df motif:',df_motif.nunique())


    df_exp.to_csv(config.path_exp_cleaned)
    df_exp_mirna.to_csv(config.path_exp_mirna_cleaned)
    df_motif.to_csv(config.path_motif,header=False, sep='\t',index=False)
    with open(config.path_mirna, "w") as f:
        f.write('\n'.join(list(set_mirna)))


    df_corr = pd.DataFrame(np.corrcoef(df_exp), index = df_exp.index, columns = df_exp.index)
    df_corr_mirna = df_exp_mirna.T.corr()

    # correlation between the two genes in COPD
    list_mirna_mrna_sbj = df_exp.columns.intersection(df_exp_mirna.columns)
    print('\nSUbjects with mRNA and miRNA',len(list_mirna_mrna_sbj))
    df_corr_mirna_mrna = pd.DataFrame(1 - cdist(df_exp[list_mirna_mrna_sbj],df_exp_mirna[list_mirna_mrna_sbj], metric='correlation'),index = df_exp.index,columns=df_exp_mirna.index).T

    list_control = df_pheno.query('finalGold_P2 == 0').index
    list_case = df_pheno.query('finalGold_P2 >=2').index

    list_control_mirna = df_pheno_mirna.query('finalGold_P2 == 0').index
    list_case_mirna = df_pheno_mirna.query('finalGold_P2 >=2').index


    dict_input = {
        'list_case' : list_case,
        'list_control' : list_control,
        'list_case_mirna' : list_case_mirna,
        'list_control_mirna' : list_control_mirna,
        'df_exp' : df_exp,
        'df_pheno': df_pheno,
        'df_pheno_mirna':df_pheno_mirna,
        'df_corr':df_corr,
        'df_corr_mirna':df_corr_mirna,
        'df_corr_mirna_mrna':df_corr_mirna_mrna,
        'df_exp_mirna':df_exp_mirna,
        'df_motif':df_motif,
        'set_mirna':set_mirna,
        'dict_mirna_merging':dict_mirna_merging

    }

    # create ranom shuffle of the phenotype file for later use
    np.random.seed(123)
    dict_input['dict_rnd_pheno_index'] = dict(list(zip(*[df_pheno.index,np.random.permutation(df_pheno.index.values)])))


    os.makedirs(os.path.join(config.path_lm_coeff,'FEV1_FVC_post_P2'),exist_ok=True)

    with open(config.path_dict_processed_input,'wb') as f:
        pickle.dump(dict_input,f)

    #save lists of mirna sbj and genes for snakemake workflow
    list_all_pz = list(df_exp.columns)
    list_mirna = list(set_mirna)
    list_genes = list(df_exp.index)
    pd.DataFrame(list_all_pz).to_csv(config.path_sbj_list, index=False, header=False)
    pd.DataFrame(list_mirna).to_csv(config.path_mirna_list, index=False, header=False)
    pd.DataFrame(list_genes).to_csv(config.path_gene_list, index=False, header=False)

    #save config used
    shutil.copyfile('./config.py', config.path_output+'config.py')


