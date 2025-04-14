import pandas as pd
import pickle,sys
import netZooPy.puma.puma as ntz
import ray,os, json, re
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objs as go
import plotly.express as px
from scipy.stats import ttest_ind,norm, gaussian_kde,ttest_1samp
from statsmodels.stats.multitest import multipletests
import scipy
from sklearn.preprocessing import StandardScaler
import pyreadr
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import pairwise_distances
from scipy.sparse import csr_matrix 
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pyranges as pr


import seaborn as sns
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression

from scipy import stats
from statsmodels.stats.multitest import multipletests
from patsy import dmatrices
from rnanorm import TMM




'''
def import_config(file_path_experiment):
    import importlib.util
    file_path = file_path_experiment+'config.py'
    # Module name (without extension) - you can customize this
    module_name = "config"
    # Create a spec from the file path
    print(file_path)
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    # Import the module
    config = importlib.util.module_from_spec(spec)
    # Load the module
    spec.loader.exec_module(config)
    return config
'''
def get_processed_copdgene(path_file_exp):

    with open(path_file_exp,'rb') as f:
        dict_input = pickle.load(f)
    #dict_input['df_exp'] = dict_input['df_exp'].loc[dict_input['df_motif'].target.unique()]
    #with open(path_file_exp,'wb') as f:
    #   pickle.dump(dict_input,f)
    return dict_input#dict_input['df_exp'], dict_input['df_pheno_master'], dict_input['df_mirna_pheno'],dict_input['df_corr'],dict_input['df_exp_mirna']



### Gene Expression

def preprocess_mrna_expression(path_mrna_pheno,path_mrna_master,path_exp,deco_bcell,path_file_annot):

    # new Phenotype file related to gene Expression ~ 10k patients
    df_exp_pheno = pd.read_csv(path_mrna_pheno, sep ='\t',na_values=' ',index_col=0,low_memory=False)
    print('Phenotype',df_exp_pheno.shape)

    # Masterfile 
    df_master_file = pd.read_csv(path_mrna_master,sep='\t')
    print('Masterfile:',df_master_file.shape)

    df_master_file = df_master_file[df_master_file['final.analyzed']==1].set_index('actual_id')
    print('final.analyzed==1', df_master_file.shape)
    df_pheno_master = df_master_file.merge(df_exp_pheno, left_index=True,right_index=True)
    df_pheno_master = df_pheno_master.groupby('Lc.Batch').filter(lambda x:len(x)>10)
    #drop finalgold -1 and missin FEV1_FVC_post_P2 and never smokers(==0)
    df_pheno_master = df_pheno_master.query('finalGold_P2>=0 and smoking_status_P2!=0 and FEV1_FVC_post_P2.notna()')
    df_pheno_master = df_pheno_master.dropna(subset=['FEV1_FVC_post_P2'])
    print('Drop batches <10pz and finalGold_P2>=0 and smoking_status_P2!=0',df_pheno_master.shape)

    df_pheno_master['super_control'] = np.all(df_pheno_master[['finalGold_P1','finalGold_P2','finalGold_P3']]==0,axis=1)
    #df_pheno_master["smoking_status_P2"]=-(df_pheno_master["smoking_status_P2"]-2) # 1 if you smoke
    

    df_exp = pd.read_csv(path_exp, sep ='\t',index_col=0)
    print('Gene expression: ',df_exp.shape)
    
    if not deco_bcell:
        df_exp.index = df_exp.index.str.split('.').str[0]
        # select only protein coding genes
        df_ann = get_df_gene_annotations(path_file_annot)
        df_ann = df_ann.loc[df_exp.index].query('gene_type=="protein_coding"')
        df_exp = df_exp.loc[df_ann.index]
        print('Gene expression only pretein coding genes: ',df_exp.shape)

    # common sampl expression and pheno
    list_common_sampl = df_pheno_master.index.intersection(df_exp.columns)
    df_pheno_master = df_pheno_master.loc[list_common_sampl]
    df_exp = df_exp[list_common_sampl]

    if deco_bcell:
        dict_ens_gene = get_dict_ens_gene()
        dict_gene_ens = {y:x for x,y in dict_ens_gene.items()} #34-zeb1
        df_exp = df_exp.rename(index = dict_gene_ens)
        #df_exp = df_exp[df_exp.quantile(.90,axis=1)>=10]
    print('Final shape intersection with phenotype:',df_exp.shape)
    return df_exp,df_pheno_master

def get_dict_ens_gene(path_to_ens_mapping = None ):
    if path_to_ens_mapping is None:
        from config import Config
        path_to_ens_mapping = Config('bulk').path_to_ens_mapping
    return  pd.read_csv(path_to_ens_mapping,sep='\t').set_index('Gene stable ID')['Gene name'].to_dict()

def map_ens_to_genesymb(df_exp):
    dict_ens_gene = get_dict_ens_gene()
    df_exp.index = list(map(lambda x: x.split('.')[0],df_exp.index))
    df_exp = df_exp.rename(index=dict_ens_gene)
    df_exp = df_exp[~df_exp.index.str.contains('ENS')]
    df_exp.index.name = 'Gene'
    print('Gene exp after mapping: ',df_exp.shape)
    
    # remove double genes
    df_exp = df_exp.loc[df_exp.index.value_counts()==1]
    print('Gene exp after duplicates removal: ',df_exp.shape)
    
    return df_exp

def get_df_gene_annotations(path_file):
    if not os.path.isfile(path_file.replace('.gtf','.csv')):
        df_ann = pr.read_gtf(path_file)
        df_ann = df_ann.df
        df_ann['gene_id_stable'] = df_ann.gene_id.str.split('.').str[0]
        df_ann = df_ann[['gene_id_stable','gene_type','Chromosome','Start','End']]

        df_ann = df_ann.drop_duplicates()
        #only gene in X and Y are duplicated. So I can use the gene_id_stable as index
        df_ann[df_ann.gene_id_stable.duplicated(keep=False)].Chromosome.value_counts()
        df_ann = df_ann.drop_duplicates(subset=['gene_id_stable'])
        df_ann = df_ann.set_index('gene_id_stable')
        df_ann.to_csv(path_file.replace('.gtf','.csv'))
    else:
        df_ann = pd.read_csv(path_file.replace('.gtf','.csv'),index_col=0)
    return df_ann



### miRNA Preprocessing

def preprocess_mirna_expression(path_mirna_pheno,path_mirna_exp):

    # Load Data

    #file_miRNA_exp =  '/proj/regeps/regep00/studies/COPDGene/analyses/rebdh/craigMiRna/mirSeqCounts_filt_blockFreeze2_20200608.rds'
    #file_miRNA_pheno =  '/proj/regeps/regep00/studies/COPDGene/analyses/rebdh/craigMiRna/mirSeqPheno_filt_blockFreeze2_20200608.rds'
    df_exp = pyreadr.read_r(path_mirna_exp)[None]
    df_pheno = pyreadr.read_r(path_mirna_pheno)[None]
    print('Input MiRNA expression: ',df_exp.shape)

    # clean index MiRNA
    df_exp.index = df_exp.index.map(mirna_str_proc)
    #remove x from subject ID
    df_exp.columns = df_exp.columns.str[1:]
    df_pheno = df_pheno.set_index('sid')

    # FEV1_FVC_utah_P2 is the same compared to the phenotype of mrna. Aligning the name
    df_pheno = df_pheno.rename(columns={'FEV1_FVC_utah_P2':'FEV1_FVC_post_P2'})  
    
    # Drop overexpressed miRNA in blood from Brian's papers
    # https://doi.org/10.1016/j.jmoldx.2021.03.006
    # https://doi.org/10.3389/fgene.2021.748356


    df_exp = df_exp.drop(['mir-191-5p','mir-486-5p','mir-92a-3p'])
    print('Drop 3 miRNA: ',df_exp.shape)

    # Phenotype 
    print('\nInput miRNA_pheno shape: ',df_pheno.shape)
    # Remove batch plate5 read Brian's papers for info
    df_pheno = df_pheno.query('mirSeqBatch!="plate5"')
    print('\tRemove plate5: ',df_pheno.shape)
    # Patients with missing clinical information
    print('Missing FEV1/FVC:',df_pheno[df_pheno.isna()['FEV1_FVC_post_P2']].index)

    df_pheno = df_pheno.dropna(subset=['FEV1_FVC_post_P2','finalGold_P2'])
    print('\tRemove missing FEV1_FVC_post_P2 or finalGold_P2: ',df_pheno.shape)

    # Intersect MiRNA exp with Phenotype 
    df_exp = df_exp[df_pheno.index]
    print('Intersect MiRNA exp with Phenotype: ',df_exp.shape)
    df_exp = df_exp[df_exp.columns[df_exp.sum(0)>200000]]
    print('\tRemove low sequencing depth <200k',df_exp.shape)
    # filtered out absent and low-variant miRNAs 
    # more than 10 reads in at least 200 subjects, 
    print('miRNA preprocessing')
    df_exp = df_exp[(df_exp>=10).sum(1)>=200]
    print('\tMore than 10 reads in at least 200:',df_exp.shape)
    # minimum standard deviation of 10 across subjects 
    df_exp = df_exp[df_exp.std(1)>=10]
    print('\tStd >=10:',df_exp.shape)

    # TMM and CPM
    tmm = TMM().set_output(transform="pandas")
    tmm.fit(df_exp.T)
    tmm.get_norm_factors(df_exp.T)
    df_exp_tmm = tmm.transform(df_exp.T).T

    return df_exp_tmm, df_pheno

def mirna_str_proc(mirna):
    mirna = mirna.lower().replace('hsa-','')#.split('.')[0].replace('mir-','mir').replace('let-','mirlet').upper()
    #mirna = re.sub(r'-\d-', '-', mirna)
    return mirna
'''
def old_process_mirna_exp():
    ### OLD. Now read file preprocess by Brian.
    # Load Data
    
    file_miRNA_exp =  '/proj/regeps/regep00/studies/COPDGene/analyses/rebdh/craigMiRna/mirSeqCounts_filt_blockFreeze2_20200608.rds'
    file_miRNA_pheno =  '/proj/regeps/regep00/studies/COPDGene/analyses/rebdh/craigMiRna/mirSeqPheno_filt_blockFreeze2_20200608.rds'
    df_miR_exp = pyreadr.read_r(file_miRNA_exp)[None]
    df_mirna_pheno_new = pyreadr.read_r(file_miRNA_pheno)[None]
    print('Input MiRNA expression: ',df_miR_exp.shape)

    # clean index MiRNA
    df_miR_exp.index = [x.replace('hsa-','') for x in df_miR_exp.index]
    df_miR_exp.columns = [x if x[0]!='X' else x[1:] for x in df_miR_exp.columns]

    # Drop overexpressed miRNA
    df_miR_exp = df_miR_exp.drop(['miR-191-5p','miR-486-5p','miR-92a-3p'])
    print('Drop 3 miRNA: ',df_miR_exp.shape)

    # Phenotype 
    print('\nInput miRNA_pheno shape: ',df_mirna_pheno_new.shape)
    # Remove batch plate5
    df_mirna_pheno_new = df_mirna_pheno_new.query('mirSeqBatch!="plate5"')
    print('\tRemove plate5: ',df_mirna_pheno_new.shape)
    # Patients with missing clinical information
    df_mirna_pheno_new = df_mirna_pheno_new.dropna(subset=['Perc15_Insp_Thirona_P2'])
    df_mirna_pheno_new = df_mirna_pheno_new.dropna(subset=['finalGold_P2'])
    print('\tRemove missing Perc15_Insp_Thirona_P2 or finalGold_P2: ',df_mirna_pheno_new.shape)

    # Intersect MiRNA exp with Phenotype 
    df_mirna_pheno_new = df_mirna_pheno_new.set_index('sid')
    df_miR_exp = df_miR_exp[df_mirna_pheno_new.index]
    print('\nIntersect MiRNA exp with Phenotype: ',df_miR_exp.shape)
    df_miR_exp = df_miR_exp[df_miR_exp.columns[df_miR_exp.sum(0)>200000]]
    print('Remove low sequencing depth <200k',df_miR_exp.shape)
    # filtered out absent and low-variant miRNAs 
    # more than 10 reads in at least 200 subjects, 
    df_miR_exp = df_miR_exp[(df_miR_exp>=10).sum(1)>=200]
    print('More than 10 reads in at least 200:',df_miR_exp.shape)
    # minimum standard deviation of 10 across subjects 
    df_miR_exp = df_miR_exp[df_miR_exp.std(1)>=10]
    print('Std >=10:',df_miR_exp.shape)

    
    df_miR_exp.to_csv('../Data/COPDGene/Processed/MicroRNA_clean.csv')


'''
### Motif Preprocessing 

def get_mirna_family(set_mirna,path_data):
    df_mirna_family = pd.read_csv(path_data + 'motif/miR_Family_Info.txt.zip',compression='zip',sep='\t')
    df_mirna_family = df_mirna_family[df_mirna_family['MiRBase ID'].str.contains('hsa')]
    df_mirna_family['MiRBase ID'] = df_mirna_family['MiRBase ID'].apply(mirna_str_proc)
    df_mirna_family.columns = ['_'.join(x.split(' ')) for x in df_mirna_family.columns]
    df_mirna_family = df_mirna_family.query('MiRBase_ID in @set_mirna')
    df_mirna_family = df_mirna_family.drop_duplicates(subset=['miR_family','MiRBase_ID'])
    df_mirna_family['weight']=1
    df_mirna_family = df_mirna_family.groupby('miR_family').filter(lambda x: len(x)>=2)
    df_mirna_family_adj = df_mirna_family.pivot(index='miR_family', columns='MiRBase_ID', values = 'weight')
    print('No mirna belong to more than one family, max n fam per mirna:',df_mirna_family_adj.sum(0).max())
    assert(df_mirna_family_adj.sum(0).max()==1)
    return df_mirna_family,df_mirna_family_adj
'''
def create_motif_prior_v6():
    df_motif_60 = pd.read_csv('../Data/Summary_Counts_6.txt',sep='\t')

    verbose=True
    score_col = 'Total context score'
    if verbose: print('Initial')
    if verbose:print(df_motif_60[['Gene Symbol','Representative miRNA']].nunique())
    # only Homo Sapiens interactions
    df_motif_60 = df_motif_60[df_motif_60['Representative miRNA'].str.contains('hsa')]
    if verbose: print('Homo Sapiens')
    if verbose: print(df_motif_60[['Gene Symbol','Representative miRNA']].nunique())
    # removing cotext++ score >-0.1
    df_motif_6 = df_motif_60.copy()
    df_motif_6 = df_motif_6[df_motif_6[score_col]<-0.1]
    if verbose: print('context++ score <=-.1')
    if verbose: print(df_motif_6[['Gene Symbol','Representative miRNA']].nunique())
    
    df_motif_6['Representative miRNA'] = df_motif_6['Representative miRNA'].apply(mirna_str_proc)
    df_motif_6 = df_motif_6.drop_duplicates(subset=['Gene Symbol','Representative miRNA'])
    if verbose: print('Dropping duplicates from different loci and small variances')
    if verbose: print(df_motif_6[['Gene Symbol','Representative miRNA']].nunique())


    df_motif_6 = df_motif_6.rename(columns={'Representative miRNA':'source','Gene Symbol':'target'})

    df_motif_6 = df_motif_6[['source','target']].drop_duplicates()
    df_motif_6['weight'] = 1

    miR_6 = list(df_motif_6.source.unique())
    return df_motif_6,miR_6

def create_motif_prior(path_targetscan, path_mapping_t_g,verbose=False):

    df_motif = pd.read_csv(path_targetscan,sep='\t',usecols=['Transcript ID','Representative miRNA','Cumulative weighted context++ score'])
    if verbose: print('Initial')
    if verbose:print(df_motif[['Transcript ID','Representative miRNA']].nunique())
    # only Homo Sapiens interactions
    #df_motif = df_motif[df_motif['Representative miRNA'].str.contains('hsa')]
    if verbose: print('Homo Sapiens')
    if verbose: print(df_motif[['Transcript ID','Representative miRNA']].nunique())
    # removing cotext++ score >-0.1
    df_motif = df_motif[df_motif['Cumulative weighted context++ score']<-0.1]
    if verbose: print('Cumulative weighted context++ score <-0.1')
    if verbose: print(df_motif[['Transcript ID','Representative miRNA']].nunique())

    #The numbered suffix in the Representative miRNA identifiers indicate
    #diverse loci that produce identical mature Representative miRNAs. We
    #therefore collapsed these Representative miRNAs in the prior by taking
    #the union of all edges.
    df_motif['Representative miRNA'] = df_motif['Representative miRNA'].apply(mirna_str_proc)
    df_motif = df_motif.drop_duplicates(subset=['Transcript ID','Representative miRNA'])

    ### Map Transcript to Gene and Save
    df_ens_mapping = pd.read_csv(path_mapping_t_g,sep='\t')
    df_ens_mapping['Gene stable ID'] = df_ens_mapping['Gene stable ID version'].str.split('.').str[0]
    dict_mapping = df_ens_mapping.set_index('Transcript stable ID version')['Gene stable ID'].to_dict()

    df_motif['Gene ID'] = df_motif['Transcript ID'].map(dict_mapping)

    # CDR1as transcript  not mapped
    df_motif = df_motif.dropna(subset=['Gene ID'])

    if verbose: print('Dropping duplicates after mapping to ensg')
    if verbose: print(df_motif[['Transcript ID','Representative miRNA']].nunique())

    df_motif = df_motif.rename(columns={'Representative miRNA':'source','Gene ID':'target'})

    df_motif = df_motif[['source','target']].drop_duplicates()
    df_motif['weight'] = 1

    miR = list(df_motif.source.unique())
    #df_motif.to_csv('../Data/Df_motif.csv',index=False)

    return df_motif,set(miR)

def create_motif_prior_fam_old(path_targetscan,set_genes,set_mirna,threshold_score,path_output):
    
    #from targetscan
    #df_motif_c = pd.read_csv(config.path_data+'motif/Conserved_Site_Context_Scores.txt.zip',compression='zip',sep='\t')
    #df_motif_c = df_motif_c[df_motif_c.miRNA.str.contains('hsa')]
    #df_motif_c.columns = ['_'.join(x.split(' ')) for x in df_motif_c.columns]
    #df_motif_c = df_motif_c.sort_values('weighted_context++_score').drop_duplicates(subset=['Gene_ID','miRNA']).copy()
    #df_motif_c.Gene_ID = df_motif_c.Gene_ID.str.split('.').str[0]
    #df_motif.to_csv(path_data + 'motif/Conserved_Site_Context_Scores_hsa.csv',index=False)

    #df_motif = pd.read_csv('/d/tmp/remge/Nonconserved_Site_Context_Scores.txt.zip',compression='zip',sep='\t')
    #df_motif = df_motif[df_motif.miRNA.str.contains('hsa')]
    #df_motif.columns = ['_'.join(x.split(' ')) for x in df_motif.columns]
    #df_motif = df_motif.sort_values('weighted_context++_score').drop_duplicates(subset=['Gene_ID','miRNA']).copy()
    #df_motif.Gene_ID = df_motif.Gene_ID.str.split('.').str[0]
    #df_motif.to_csv(path_targetscan,index=False)

    df_motif = pd.read_csv(path_targetscan)
    df_motif.Gene_ID = df_motif.Gene_ID.str.split('.').str[0]
    df_motif.miRNA = df_motif.miRNA.map(mirna_str_proc)
    df_motif = df_motif.query('Gene_ID in @set_genes and miRNA in @set_mirna')
    df_mirna_family,df_mirna_family_adj = get_mirna_family(set_mirna)


    df_motif_adj_binary = (df_motif.set_index(['miRNA','Gene_ID'])['weighted_context++_score'].unstack()<threshold_score).astype(int)

    df_mirna_motif_jaccard,df_mirna_motif_jaccard_list = get_jaccard(df_motif_adj_binary)

    ### N common and different edges (mirna-mrna) between mirnas
    tmp_sparse = csr_matrix(df_motif_adj_binary)
    df_count_co = pd.DataFrame((tmp_sparse*tmp_sparse.T).todense(), index = df_motif_adj_binary.index, columns=df_motif_adj_binary.index)
    df_count_diff = df_count_co*(1/df_mirna_motif_jaccard-1)
    df_count_diff_list = df_count_diff.unstack().dropna()
    df_count_diff_list.index.rename(['mirna1','mirna2'],inplace=True)
    df_count_diff_list = df_count_diff_list.reset_index().rename(columns={0:'n_diff_edges'})

    ## MERGE MIRNA if same family AND Jaccard>.75 AND n_diff_edges<50
    # get info on pairwise mirna: Jaccard similarity and number of different edges and same family
    df_merge_sel = df_mirna_motif_jaccard_list.query('Jaccard>.75').merge(df_count_diff_list.query('n_diff_edges<50'),how='inner')
    df_merge_sel = df_merge_sel.merge(df_mirna_family[['miR_family','MiRBase_ID']],left_on='mirna1',right_on='MiRBase_ID',how='inner').drop('MiRBase_ID',axis=1)
    df_merge_sel = df_merge_sel.merge(df_mirna_family[['miR_family','MiRBase_ID']],left_on='mirna2',right_on='MiRBase_ID',how='inner').drop('MiRBase_ID',axis=1)
    print('All similar mirna belong to same family:',df_merge_sel.query('miR_family_x!=miR_family_y').shape)
    assert(np.all(df_merge_sel.miR_family_x==df_merge_sel.miR_family_y))

    # merging mirnas in the same family
    dict_mirna_merging = list(map(lambda x: list(set(x[1][['mirna1','mirna2']].values.flatten())),df_merge_sel.groupby('miR_family_x')))
    #pick first mirna as mirna name used when merged
    dict_mirna_merging = {y:x[0] for x in dict_mirna_merging for y in x } 
    
    #merge mirnas
    df_motif_adj = df_motif.set_index(['miRNA','Gene_ID'])['weighted_context++_score'].unstack()
    df_motif_adj = df_motif_adj.rename(index=dict_mirna_merging)
    # if mirna are merged keep the lower score.
    df_motif_adj = df_motif_adj.groupby('miRNA').min()
    df_motif_adj = (df_motif_adj<threshold_score)
    # in case some mirna has no edge, drop it
    df_motif_adj = df_motif_adj[np.any(df_motif_adj,axis=1)]


    # Plots
    #plot_jaccard_vs_count(df_motif_adj_binary,df_mirna_family_adj,config.path_figures + 'Jaccard_motif.pdf')
    #plot_jaccard_vs_count(df_motif_adj.astype(int),df_mirna_family_adj_new,config.path_figures + 'Jaccard_motif_clean.pdf')

    df_motif_final = df_motif_adj.unstack()
    df_motif_final = df_motif_final[df_motif_final].astype(int).reset_index().rename(columns={'Gene_ID':'target','miRNA':'source',0:'weight'})
    df_motif_final = df_motif_final[['source','target','weight']]

    # rename  the mirna family with the name of one of the mirnas
    pd.DataFrame(dict_mirna_merging.items(),columns=['mirna','merged_mirna']).to_csv(path_output+'miRNA_merging_mapping.csv',index=False)

    return df_motif_final,set(df_motif_adj.index),dict_mirna_merging

'''
def create_motif_prior_fam(path_targetscan,set_genes,set_mirna,threshold_score,path_data,path_output):
    
    #from targetscan
    #df_motif_c = pd.read_csv(config.path_data+'motif/Conserved_Site_Context_Scores.txt.zip',compression='zip',sep='\t')
    #df_motif_c = df_motif_c[df_motif_c.miRNA.str.contains('hsa')]
    #df_motif_c.columns = ['_'.join(x.split(' ')) for x in df_motif_c.columns]
    #df_motif_c = df_motif_c.sort_values('weighted_context++_score').drop_duplicates(subset=['Gene_ID','miRNA']).copy()
    #df_motif_c.Gene_ID = df_motif_c.Gene_ID.str.split('.').str[0]
    #df_motif.to_csv(path_data + 'motif/Conserved_Site_Context_Scores_hsa.csv',index=False)

    #df_motif = pd.read_csv('/d/tmp/remge/Nonconserved_Site_Context_Scores.txt.zip',compression='zip',sep='\t')
    #df_motif = df_motif[df_motif.miRNA.str.contains('hsa')]
    #df_motif.columns = ['_'.join(x.split(' ')) for x in df_motif.columns]
    # drop duplicated because the miRNA sequence can overlap to multiple places in the same trascript. We keep the overlap with best (lowest) score
    #df_motif = df_motif.sort_values('weighted_context++_score').drop_duplicates(subset=['Gene_ID','miRNA']).copy()
    #df_motif.Gene_ID = df_motif.Gene_ID.str.split('.').str[0]
    #df_motif.to_csv(path_targetscan,index=False)

    print('\nMotif prior preprocessing')
    df_motif = pd.read_csv(path_targetscan)
    #df_motif.Gene_ID = df_motif.Gene_ID.str.split('.').str[0]
    df_motif.miRNA = df_motif.miRNA.map(mirna_str_proc)
    print('hsa motif\n', df_motif[['Gene_ID','miRNA']].nunique(), len(df_motif))
    df_motif = df_motif.query('Gene_ID in @set_genes and miRNA in @set_mirna')
    df_mirna_family,df_mirna_family_adj = get_mirna_family(set_mirna,path_data)

    print('With expression \n', df_motif[['Gene_ID','miRNA']].nunique(), len(df_motif))

    # merging mirnas in the same family
    dict_mirna_merging = df_mirna_family.set_index('MiRBase_ID').miR_family.to_dict()
    tmp_dict = {y:x for x,y in dict_mirna_merging.items()}
    dict_mirna_merging = {x:tmp_dict[y] for x,y in dict_mirna_merging.items()}#pick first mirna as mirna name used when merged

    #merge mirnas
    df_motif_adj = df_motif.set_index(['miRNA','Gene_ID'])['weighted_context++_score'].unstack()
    df_motif_adj = df_motif_adj.rename(index=dict_mirna_merging)
    # if mirna are merged keep the lower score.
    df_motif_adj = df_motif_adj.groupby('miRNA').min()
    print('merging families motif\n', df_motif_adj.shape)

    if threshold_score =='cont':
        df_motif_final = df_motif_adj.unstack().reset_index().rename(columns={'Gene_ID':'target','miRNA':'source',0:'weight'})
        df_motif_final = df_motif_final.dropna(subset=['weight'])
        df_motif_final['weight'] = -df_motif_final['weight']
    else:
        df_motif_adj = (df_motif_adj<threshold_score)
        # in case some mirna has no edge, drop it
        df_motif_adj = df_motif_adj[np.any(df_motif_adj,axis=1)]
        df_motif_final = df_motif_adj.unstack()
        df_motif_final = df_motif_final[df_motif_final].astype(int).reset_index().rename(columns={'Gene_ID':'target','miRNA':'source',0:'weight'})

    df_motif_final = df_motif_final[['source','target','weight']]
    print('Final motif\n', df_motif_final.shape, df_motif_final.nunique())


    # Plots
    #plot_jaccard_vs_count(df_motif_adj_binary,df_mirna_family_adj,config.path_figures + 'Jaccard_motif.pdf')
    #plot_jaccard_vs_count(df_motif_adj.astype(int),df_mirna_family_adj_new,config.path_figures + 'Jaccard_motif_clean.pdf')


    # rename  the mirna family with the name of one of the mirnas
    pd.DataFrame(dict_mirna_merging.items(),columns=['mirna','family name']).to_csv(path_output+'miRNA_merging_mapping.csv',index=False)

    return df_motif_final,set(df_motif_adj.index),dict_mirna_merging
'''
def get_jaccard(df_adj_binary,cols_names = ['mirna1','mirna2']):
    """
    Get Jaccard similarity between rows of a binary matrix
    """
    df_jaccard = 1 - pairwise_distances(df_adj_binary.astype(bool).values, metric='jaccard')
    df_jaccard = pd.DataFrame(df_jaccard,index=df_adj_binary.index, columns=df_adj_binary.index)

    df_jaccard_list = df_jaccard.where(np.triu(np.ones(df_jaccard.shape),k=1).astype(bool))
    df_jaccard_list.index.name = 'mirna1'
    df_jaccard_list.columns = df_jaccard_list.columns.set_names('mirna2')
    df_jaccard_list = df_jaccard_list.stack().reset_index()
    df_jaccard_list.columns = cols_names+['Jaccard']
    return df_jaccard,df_jaccard_list

def plot_jaccard_vs_count(df_motif_adj_binary,df_mirna_family_adj,path_save_fig=None):
    df_mirna_motif_jaccard,df_mirna_motif_jaccard_list = get_jaccard(df_motif_adj_binary)
    tmp_sparse = csr_matrix(df_motif_adj_binary)
    df_count_co = pd.DataFrame((tmp_sparse*tmp_sparse.T).todense(), index = df_motif_adj_binary.index, columns=df_motif_adj_binary.index)
    df_count_diff = df_count_co*(1/df_mirna_motif_jaccard-1)
    df_count_diff_list=df_count_diff.unstack().dropna()
    df_count_diff_list.index.rename(['mirna1','mirna2'],inplace=True)
    df_count_diff_list = df_count_diff_list.reset_index().rename(columns={0:'n_diff_edges'})

    x = df_mirna_motif_jaccard.values[np.triu_indices(df_mirna_motif_jaccard.shape[0],k=1)]
    df_mirna_cofam = pd.DataFrame((np.dot(df_mirna_family_adj.T.fillna(0),df_mirna_family_adj.fillna(0))>0).astype(int),index=df_mirna_family_adj.columns, columns=df_mirna_family_adj.columns)
    y = df_mirna_cofam.loc[df_mirna_motif_jaccard.index,df_mirna_motif_jaccard.index].values[np.triu_indices(df_mirna_motif_jaccard.shape[0],k=1)]
    count_co =df_count_co.loc[df_mirna_motif_jaccard.index,df_mirna_motif_jaccard.index].values[np.triu_indices(df_mirna_motif_jaccard.shape[0],k=1)]


    plt.figure(figsize=[15,10])
    
    plt.subplot(2,3,1)
    plt.boxplot([x[y==0],x[y==1]])
    plt.title('Jaccard in cofam mirna and not')
    
    plt.subplot(2,3,2)
    plt.hist(x[y==1],bins=100)
    plt.title('Jaccard similarity between mirna in same families')
    
    plt.subplot(2,3,3)
    plt.scatter(x,count_co,c=y)
    plt.xlabel('Jaccard')
    plt.ylabel('N_common edges')
    plt.axhline(20)
    plt.title('Distr jaccard')
    
    plt.subplot(2,3,4)
    plt.scatter(x,count_co*(1/x-1),c=y)
    plt.xlabel('Jaccard')
    plt.ylabel('N_different edges')
    plt.axhline(20)

    plt.subplot(2,3,5)
    plt.title('Focus Jaccard>.5')
    plt.scatter(x[x>.5],(count_co*(1/x-1))[x>.5],c=y[x>.5])
    plt.xlabel('Jaccard')
    plt.ylabel('N_different edges')
    plt.axhline(50)
    if path_save_fig is not None:
        plt.savefig(path_save_fig)

def get_list_mirna(file_path):
    if not os.path.isfile(file_path):
        return []
    with open(file_path, "r") as f:
            list_miR = f.read().splitlines()
    return list_miR

def get_df_motif(file_path):
    df_motif = pd.read_csv(file_path, names=['source', 'target', 'weight'], sep='\t')
    return df_motif





### Pipeline Utils

def get_subjects_id(path_dict_processed_input):
    if not os.path.isfile(path_dict_processed_input):
        return []
    dict_input = get_processed_copdgene(path_dict_processed_input)
    df_exp = dict_input['df_exp']
    list_all = list(df_exp.columns)
    return sorted(list_all)[:10]
'''    

def run_lio_puma(subject_id,df_exp,np_corr_full,df_motif,list_mirna,path_file_output):
    # Compute Lioness Correlation
    np_corr_subset = np.corrcoef(df_exp.drop(subject_id, axis=1))
    n_conditions = df_exp.shape[1]
    np_lioness = n_conditions * (np_corr_full - np_corr_subset) + np_corr_subset
    df_corr = pd.DataFrame(np_lioness,index=df_exp.index,columns=df_exp.index)

    # the order of the genes is given by the expression file on which we computed the correlation before
    df_puma = run_puma_from_exp(None, df_motif, list_mirna,df_corr=df_corr)
    #later on will access following list_mirna id order
    df_puma.loc[list_mirna,df_exp.index].to_csv(path_file_output)

def run_puma_from_exp(df_exp, df_motif, miR,df_corr=None):
    puma = ntz.Puma(expression_file=df_exp, motif_file=df_motif, ppi_file=None, mir_file=miR,
                    precision='double', modeProcess='intersection',save_tmp=False,df_correlation_matrix=df_corr)
    df_puma = pd.DataFrame(puma.puma_network, index=puma.motif_tfs, columns=puma.motif_genes)
    return df_puma


def get_lm_coef_file_name(path_lm_coeff,target_pheno,mir):
    return os.path.join(path_lm_coeff,target_pheno,mir,'lm_coefficients.csv')   
        
'''
def _get_gsea_folder_name(path_lm_coeff,lm_regression,mir,target_pheno):
    return path_lm_coeff + '%s/%s/%s/' %(lm_regr_to_string(lm_regression),mir,target_pheno)    

def get_all_lm_coef_file_names(path_lm_coeff,lm_regressions,file_path_mirna):
    list_miR = get_list_mirna(file_path_mirna)
    list_file_name = [get_lm_coef_file_name(path_lm_coeff,lm_regression,mir) for mir in list_miR for lm_regression in lm_regressions]
    return   list_file_name

'''
# ray funtion to have progress bar
def to_iterator(obj_ids):
    while obj_ids:
        done, obj_ids = ray.wait(obj_ids)
        yield ray.get(done[0])
'''

# OLD: Lioness save triangular
def save_corr_triangular(matrix, file_output):
    # Extract upper triangular indices and values
    indices = np.triu_indices_from(matrix,1)
    values = matrix[indices]

    np.save(file_output,values)

def reconstruct_corr_matrix(vector):
    # Determine the size of the matrix
    offset = 1
    n = int((-1 + np.sqrt(1 + 8 * len(vector)))/ 2) + offset

    # Create an empty matrix, diagonal will be 1 because is correlation with itself
    matrix = np.ones((n, n))

    # Reconstruct the matrix from the vector
    matrix[np.triu_indices(n,1)] = vector
    matrix = matrix.T
    matrix[np.triu_indices(n,1)] = vector

    return matrix

def load_corr_matrix(filename):
    # Load the sparse matrix from file
    vector = np.load(filename)

    # Reconstruct the full matrix
    matrix = reconstruct_corr_matrix(vector)

    return matrix

'''
### DESEQ

def run_deseq(df_exp_raw,df_pheno,target_pheno,confounding_factors,continuous_factors,path_file_output=None):
    
    df_pheno = df_pheno[confounding_factors+[target_pheno]].dropna()
    # remove genes with 1 or more 0
    #df_exp_raw = df_exp_raw.loc[(df_exp_raw==0).sum(1)==0]
    print('confounding factors',confounding_factors)
    print('continuous_factors',continuous_factors)
    print('Target pheno',target_pheno)
    print('Path output',path_file_output)

    dds = DeseqDataSet(
        counts=df_exp_raw[df_pheno.index].T.astype(int),
        # set target pheno as the lastone to be used as contrast
        metadata=df_pheno,
        design_factors= confounding_factors+[target_pheno],
        continuous_factors = continuous_factors,
        refit_cooks=True,
        )
    dds.deseq2()
    import psutil
    stat_res = DeseqStats(dds)#,n_cpus=psutil.cpu_count(logical=False) - 1)
    
    if path_file_output is not None:
        with open(path_file_output+'.pckl', "wb") as f:
            pickle.dump(stat_res, f)
        
    stat_res.summary()
    if path_file_output is not None: 
        stat_res.results_df.to_csv(path_file_output + '_table.csv')
    return stat_res.results_df

def run_deseq_mrna_mirna(config, save_file=True,filter_genes=True):
    dict_input = get_processed_copdgene(config.path_dict_processed_input)
    df_exp, df_pheno,df_mirna_pheno,df_exp_mirna = dict_input['df_exp'], dict_input['df_pheno'], dict_input['df_pheno_mirna'],dict_input['df_exp_mirna']
    
    if config.deco_bcell: 
        df_exp_raw=dict_input['df_exp'].astype(int)

    else:
        df_exp_raw = pd.read_csv(config.path_file_exp_raw,index_col=0,sep='\t')
        df_exp_raw.index = df_exp_raw.index.str.split('.').str[0]

    if filter_genes:
        df_exp_raw = df_exp_raw.loc[df_exp.index,df_pheno.index]
    else:
        df_exp_raw = df_exp_raw[df_pheno.index]
    #df_exp_raw = df_exp_raw[df_exp_raw.min(axis=1)>0]
    confounding_factors = config.lm_confounder_columns 
    if save_file:
        path_file_output = config.path_deseq+'DESEQ_'+config.target_variable
        path_file_output_mirna = config.path_deseq+'DESEQ_mirna_'+config.target_variable
        
    else:
        path_file_output=path_file_output_mirna=None
    
    
    df_deseq = run_deseq(df_exp_raw,df_pheno,config.target_variable,confounding_factors+['Lc.Batch'],continuous_factors= config.continuous_factors,path_file_output = path_file_output)

    if config.cell_count_cols!=[]:
        df_deseq_cc = run_deseq(df_exp_raw,df_pheno,config.target_variable,confounding_factors+['Lc.Batch']+config.cell_count_cols,continuous_factors= config.continuous_factors+config.cell_count_cols,path_file_output = path_file_output+'_cell_count')

    if config.deco_bcell:
        # don'r run the mirna already run in the bulk
        return  df_deseq,df_deseq_cc,None,None


    df_exp_mirna_raw = pyreadr.read_r(config.path_mirna_exp)[None]
    df_exp_mirna_raw.index = df_exp_mirna_raw.index.map(mirna_str_proc)
    df_exp_mirna_raw.columns = df_exp_mirna_raw.columns.str[1:]
    # use only mirna that survived filtering on the mirna preproc.
    #print(df_exp_mirna_raw.head(),df_exp_mirna.head(),df_mirna_pheno.head())
    df_exp_mirna_raw = df_exp_mirna_raw.loc[df_exp_mirna.index,df_mirna_pheno.index]
    
    print('miRNA shape', df_exp_mirna_raw.shape)

    df_deseq_mirna = run_deseq(df_exp_mirna_raw,df_mirna_pheno,config.target_variable,confounding_factors +["mirSeqBatch"],continuous_factors= config.continuous_factors,path_file_output = path_file_output_mirna)
    if config.cell_count_cols!=[]:
        df_deseq_mirna_cc = run_deseq(df_exp_mirna_raw,df_mirna_pheno,config.target_variable,confounding_factors +["mirSeqBatch"]+config.cell_count_cols,continuous_factors= config.continuous_factors+config.cell_count_cols,path_file_output = path_file_output_mirna+'_cell_count')
        

    return df_deseq,df_deseq_cc,df_deseq_mirna,df_deseq_mirna_cc


def compute_pvalue_l2fc_groups(df_deseq,list_mirnas1,list_mirnas2,column):
    n1= len(list_mirnas1)
    n2= len(list_mirnas2)
    s0 = df_deseq[column]
    diff_abs = []

    for _ in range(10000):
        s1 = np.random.choice(s0,n1)
        s2 = np.random.choice(s0,n2)
        diff_abs.append(np.abs(s1.mean()-s2.mean()))

    obs_diff = np.abs(df_deseq.loc[list_mirnas1,column].mean()-df_deseq.loc[list_mirnas2,column].mean())
    diff_abs = np.array(diff_abs)
    p_val_emp = np.mean(obs_diff<=diff_abs)
    t_statistic, p_value = ttest_ind(df_deseq.loc[list_mirnas1,column], df_deseq.loc[list_mirnas2,column])
    
    print('Empirical T-test',p_val_emp,'T-test',p_value)
    _,p_value_1l = ttest_1samp(df_deseq.loc[list_mirnas1,column], 0,alternative='greater')
    _,p_value_1g = ttest_1samp(df_deseq.loc[list_mirnas1,column], 0,alternative='less')
    
    print('group 1 >0',p_value_1l,'<0',p_value_1g)

    _,p_value_2l = ttest_1samp(df_deseq.loc[list_mirnas2,column], 0,alternative='greater')
    _,p_value_2g = ttest_1samp(df_deseq.loc[list_mirnas2,column], 0,alternative='less')
    print('group 2 >0',p_value_2l,'<0',p_value_2g)

    return p_val_emp,p_value

def plot_hist_l2fc(df_deseq,list_mirnas1,list_mirnas2, column, name,with_all = False, deco=False):
    if name=='mRNA':
        group_name = 'G'
    elif name=='miRNA':
        group_name = 'm'
    if deco:
        group_name = 'b'+group_name 
    list_mirnas1_filt = [x for x in list_mirnas1 if x in df_deseq.index]    
    list_mirnas2_filt = [x for x in list_mirnas2 if x in df_deseq.index]   

    print("missing from %s1 %i %s and from %s2 %i %s" %(group_name,len(set(list_mirnas1).difference(list_mirnas1_filt)),name,group_name,len(set(list_mirnas2).difference(list_mirnas2_filt)),name))
    pvalue, pval_ttest= compute_pvalue_l2fc_groups(df_deseq,list_mirnas1_filt,list_mirnas2_filt,column)
 
    if with_all:
        plt.boxplot([df_deseq.reindex(list_mirnas1_filt).dropna()[column],df_deseq.reindex(list_mirnas2_filt).dropna()[column],df_deseq.dropna()[column]])
        plt.xticks([1,2,3],[group_name+'1',group_name+'2','all'])

    if not with_all:
        plt.boxplot([df_deseq.reindex(list_mirnas1_filt).dropna()[column],df_deseq.reindex(list_mirnas2_filt).dropna()[column]])
        plt.xticks([1,2],[group_name+'1',group_name+'2'])

    plt.title('%s fold change Boxplot\nSignificance:%.1E'%(name,pvalue),fontsize=15)
    plt.ylabel('log2FoldChange',fontsize=15)
    plt.xlabel(name+' group',fontsize=15)
    if with_all:
        plt.axhline(df_deseq[column].median(),c='orange',linestyle='--',linewidth=1)
    else:
        plt.axhline(df_deseq[column].median(),c='grey',linestyle='--',linewidth=1)
        
### Linear model coefficient analysis
def get_stats_all_edges(X,y,target_col):
    lm = LinearRegression()
    lm.fit(X,y)
    params = lm.coef_
    predictions = lm.predict(X)

    MSE = (sum((y-predictions)**2))/(X.shape[0]-X.shape[1])
    var_b = MSE.reshape(-1,1)*(np.linalg.inv(np.dot(X.T,X)).diagonal())
    sd_b = np.sqrt(var_b)
    ts_b = params/ sd_b

    p_values =[2*(1-stats.t.cdf(np.abs(i),(X.shape[0]-X.shape[1]))) for i in ts_b]
    return [params[:,target_col],sd_b[:,target_col],ts_b[:,target_col],np.array(p_values)[:,target_col]]
'''

def get_stats_all_edges_new(X,y,target_col):
    lm = LinearRegression()
    lm.fit(X,y)
    params = lm.coef_
    y_pred = lm.predict(X)
    #r2 = 1-((y - y_pred)** 2).sum(0)/ ((y - y.mean(0)) ** 2).sum(0)
    MSE = (sum((y-y_pred)**2))/(X.shape[0]-X.shape[1])
    var_b = MSE.reshape(-1,1)*(np.linalg.inv(np.dot(X.T,X)).diagonal())
    sd_b = np.sqrt(var_b)
    ts_b = params/ sd_b

    p_values =[2*(1-stats.t.cdf(np.abs(i),(X.shape[0]-X.shape[1]))) for i in ts_b]
    return np.array([params[:,target_col],sd_b[:,target_col],ts_b[:,target_col],np.array(p_values)[:,target_col]])#,r2


'''
### BCR DATA

def get_stats_all_edges_bcr(df_input,y_columns,formula_covariates):

    df_stat = pd.DataFrame()
    for y_column in y_columns:

        tmp = df_input.dropna(subset=[y_column]).copy()
        tmp[y_column] = tmp[y_column].astype(float)
        
        formula = 'standardize(%s) ~ '%y_column + formula_covariates
        Y, X = dmatrices(formula, tmp,NA_action='drop')
        # Get the slice names from patsy
        try:
            result = sm.OLS(Y, X).fit()
        except:
            print('OLS error')
            return pd.DataFrame()
        summary = result.summary().tables[1].as_html()
        summary = pd.read_html(summary, header=0, index_col=0)[0]
        summary = summary.reset_index().rename(columns={'index':'Variable','P>|t|':'pvalue'})
        summary.index = [y_column]*summary.shape[0]
        summary.pvalue = result.pvalues
        df_stat = pd.concat([df_stat,summary], axis=0)
    df_stat = df_stat.reset_index()
    
    adj_pval = df_stat.groupby('Variable').apply(lambda x: pd.Series(multipletests(x.pvalue,alpha=0.01,method='fdr_bh')[1], index = x.index))
    adj_pval = adj_pval.reset_index().set_index('level_1')[0]

    df_stat['adj_pval'] = adj_pval
    df_stat = df_stat.set_index('index')
    #df_stat.pvalue = df_stat.pvalue+10**(-15)
    #df_stat.adj_pval = df_stat.adj_pval+10**(-15)
    df_stat['neglogpval'] = df_stat.pvalue.map(lambda x: -np.log10(x))
    df_stat['neglogadjpval'] = df_stat.adj_pval.map(lambda x: -np.log10(x))

    
    return df_stat
'''
def get_stats_all_edges_bcr_old(df,confounders,y_columns, x_columns): 
    # Add constant
    df['Constant'] = 1
    confounders = confounders+['Constant']
    df = df[x_columns+confounders+y_columns].dropna()
    y_scaler = StandardScaler().fit(df[y_columns])
    y = y_scaler.transform(df[y_columns])
    
    df = df[x_columns+confounders]
    X = df[x_columns+confounders].copy()
    X = StandardScaler().fit_transform(X)
    X = pd.DataFrame(X, columns=df.columns)
    
    X.loc[:,df.columns[df.nunique()<=2]] = df[df.columns[df.nunique()<=2]].values
    
    lm = LinearRegression()
    
    lm.fit(X,y)
    params = lm.coef_ #np.append(lm.intercept_,lm.coef_)#
    #print(pd.DataFrame(params.T, index=x_columns+confounders))
    predictions = lm.predict(X)

    MSE = (sum((y-predictions)**2))/(X.shape[0]-X.shape[1])
    var_b = MSE.reshape(-1,1)*(np.linalg.inv(np.dot(X.T,X)).diagonal())
    sd_b = np.sqrt(var_b)
    ts_b = params/ sd_b
    predictions = pd.DataFrame(y_scaler.inverse_transform(predictions),index=df.index,columns =y_columns)
    p_values =np.array([2*(1-stats.t.cdf(np.abs(i),(X.shape[0]-X.shape[1]))) for i in ts_b])
    adj_p_values = np.array([multipletests(x,alpha=0.01,method='fdr_bh')[1] for x in p_values.T]).T
    
    
    lm_stats = [params,sd_b,ts_b,p_values,adj_p_values]
    # remove constant
    lm_stats = np.swapaxes(lm_stats, 0, 2)
    m,n,r = lm_stats.shape
    out_arr = np.column_stack((np.repeat(np.arange(m),n),lm_stats.reshape(m*n,-1)))
    df_stats = pd.DataFrame(out_arr,index = list(y_columns)*m,columns = ['Variable','Coefficient','Standard_Error','t','pvalue','adj_pval'])
    df_stats.Variable = df_stats.Variable.replace(dict(enumerate(x_columns+confounders)))
    df_stats = df_stats.sort_values('Coefficient')#,key=abs,ascending=False) 
    df_stats.pvalue = df_stats.pvalue+10**(-15)
    df_stats.adj_pval = df_stats.adj_pval+10**(-15)
    df_stats['neglogpval'] = df_stats.pvalue.map(lambda x: -np.log10(x))
    df_stats['neglogadjpval'] = df_stats.adj_pval.map(lambda x: -np.log10(x))

    df_stats_x = df_stats.query('Variable in @x_columns')
    df_stats_conf = df_stats.query('Variable in @confounders')
    
    return df_stats_x,df_stats_conf,predictions


def plot_bcr_term(tmp_var,df_bcr, ylabel_rename=None,gene='PAX5'):
    tmp = df_bcr.dropna(subset=[tmp_var]).copy()
    #sns.regplot(data=tmp, x='PAX5',y=tmp_var,label='all',ci=None)
    
    sns.regplot(data=tmp.query('fev1fvc_thr==0'), x=gene,y=tmp_var,label='Control',ci=None)
    sns.regplot(data=tmp.query('fev1fvc_thr==1'), x=gene,y=tmp_var,label='COPD',ci=None)
    if ylabel_rename is not None:
        plt.ylabel(ylabel_rename,fontsize=20)#(fontsize=30)
    plt.xlabel(gene,fontsize=20)#(fontsize=30)
    plt.tick_params(labelsize=15)
    plt.legend(fontsize=15)
    plt.tight_layout()
    #plt.savefig('../Tables/PAX5_BCR.jpeg', dpi=500)
    plt.show()
    df_bcr[gene+'_bin'] = pd.qcut(df_bcr[gene],4,labels=range(1,5))
    sns.boxplot(data=df_bcr, y=tmp_var,x=gene+'_bin',hue='fev1fvc_thr')
    if ylabel_rename is not None:
        plt.ylabel(ylabel_rename)
    #sns.regplot(data=df_bcr.query('fev1fvc_thr==1'), x='PAX5',y='totCount')
    plt.show()
    sns.boxplot(data=df_bcr, y=tmp_var,hue=gene+'_bin',x='fev1fvc_thr')
    #sns.regplot(data=df_bcr.query('fev1fvc_thr==1'), x='PAX5',y='totCount')
    if ylabel_rename is not None:
        plt.ylabel(ylabel_rename)
    plt.show()

def plot_bcr_term_reg(tmp_var,df_bcr, ylabel_rename=None,gene='PAX5',hue_col = 'gold'):
    tmp = df_bcr.dropna(subset=[tmp_var]).copy()
    #sns.regplot(data=tmp, x='PAX5',y=tmp_var,label='all',ci=None)
    for v in tmp[hue_col].unique():
        sns.regplot(data=tmp.query('%s==@v'%hue_col), x=gene,y=tmp_var,label=v,ci=None)
        
    if ylabel_rename is not None:
        plt.ylabel(ylabel_rename,fontsize=20)#(fontsize=30)
    else:
        plt.ylabel(tmp_var,fontsize=20)#(fontsize=30)
    plt.xlabel(gene,fontsize=20)#(fontsize=30)
    plt.tick_params(labelsize=15)
    plt.legend(fontsize=15)
    plt.tight_layout()
    #plt.savefig('../Tables/PAX5_BCR.jpeg', dpi=500)

'''
### Association to miRNA/mRNA  
def get_stats_lm_from_dict(dict_lm_stats,list_1,list_2,list_other,name_lm,group_name='m'):
    def get_list_coef(list_mirna,group_name):
        list_coeff = []
        for mirna in list_mirna:
            tmp_df = dict_lm_stats[mirna].reset_index()
            tmp_df.Variable = tmp_df.Variable.str.replace(mirna,name_lm) 
            tmp_df = tmp_df.set_index(['index',"Variable"])#[variable].rename(mirna)    
            tmp_df['reg'] = mirna
            list_coeff.append(tmp_df)
        
        df_coeff = pd.concat(list_coeff)#.T.reset_index()
        df_coeff['group'] = group_name
        
        return df_coeff

    df_coef_mirna1 = get_list_coef(list_1,group_name=group_name+'1')
    df_coef_mirna2 = get_list_coef(list_2,group_name=group_name+'2')
    df_coef_mirna_all = get_list_coef(list_other,group_name='other')
    df_coef_mirna = pd.concat([df_coef_mirna1,df_coef_mirna2,df_coef_mirna_all])
    return df_coef_mirna


### Phenotype Network
def get_pheno_specific_nt(config,set_mirna,thr_pval=1,thr_adj_pval=1,recompute=False):
    """
    Read through the results of the linear model, get the coefficients for the target phenotype with p_val>threshold_pvalue
    """
    path_outputfile = os.path.join(config.path_lm_coeff,config.target_variable, 'df_lm_coef_%s.csv'%(str(thr_adj_pval)))
    if not recompute and os.path.isfile(path_outputfile):
        df_lm_coef = pd.read_csv(path_outputfile)
        df_lm_adj = get_adjacency_matrix(df_lm_coef)
        return df_lm_coef,df_lm_adj
 
    df_lm_coef = pd.DataFrame()
    for mir in tqdm(set_mirna):
        df_stats = pd.read_csv(config.get_lm_coef_file_name(config.lm_exp,config.target_variable,mir),index_col=0)
        df_stats = df_stats.query('Variable==@config.target_variable and pvalues<=@thr_pval and adj_pval<=@thr_adj_pval').copy()
        df_stats['mirna'] = mir
        df_lm_coef = pd.concat([df_lm_coef,df_stats])
    ### If no results return them empty
    if df_lm_coef.shape[0]==0:
        return pd.DataFrame(),pd.DataFrame()
    df_lm_coef = df_lm_coef.reset_index().rename(columns={'index':'gene'})
    df_lm_coef['specific'] = (df_lm_coef.Coefficients>0).astype(int)
    dict_ens_gene = get_dict_ens_gene()
    df_lm_coef['gene_symbol'] = df_lm_coef.gene.apply(lambda x: dict_ens_gene[x])
    df_lm_coef.to_csv(path_outputfile,index=False)
    df_lm_adj = get_adjacency_matrix(df_lm_coef)

    return df_lm_coef,df_lm_adj

def plot_degrees(df_lm_pheno_spec,names={0:'control',1:'case'}, path_tables='./'):
    """
    Plot degree of bipartite network mirna-gene. Specific is 1 COPD (coeff>0) and 0 controls
    """
    df_degree_mirna = df_lm_pheno_spec.apply(lambda x: x.value_counts(),axis=1).rename(columns=names).fillna(0)
    df_degree_gene = df_lm_pheno_spec.apply(lambda x: x.value_counts()).T.rename(columns=names).fillna(0)
    
    # scatter plot degree of mirna and genes in Case and Control network
    _, ax = plt.subplots(1,3,figsize = [20,5])
    _, bins = np.histogram(df_degree_mirna.values.flatten(),bins=20)
    for tmp_value in names.values():
        df_degree_mirna[tmp_value].hist(bins=bins,label=tmp_value,alpha=.75,ax=ax[0])
    ax[0].legend()
    ax[0].set_title('miRNA distribution')

    
    def plot_scatter_density(df_degree,ax,name_subfig,names):
        case_name,control_name = names.values()
        x= df_degree[control_name].fillna(0).values
        y= df_degree[case_name].fillna(0).values
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        ax.scatter(x=x,y=y,c=z,label = 'Node degree')
        ax.set_xlabel(control_name,fontsize=15)
        ax.set_ylabel(case_name,fontsize=15)
        tmp_lim = max(ax.get_xlim()[1],ax.get_ylim()[1])
        ax.plot([0,tmp_lim],[0,tmp_lim],label='x=y',linestyle='--')
        ax.legend()
        ax.set_title(name_subfig,fontsize=15)
        
    plot_scatter_density(df_degree_mirna,ax[1],'miRNAs degree',names)
    plot_scatter_density(df_degree_gene,ax[2],'Genes degree',names)

    dict_ens_gene = get_dict_ens_gene()
    df_degree_gene.loc[df_degree_gene.sum(1).sort_values(ascending=False).index].to_csv(os.path.join(path_tables,'Degree_bulk_mrna.csv'))
    df_degree_mirna.loc[df_degree_mirna.sum(1).sort_values(ascending=False).index].to_csv(os.path.join(path_tables,'Degree_bulk_mirna.csv'))


    plt.show()
    return df_degree_mirna,df_degree_gene

def plot_adjacency_matrix(df_adj,zmin_max = (None,None),figsize = [1400,500],all_labels=False,color_continuous_scale = 'rdbu'):
    layout = go.Layout(
        title='some title',
        yaxis=dict(
            title='Xaxis Name',
            tickmode='linear'))#height=2000

    img = px.imshow(df_adj,aspect='auto',color_continuous_scale = color_continuous_scale,zmin=zmin_max[0],zmax=zmin_max[1],width = figsize[0],height=figsize[1])#

    if all_labels :
        img.update_layout(yaxis={"dtick":1},xaxis={"dtick":1},margin={"t":0,"b":0})#,height=1000)
    img.show()

    # 67 miRNAs and 2470 genes were associated with at least one of these disease-specific edges
    print(df_adj.shape)  
    return img

def get_adjacency_matrix(df_lm_coef):
    df_lm_adj = df_lm_coef.set_index(['mirna','gene']).specific.unstack() #.query('adj_pval<0.1')
    df_lm_adj = df_lm_adj.loc[df_lm_adj.notna().sum(1).sort_values(ascending=False).index,df_lm_adj.notna().sum(0).sort_values(ascending=False).index]
    return df_lm_adj
'''
def tsnr(x,y):
        return (np.mean(x)-np.mean(y))/(np.sqrt(np.std(x,ddof=1)/len(x) + np.std(y,ddof=1)/len(y)))  

def t_test_mirna_target(tmp_mirna,df_stat,df_lm_adj,df_motif_adj,label0='',label1='',plot=False,path_figure=None):
    """
    Positive tsnr meand that genes target by that mirna are on average more expressed in COPD. 
    """
    
    degree = df_lm_adj.loc[tmp_mirna].value_counts()
    # in case there are no edges for label0 or label1
    degree = degree.reindex([0,1]).fillna(0)
    tmp_edges1 = df_lm_adj.loc[tmp_mirna].pipe(lambda x: x[x==1]).index
    tmp_edges0 = df_lm_adj.loc[tmp_mirna].pipe(lambda x: x[x==0]).index
    tmp_prior_edges_genes = df_motif_adj.loc[tmp_mirna,df_stat.index].pipe(lambda x: x[x==1]).index

    t_test1 = df_stat.loc[tmp_edges1]
    t_test0 = df_stat.loc[tmp_edges0]
    t_test_prior = df_stat.loc[tmp_prior_edges_genes]

    t_test_t_test_stat,t_test_t_test = ttest_ind(t_test1,t_test0)
    tsnr_score = tsnr(t_test1,t_test0)
    if plot:
        plt.hist(t_test1,label='Genes from %s edges'%label0,alpha=.7,bins=30)
        plt.hist(t_test0,label='Genes from %s edges'%label0,alpha=.7,bins=30)
        plt.hist(t_test_prior,label='Genes from Prior edges',alpha=.7,bins=30)
        
        plt.title('Differential expression t-test distributions\nGene targets of '+tmp_mirna +', pval:' +"{:.1e}".format(t_test_t_test)+', tsnr:' +"{:.1e}".format(tsnr_score)+
                 '\n N edges: %i %s, %i %s, %i Motif'%(degree[0],label0,degree[1],label1,len(tmp_prior_edges_genes)))
        plt.xlabel('T-statistic')
        plt.ylabel('N genes')
        plt.legend()
        if path_figure is not None:
            plt.savefig(path_figure,bbox_inches='tight')
        plt.show()
        plt.clf()
        
    return t_test_t_test_stat,t_test_t_test,degree[0],degree[1],len(tmp_prior_edges_genes),tsnr_score

@ray.remote
def ray_t_test_mirna_target(df_exp,df_lm_adj,df_motif_adj,tmp_mirna,list_case,list_control,plot=False):
    return t_test_mirna_target(df_exp,tmp_mirna,list_case,list_control,plot=plot)

def genes_targeted_diff_exp(df_exp,df_motif_adj,df_lm_adj,list_mirna,list_case,list_control,output_file=None,plot=False):

    #df_exp = map_ens_to_genesymb(df_exp)
    #df_motif_adj = map_ens_to_genesymb(df_motif_adj)
    #df_lm_adj = map_ens_to_genesymb(df_lm_adj)

    id_df_exp = ray.put(df_exp)
    id_df_motif_adj = ray.put(df_motif_adj)
    id_df_lm_adj = ray.put(df_lm_adj)

    res_id = [ray_t_test_mirna_target.remote(id_df_exp,id_df_motif_adj,id_df_lm_adj,tmp_mirna,list_case,list_control,
                                             id_df_exp,plot=plot) for  tmp_mirna in tqdm(list_mirna)]

    df_results = []
    for x in tqdm(to_iterator(res_id), total=len(res_id)):
        df_results.append(x)
    
    df_results = pd.DataFrame(df_results,columns=['mirna','T-stat','Pvalue','degree_Control','degree_COPD','degree_prior','tsnr']).set_index('mirna')
    try:
        df_results['pval_adj'] = multipletests(df_results.Pvalue,alpha=0.01,method='fdr_bh')[1]
        df_results['Pvalue_log'] = -np.log10(df_results.Pvalue)
        df_results['pval_adj_log'] = -np.log10(df_results.pval_adj)
    except:
        pass
    if output_file is not None:
        df_results.to_csv(output_file)
    return df_results

'''
def boxplot_corr(df,list1,list2,list1b=None,list2b=None,list_rnd1=None,list_rnd2=None,title='',omics='',df_not_filt=None):
    
    def get_triu(df,row_group1,col_group1):
        a = df.loc[row_group1,col_group1].values

        if len(row_group1)== len(col_group1) and row_group1==col_group1:
            a = a[np.triu_indices_from(a,k=1)]
        else:
            a = a.flatten()
        return a
    
    
    if list1b is None and list2b is None:
        list1b = list1
        list2b = list2
    list_plot = [] 
    list_plot.append(get_triu(df,list1,list1b))
    list_plot.append(get_triu(df,list2,list2b))
    list_plot.append(get_triu(df,list1,list2b))
    if len(list1b)!= len(list1) or list1b != list1:
        list_plot.append(get_triu(df,list2,list1b))
    if list_rnd1 is None: 
        #list_rnd1=list_rnd2 = list(np.random.choice(df.index,max(len(list1),len(list2)),replace=False))
        list_rnd1=list_rnd2 = list(df.index.difference(list1+list2))
        if len(omics)==2:
            list_rnd2 = list(df.columns.difference(list1+list2))

    if df_not_filt is not None:
        # when df input is filtered (maybe only correlation of pairs with an edge in the network), for the comparison with the "other" group you want to use the unfiltered
        e = get_triu(df_not_filt,list_rnd1,list_rnd2)
    else:
        e = get_triu(df,list_rnd1,list_rnd2)

    list_plot.append(e)

    
    ax = sns.boxplot(data=list_plot)
    if len(list1b)!= len(list1) or list1b != list1:
        ax.set_xticklabels(['%s1-%s1*'%(omics[0],omics[-1]),'%s2-%s2*'%(omics[0],omics[-1]),'%s1-%s2*'%(omics[0],omics[-1]),'%s2-%s1*'%(omics[0],omics[-1]),'Other'])
    else:
        ax.set_xticklabels(['%s1-%s1*'%(omics[0],omics[-1]),'%s2-%s2*'%(omics[0],omics[-1]),'%s2-%s1*'%(omics[0],omics[-1]),'Other'])
    ax.tick_params(axis='both', labelsize=15)
    plt.axhline(0,c='grey',linestyle='--')
    plt.title(title,fontsize=15)


    
'''
def pers_pagerank(network,personalization = None, max_iter=100, tol=1.0e-6, dangling=None):

    M = nx.to_scipy_sparse_array(network)
    nodelist = network.nodes()

    S = np.array(M.sum(axis=1), dtype=int).flatten()
    S = np.divide(1, S, out=np.zeros(network.number_of_nodes()), where= S!=0)
    is_dangling = np.where(S == 0)[0]
    Q = scipy.sparse.spdiags(S.T, 0, *M.shape, format='csr')
    M = Q * M
    M = M.transpose()
    N =len(nodelist)

    # initial vector
    x = np.repeat(1.0 / N, N)

    # Personalization vector
    if personalization is None:
        p = np.repeat(1.0 / N, N)
    else:
        missing = set(nodelist) - set(personalization)
        if missing:
            raise nx.exception.NetworkXError('Personalization vector dictionary must have a value for every node. '
                                             'Missing nodes %s' % missing)
        p = scipy.array([personalization[n] for n in nodelist], dtype=float)
        p = p / p.sum()

    # Dangling nodes
    if dangling is None:
        dangling_weights = p
    else:
        missing = set(nodelist) - set(dangling)
        if missing:
            raise nx.exception.NetworkXError('Dangling node dictionary must have a value for every node. Missing '
                                             'nodes %s' % missing)
        # Convert the dangling dictionary into an array in nodelist order
        dangling_weights = scipy.array([dangling[n] for n in nodelist], dtype=float)
        dangling_weights /= dangling_weights.sum()
    # print(x,p,M.sum())
    # power iteration: make up to max_iter iterations
    for _ in range(max_iter):
        xlast = x
        x = alpha * (M@x + sum(x[is_dangling]) * dangling_weights) + (1 - alpha) * p
        # check convergence, l1 norm
        err = np.absolute(x - xlast).sum()
        if err < N * tol:

            return dict(zip(nodelist, map(float, x))) # ,x,M)
    raise nx.exception.NetworkXError('power iteration failed to converge in %d iterations'% max_iter )
'''
def get_nt_groups(df_lm_adj,thr_mirna=1,thr_gene=1,path_nt_groups=None, plot=True, fisher_test=True,m1_bulk=None):
    """
    Compute the 4 groups, m1,m2,g1,g2. 
    If they are not perfected the nodes are assigned based on the highest number of edges in common
    
    """
    df = df_lm_adj.copy()
    while (df.sum()<thr_gene).sum() + (df.sum(1)< thr_mirna).sum()>0:
        df = df.loc[df.sum(1)>=thr_mirna, df.sum(0)>=thr_gene]

    list_genes1,list_genes2,list_mirnas1,list_mirnas2 = get_mirna_gene_groups(df)
    dict_ens_gene = get_dict_ens_gene()
    dict_gene_ens = {y:x for x,y in dict_ens_gene.items()} 
    
    if dict_gene_ens['PAX5'] not in list_genes2:
        # Fix pax5 in G2 to be coherent with the paper draft
        a =list_genes2
        list_genes2 = list_genes1
        list_genes1 = a
    if m1_bulk is not None and len(set(list_mirnas1).intersection(m1_bulk))<len(set(list_mirnas2).intersection(m1_bulk)):
        # If there's a reference group for miRNAs group (like the bulk analysis when doing the b-cell), make sure the groups are aligned
        a = list_mirnas2.copy()
        list_mirnas2 = list_mirnas1
        list_mirnas1 = a
        
    df_lm_adj = df_lm_adj.loc[list_mirnas1 + list_mirnas2,list_genes1 + list_genes2]

    if path_nt_groups is not None:

        dict_groups = {'list_genes1':list_genes1,
                       'list_genes2':list_genes2,
                       'list_mirnas1':list_mirnas1,
                       'list_mirnas2':list_mirnas2}
        with open(path_nt_groups,'w') as f:
            json.dump(dict_groups,f)

        
    if plot:
        tmp_plot = df_lm_adj*2-1
        #tmp_plot = tmp_plot.append(df_deseq.loc[tmp_plot.columns].log2FoldChange)
        #tmp_plot = tmp_plot.merge(df_deseq_mirna.log2FoldChange,left_index=True, right_index=True,how='left')
        tmp_plot = tmp_plot.rename(columns=dict_ens_gene).dropna(how='all').dropna(how='all',axis=1)

        print('Number of interactions:',df_lm_adj.notna().sum().sum())
        print('Number of mirnas and genes',tmp_plot.shape)
        print('m1:',len(list_mirnas1),', m2:',len(list_mirnas2),', G1:',len(list_genes1),', G2:',len(list_genes2))

        plot_adjacency_matrix(tmp_plot,figsize = [1000,500])

    if fisher_test:
        # Fisher Test
        n_edg_m1_g1= df_lm_adj.loc[list_mirnas1,list_genes1].notna().sum().sum()
        n_edg_m1_g2= df_lm_adj.loc[list_mirnas1,list_genes2].notna().sum().sum()
        n_edg_m2_g1= df_lm_adj.loc[list_mirnas2,list_genes1].notna().sum().sum()
        n_edg_m2_g2= df_lm_adj.loc[list_mirnas2,list_genes2].notna().sum().sum()

        table = [[n_edg_m1_g1,n_edg_m1_g2],[n_edg_m2_g1,n_edg_m2_g2]]

        from scipy.stats import fisher_exact
        res = fisher_exact(table, alternative='greater')
        print('Fisher Test on the number of edges between the 4 groups')
        print(res)
        print(pd.DataFrame(table,index = ['M1','M2'],columns = ['G1','G2']))
        
    return df_lm_adj,list_genes1,list_genes2,list_mirnas1,list_mirnas2

    
def pair_wise_hamming_mask_nan(i):
    """
    Compute the hamming distance between two vectors. Drop nans.
    [:,None] allows you to cast the operation ;)
    """

    mask_pairwise = ~np.isnan(i) & ~np.isnan(i)[:, None]

    sum_mask = np.nansum(mask_pairwise,axis=2)
    #sum_mask[sum_mask==0] = 2

    d = np.nansum(~(i == i[:, None]) & mask_pairwise, axis=2)
    #if no intersection, sum_mask will be 0. Thus division will give nan. replace with 1: max distance between those two entries because no intersection
    d = d/sum_mask
    d = np.nan_to_num(d,nan=1)
    
    return d 

def get_mirna_gene_groups(df):
    #Get genes groups
    agg_clustering = AgglomerativeClustering(n_clusters=2,affinity='precomputed', linkage='average')
    cluster = agg_clustering.fit_predict(pair_wise_hamming_mask_nan(df.values))

    list_mirnas1 = list(df.index[cluster==0])
    list_mirnas2 = list(df.index[cluster==1])
    list_mirnas1 = list(sorted(list_mirnas1, key = lambda x: -(df==1).sum(1)[x]))
    list_mirnas2 = list(sorted(list_mirnas2, key = lambda x: -(df==0).sum(1)[x]))
        
    cluster = agg_clustering.fit_predict(pair_wise_hamming_mask_nan(df.T.values))
    list_genes1 = list(df.columns[cluster==0])
    list_genes2 = list(df.columns[cluster==1])
    list_genes1 = list(sorted(list_genes1, key = lambda x: -(df==1).sum(0)[x]))
    list_genes2 = list(sorted(list_genes2, key = lambda x: -(df==0).sum(0)[x]))

    return list_genes1,list_genes2,list_mirnas1,list_mirnas2



def independent_corr(xy, ab, n, n2 = None, twotailed=True):
    """
    Calculates the statistic significance between two independent correlation coefficients
    @param xy: correlation coefficient between x and y
    @param xz: correlation coefficient between a and b
    @param n: number of elements in xy
    @param n2: number of elements in ab (if distinct from n)
    @param twotailed: whether to calculate a one or two tailed test, only works for 'fisher' method
    @return: z and p-val
    """
    xy_z = 0.5 * np.log((1 + xy)/(1 - xy))
    ab_z = 0.5 * np.log((1 + ab)/(1 - ab))
    if n2 is None:
        n2 = n

    se_diff_r = np.sqrt(1/(n - 3) + 1/(n2 - 3))
    diff = xy_z - ab_z
    z = abs(diff / se_diff_r)
    p = (1 - norm.cdf(z))
    if twotailed:
        p *= 2

    return z, p

def par_corr_between(data,x, y, z):
    """
    Partial correlation coefficients between two variables x and y with respect to
    the control variable z.
    return the partial correlation coefficient
    """
    x = data[x]
    y = data[y]
    z = data[z]

    if z.ndim == 1:
        z = np.reshape(z, (-1, 1))
    # solve two linear regression problems Zw = x and Zw = y
    Z = np.hstack([z, np.ones((z.shape[0], 1))])  # bias
    wx = np.linalg.lstsq(Z, x, rcond=None)[0]
    rx = x - Z @ wx # residual
    wy = np.linalg.lstsq(Z, y, rcond=None)[0]
    ry = y - Z @ wy
    # compute the Pearson correlation coefficient between the two residuals
    if rx.shape[1]==1:
        #IDKW but if it has only 1 column you have to give it as series not DataFrame
        rx=rx.iloc[:,0]
    pcorr = pd.DataFrame(ry.corrwith(rx))
    return pcorr 


'''

def import_seaborn():
    # https://notebooks.githubusercontent.com/view/ipynb?color_mode=auto&commit=285350bb44d9ab1c61f829c1664f59279ff2c774&enc_url=68747470733a2f2f7261772e67697468756275736572636f6e74656e742e636f6d2f676973742f616672656e646569726f2f35343135383532626365646339333134623130353865356439333266356332312f7261772f323835333530626234346439616231633631663832396331363634663539323739666632633737342f736561626f726e5f636c75737465726d61705f776974685f6566666f72746c6573735f636f6c6f72626172735f6465636f7261746f722e6970796e62&logged_in=false&nwo=afrendeiro%2F5415852bcedc9314b1058e5d932f5c21&path=seaborn_clustermap_with_effortless_colorbars_decorator.ipynb&repository_id=103071526&repository_type=Gist
    from typing import List, Union, Optional

    import numpy as np
    import pandas as pd
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns


    SEQUENCIAL_CMAPS = [
        "Greys", "Purples", "Blues", "Greens", "Oranges", "Reds",
        "YlOrBr", "YlOrRd", "OrRd", "PuRd", "RdPu", "BuPu",
        "GnBu", "PuBu", "YlGnBu", "PuBuGn", "BuGn", "YlGn",
        "binary", "gist_yarg", "gist_gray", "gray", "bone", "pink",
        "spring", "summer", "autumn", "winter", "cool", "Wistia",
        "hot", "afmhot", "gist_heat", "copper","coolwarm"]

    def minmax_scale(x):
        return (x - x.min()) / (x.max() - x.min())
    def get_cmap_cent_0(x,cmap):
        abs_max = np.max(np.abs(x))
        list_cmap = plt.get_cmap(cmap)([-abs_max]+list(x)+[abs_max]).tolist()
        return list_cmap[1:-1]


    def to_color_series(x: pd.Series, cmap: str = "Greens") -> pd.Series:
        # return pd.Series(plt.get_cmap(cmap)(minmax(x)).tolist(), index=x.index, name=x.name).apply(matplotlib.colors.to_hex)
        return pd.Series(
            get_cmap_cent_0(x,cmap), index=x.index, name=x.name
            # plt.get_cmap(cmap)(x).tolist(), index=x.index, name=x.name
        )  # .apply(tuple)

    def to_color_dataframe(x: Union[pd.Series, pd.DataFrame], cmaps: Optional[Union[str, List[str]]] = None) -> pd.DataFrame:
        if isinstance(x, pd.Series):
            x = x.to_frame()
        if cmaps is None:
            cmaps = [plt.get_cmap(cmap) for cmap in SEQUENCIAL_CMAPS[: x.shape[1]]]
        return pd.concat([to_color_series(x[col], cmap) for col, cmap in zip(x, cmaps)], axis=1)

    def add_extra_colorbars_to_clustermap(
        grid: sns.matrix.ClusterGrid,
        datas: Union[pd.Series, pd.DataFrame],
        cmaps: Optional[Union[str, List[str]]] = None,
        location="cols",
        **kwargs
    ) -> None:
        def add(data: pd.Series, cmap: str, bbox: List[List[int]], orientation: str) -> None:
            ax = grid.fig.add_axes(matplotlib.transforms.Bbox(bbox))
            norm_max = data.abs().max()
            norm = matplotlib.colors.Normalize(vmin=-norm_max, vmax=norm_max)
            cb1 = matplotlib.colorbar.ColorbarBase(
                ax, cmap=plt.get_cmap(cmap), norm=norm, orientation=orientation, label=None#data.name
            )

        if isinstance(datas, pd.Series):
            datas = datas.to_frame()
        if cmaps is None:
            cmaps = SEQUENCIAL_CMAPS[: datas.shape[1]]
        if isinstance(cmaps, str):
            cmaps = [cmaps]

        # get position to add new axis in existing figure
        # # get_position() returns ((x0, y0), (x1, y1))
        heat = grid.ax_heatmap.get_position()
        cbar_spacing = .05
        cbar_size = 0.2
        if location in ["col", "cols", "columns"]:
            orientation = "horizontal"
            dend = grid.ax_col_dendrogram.get_position()
            y0 = dend.y0
            y1 = dend.y1
            for i, (data, cmap) in enumerate(zip(datas, cmaps)):
                if i == 0:
                    x0 = heat.x1 
                    x1 = heat.x1 + cbar_size 
                else:
                    x0 += cbar_size + cbar_spacing
                    x1 += cbar_size + cbar_spacing
                add(datas[data], cmap, [[x0, y0+.06], [x1, y0+.08]], orientation)
        else:
            orientation = "vertical"
            dend = grid.ax_row_dendrogram.get_position()
            x0 = dend.x0
            x1 = dend.x1
            for i, (data, cmap) in enumerate(zip(datas, cmaps)):
                if i == 0:
                    y0 = dend.y0 - cbar_size
                    y1 = dend.y0 
                else:
                    y0 -= cbar_size + cbar_spacing
                    y1 -= cbar_size + cbar_spacing
                add(datas[data], cmap, [[x1-.08, y0], [x1-0.06, y1]], orientation)
                
    def add_colorbars(
            grid,
            rows: pd.DataFrame = None,
            cols: pd.DataFrame = None,
            row_cmaps: Optional[List[str]] = None,
            col_cmaps: Optional[List[str]] = None
    ) -> None:
        if rows is not None:
            add_extra_colorbars_to_clustermap(
                grid, rows, location="row", cmaps=row_cmaps)
        if cols is not None:
            add_extra_colorbars_to_clustermap(
                grid, cols, location="cols", cmaps=col_cmaps)
            
    def colorbars(f):
        """The actual decorator of seaborn.clustermap"""
        def clustermap(*args, **kwargs):
            cmaps = {"row": None, "col": None}
            # capture "row_cmaps" and "col_cmaps" out of the kwargs
            for arg in ['row', 'col']:
                if arg + "_colors_cmaps" in kwargs:
                    cmaps[arg] = kwargs[arg + '_colors_cmaps']
                    del kwargs[arg + '_colors_cmaps']
            # get dataframe with colors and respective colormaps for rows and cols
            # instead of the original numberical values
            _kwargs = dict(rows=None, cols=None)
            for arg in ['row', 'col']:
                if arg + "_colors" in kwargs:
                    if isinstance(kwargs[arg + "_colors"], (pd.DataFrame, pd.Series)):
                        _kwargs[arg + "s"] = kwargs[arg + "_colors"]
                        kwargs[arg + "_colors"] = to_color_dataframe(kwargs[arg + "_colors"], cmaps[arg])
            grid = f(*args, **kwargs)
            add_colorbars(grid, **_kwargs, row_cmaps=cmaps['row'], col_cmaps=cmaps['col'])
            return grid
        return clustermap

    from importlib import reload

    reload(sns)
    sns.clustermap = colorbars(sns.clustermap)
    return sns

'''