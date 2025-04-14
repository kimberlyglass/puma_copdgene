
import os


folder_name='PUMA_COPDGene'
user='remge'


class Config:

    def __init__(self,expression,folder_name=folder_name, user=user) -> None:
        if expression not in ['bulk','bcell']:
            raise ValueError(f"Invalid input value: {expression}. Must be one of ['bulk','bcell'] ")

        self.folder_name = folder_name
        self.user = user

        self.expression  = expression
        # added filter on bcell expression each gene has to have exp>10 in at least 10% pz, also updated deco from min 
        self.path_output = '/proj/regeps/regep00/studies/COPDGene/analyses/%s/%s/output_%s/'%(user,self.folder_name,expression)
        self.path_data = '/proj/regeps/regep00/studies/COPDGene/analyses/%s/%s/data/'%(user,self.folder_name)
        self.path_lio_puma = '/d/tmp/%s/%s/%s/Lioness_Puma/'%(user,self.folder_name,expression) 
        
        # if you want to binarize the miRNA-mRNA motif interactions based on a threshold of the TargetScan score. Cont: means keep the score as is.
        self.motif_threshold_score = 'cont'#-.15# 
        # scanner type and smoking for emphyzema
        self.lm_confounder_columns = ["gender", "Age_P2", "race"]#,"BMI_P2"] #,,'' #ATS_PackYears_P2 "BMI_P2","ATS_PackYears_P2","BMI_P2","ATS_PackYears_P2"
        #lm_confounder_columns =lm_confounder_columns + ['ATS_PackYears_P2','wbc_P2','neutrophl_pct_P2','lymphcyt_pct_P2','monocyt_pct_P2']#,'eosinphl_pct_P2']
        self.cell_count_cols  = ['wbc_P2','lymphcyt_pct_P2']#['wbc_P2','neutrophl_pct_P2']#[','wbc_P2','neutrophl_pct_P2']#['wbc_P2']#,'eosinphl_pct_P2']'neutrophl_pct_P2'
        self.target_variable = "smoking_status_P2"#"FEV1_FVC_post_P2"
        self.continuous_factors = ['Age_P2','FEV1_FVC_post_P2',"smoking_status_P2","BMI_P2"] #,'BMI_P2' ,"ATS_PackYears_P2","BMI_P2","ATS_PackYears_P2"
        self.lm_randomize = False


        #Expression Files
        path_copdgene_freeze4 = '/proj/regeps/regep00/studies/COPDGene/data/rna/mrna/blood_shortread/data/freezes/freeze4/'
        path_copdgene_freeze4_salmon = path_copdgene_freeze4 + 'salmon/GRCh38.p13_GENCODE.v37/salmon_tximport/salmon_gene/matrix/'
        self.path_mrna_master = path_copdgene_freeze4 + 'masterfile/master.file.freeze4.txt'
        self.path_exp = path_copdgene_freeze4_salmon+ 'log-CPM_from_scaledTPM_LC.tsv'
        path_exp_bcell = '/proj/edith/regeps/regep00/studies/COPDGene/analyses/remhr/RNAseq_Deconvolution/Results/COPDGene_Phase2_Freeze5_HiRes_LM22/20240328/CIBERSORTxHiRes_NA_Bcells_Window40.txt'
        self.path_file_exp_raw = path_copdgene_freeze4_salmon + 'counts_raw_from_lengthScaledTPM.tsv'
        self.path_mrna_pheno =  '/proj/regeps/regep00/studies/COPDGene/data/pheno/COPDGene_P1P2P3_Flat_SM_NS_Oct23.txt'


        #not using Brian anymore, implemented tmm and cpm with pz and mirna filtering
        #path_mirna_brian_proc = '/udd/remge/Projects/PUMA/Data/miRNA_from_hobbs_rmd_cpm_from_tmm.csv'
        self.path_mirna_exp =  '/proj/regeps/regep00/studies/COPDGene/analyses/rebdh/craigMiRna/mirSeqCounts_filt_blockFreeze2_20200608.rds'
        self.path_mirna_pheno =  '/proj/regeps/regep00/studies/COPDGene/analyses/rebdh/craigMiRna/mirSeqPheno_filt_blockFreeze2_20200608.rds'


        self.path_targetscan = '/proj/rerefs/reref00/TargetScan/Nonconserved_Site_Context_Scores_hsa.csv' #self.path_data + 'motif/Nonconserved_Site_Context_Scores_hsa.csv'#'../Data/Summary_Counts.all_predictions.sel.cols.hsa.txt'


        if expression=='bcell':
            self.path_exp = path_exp_bcell
            self.deco_bcell= True
        else:
            self.deco_bcell= False


        
        self.path_puma_sbj = self.path_lio_puma +'by_subj/'
        self.path_puma_mirna = self.path_lio_puma + 'by_mirna/'
        self.path_all_puma = self.path_lio_puma + 'Puma_Network_all.csv'


        self.path_motif = self.path_output+ "df_motif_cleaned.csv"
        self.path_exp_cleaned = self.path_output+ "df_exp_cleaned.csv"
        self.path_exp_mirna_cleaned = self.path_output+ "df_exp_mirna_cleaned.csv"
        self.path_mirna = self.path_output+ "List_mirna_cleaned.csv"
        self.path_file_co_exp  = self.path_output+ "Df_corr_all.csv"
            

        self.path_dict_processed_input = self.path_output+'dict_input_COPDGene.pickle'
        self.path_sbj_list = self.path_output + 'list_sbj.csv'
        self.path_mirna_list = self.path_output + 'list_mirna.csv'
        self.path_gene_list = self.path_output + 'list_gene.csv'
        self.path_sbj_list = self.path_output + 'list_sbj.csv'
        self.path_figures = self.path_output + 'Figures/'
        self.path_tables = self.path_output + 'Tables/'#'/udd/remge/Projects/PUMA/Tables/'
        self.path_lm_dir = self.path_output +'LM_analysis/'
        self.lm_exp = 'base'
        self.path_lm_coeff = self.get_lm_folder(self.lm_exp)
        self.path_nt_groups = os.path.join(self.path_lm_coeff,self.target_variable,'Dict_mirna_gene_group.json')


        self.path_to_ens_mapping = self.path_data+'mapping/mart_export_grch38.txt'
        self.path_mapping_t_g = self.path_data+'mapping/mapping_ENSG_ENST.txt'
        self.path_file_annot = self.path_data+'mapping/gencode.v37.annotation.csv'

        self.path_deseq = os.path.join(self.path_output,'DESEQ', '')

        os.makedirs(self.path_output,exist_ok=True)
        os.makedirs(self.path_figures,exist_ok=True)
        os.makedirs(self.path_tables,exist_ok=True)
        os.makedirs(self.path_puma_mirna,exist_ok=True)
        os.makedirs(self.path_puma_sbj,exist_ok=True)
        os.makedirs(self.path_lm_coeff,exist_ok=True)
        os.makedirs(self.path_deseq,exist_ok=True)

        if self.expression =='bcell':
            # add path for BCR analysis
            self.path_data_bcr_file = os.path.join(self.path_data,'BCR', 'pheno_airr2.Rdata')
            self.path_data_bcr_lm = os.path.join(self.path_output,'BCR')
            os.makedirs(self.path_data_bcr_lm,exist_ok=True)

    def get_lm_folder(self,lm_exp):
        path_lm_coeff = os.path.join(self.path_lm_dir,lm_exp)
        return path_lm_coeff

    def update_lm(self,lm_exp):
        self.path_lm_coeff = self.get_lm_folder(lm_exp)
        if lm_exp=='base':
            return
        if 'cc' in lm_exp:
            self.lm_confounder_columns = self.lm_confounder_columns + self.cell_count_cols
        if 'py' in lm_exp:
            self.lm_confounder_columns = self.lm_confounder_columns + ['ATS_PackYears_P2']
        if lm_exp=='random':
            self.lm_randomize = True
        self.lm_exp=lm_exp

    
    def get_lm_coef_file_name(self,lm_exp,target_variable,mir):
        return os.path.join(self.get_lm_folder(lm_exp),target_variable,mir,'lm_coefficients.csv')   
    
    def get_lm_bcr_gene(self,gene):
        return os.path.join(self.path_data_bcr_lm,gene+'.csv')