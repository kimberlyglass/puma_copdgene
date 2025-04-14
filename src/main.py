from modules.Process_data import run_preprocess
from modules.Run_deseq import run_deseq
#from modules.BCR_mRNA_groups_lm import 
from config import Config


# create files for the bulk analysis
config_bulk = Config('bulk')
run_preprocess(config_bulk)

# create files for the bcell analysis
config_bcell = Config('bcell')
run_preprocess(config_bcell)

# let these run last, they take a while. So in the meantime you can run the workflow
run_deseq(config_bulk)
run_deseq(config_bcell)

