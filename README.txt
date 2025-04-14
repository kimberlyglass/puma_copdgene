To re-run the analysis done for this project, there are 3 main steps to follow:
0. Create your folder
1. Run the preprocessing
2. Run the snakemake files to generate Lioness/PUMA files, and the association of each PUMA edge to FEV1/FVC. 
3. Run the notebooks to create figures and tables

0. create a folder "/PUMA" in your folder /
	
1. Run the preprocessing: python main.py 
	output: will create all the dict_input.pickle in the PUMA/output_* folder, containing all the processed input files.

2. Run workflows:
	a. log into the /src/ folder
	b. load conda and activate env.
	c. run snakemake files
		snakemake -s workflow.smk  -q --profile=standard_new --jobs=100
		snakemake -s workflow_bcell.smk  -q --profile=standard_new --jobs=100
		snakemake -s workflow_PAX5.smk  -q --profile=standard_new --jobs=100

3. Run Jupyter notebooks following this order:
	1. Phenotype.ipynb  
	2. Network Analysis Bulk.ipynb  
	3. Network Analysis BCells.ipynb  
	4. PAX5.ipynb
	
