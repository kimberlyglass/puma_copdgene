To re-run the analysis done for this project, there are 3 main steps to follow:
0. Create your folder
1. Run the preprocessing
2. Run the snakemake files to generate Lioness/PUMA files, and the association of each PUMA edge to FEV1/FVC. 
3. Run the notebooks to create figures and tables

0. create a folder "/PUMA" in your folder /proj/regeps/regep00/studies/COPDGene/analyses/*user*/PUMA/
	a. copy in it /src and /data from /proj/regeps/regep00/studies/COPDGene/analyses/remge/PUMA_final/
		rsync -av --exclude='.snakemake' --exclude='__pycache__/' --exclude='modules/__pycache__/' /proj/regeps/regep00/studies/COPDGene/analyses/remge/PUMA_final/src /proj/regeps/regep00/studies/COPDGene/analyses/*user*/PUMA/
		cp -rf /proj/regeps/regep00/studies/COPDGene/analyses/remge/PUMA_final/data /proj/regeps/regep00/studies/COPDGene/analyses/*user*/PUMA/
		
	b. make sure you create a virtual environment similar to mine pypandaenv1 
		module load conda
		conda activate /udd/remge/.conda/envs/new_marghe/
		
	c. copy my snakemake profile into your folder
		cp -rf /udd/remge/.config/snakemake/standard_new /udd/*user*/.config/snakemake/standard_new
	
1. Run the preprocessing:

	a. Log into a cluster node: salloc. 
		If you wan to be more efficient log into the chandl:
		/bin/srun  --mem=28G --mincpus=10 --partition=gpu -t 24:00:00 --pty /bin/bash
	b. then log into the directory with the code:
		cd /proj/regeps/regep00/studies/COPDGene/analyses/*user*/PUMA/src
		make sure that the config.py had the right parameter: folder_name = 'PUMA', user=*user*
	c. load conda and activate env:
	 	# if in GPU node: 
	 		. /modules/tcl/init/bash
			module load conda/22.11.1
		# if not:
			module load conda
		
		conda activate pypandaenv1
	d. run the code:
		python main.py 
		
	output:
	 will create all the dict_input.pickle in the /COPDGene/analyses/*user*/PUMA/output_* folder, containing all the processed input files.
	 
		
	NB.
	Make sure that in the config.py you choose your user name and the folder name 
	If you have any problem with the folders, check the config.py:
        self.path_data = folder with input data
        self.path_lio_puma = temporary PUMA files saved, the one used for the association with FEV1/FVC
	    self.path_output = where the final networks + table and figures are saved        



2. Run workflows:

	a. from cluster headnode, log into the /src/ folder
	b. load conda and activate env.
		module load conda
		conda activate /udd/remge/.conda/envs/new_marghe/
	c. run snakemake files
		snakemake -s workflow.smk  -q --profile=standard_new --jobs=100
		snakemake -s workflow_bcell.smk  -q --profile=standard_new --jobs=100
		snakemake -s workflow_PAX5.smk  -q --profile=standard_new --jobs=100
		
	NB.
	useful logs are saved in /src/.snakemake/log/		

3. Run Jupyter notebooks following this order:
	1. Phenotype.ipynb  
	2. Network Analysis Bulk.ipynb  
	3. Network Analysis BCells.ipynb  
	4. PAX5.ipynb
	
	If needed add the env to your jupyter notebook envs:
	module load conda
	conda activate /udd/remge/.conda/envs/new_marghe/
	python -m ipykernel install --user --name=new_marghe --display-name "Python (env_name)"
	
	
	
	
	
printshellcmds: true
cluster: "sbatch -p linux12h --export=PATH -o ./.snakemake/log/%x.o%j.txt -e ./.snakemake/log/%x.e%j.txt --mem=15G --cpus-per-task=1" #--mem=15G --cpus-per-task=1 
latency-wait: 1200
use-conda: true
conda-prefix: "/proj/relibs/relib00/smk-conda-cache/envs/"
reason: true
notemp: true
jobs: 1








git init
git add ./
git commit -m "first commit"
git branch -M main
#git remote add origin git@changit.bwh.harvard.edu:remge/PUMA_COPDGene.git
git push -u origin main