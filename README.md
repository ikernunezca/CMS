###### README for CMS Paper Multilayer network analysis scripts. Please, feel free to contact the author for any information and doubts at: iker.nunez@bsc.es

##The 3 network files used in the paper are represented here with the Gene Entrez ID. Internal processing of the files by the scripts change it to the gene name (e.g., "375790" to "AGRN"). 

This Readme file contains the instructions to launch the scripts generated for the different network analysis.

##### Before launching any of the scripts, please take a moment to check if you have installed the library dependencies of these scripts:

	R libraries
		igraph
		AnnotationDbi
		biomaRt
		org.Hs.eg.db
		gplots

	Python libraries
		IPython
		pandas
		numpy
		itertools
		PyQt4
		matplotlib
		seaborn
		sklearn
		scipy
		random

R version used: 3.5.1

Python version used: 2.7.15


#### First -MUST- steps to do before running any script:
    1. Install Jupyter Notebook. Most of the analysis are presented in Jupyter Notebook format in order to simplify reproducibility.
	2. Download the full repository and remember the destination directory where you have saved it: We will work from that very repository in most of the cases.

### A) Monolayer louvain clustering analysis (Supplementary Figures 3 & 4)
	1. Jupyter notebook in /Scripts folder: LouvainClustering_Supplementary_Figure_3_&_4.ipynb
	2. Outputs are:
		2.1 ->  3 edgelists, in .csv format, the 3 community clustering structures, from the phenotype selected, obtained from applying louvain algorithm to the reactome, metabolome and interactome networks. Saved as: /Plots/SupplementaryFigure_3NETWORK.csv
		2.2 -> Supplementary Figure 4 heatmap, depending on the phenotype selected.

### B) Figure 4 modules:
	1. Jupyter notebook in /Scripts folder: Figure_4_Plots.ipynb
	2. Outputs: 
        2.1. CSV file with the edgelist to produce the shared community pertinance graph in Cytoscape (Actual Figure 4)
        2.2. Shared community pertinance graph plotted with the jupyter notebook as output from igraph.


### C) Wilcoxon-Test (Suppl. Figure 5):
	* Splicit files generated during the study's run are available in ~/CMS/data/DisGeNET/random/RandomizationsWith12to30/. Each file correspond to a randomized version of the original downloaded database from DisGeNET.
	1. Jupyter notebook in /Scripts folder: DisGeNET_Preprocessing_WilcoxonTests.ipynb
    2. Outputs: Supplementary Figures 5 A, B & C.


### D) Module significance test:
	1. Jupyter Notebook in /Scripts folder: SignificanceShufflingTest.ipynb
    2. p-values indicate the probability of finding the given number of genes, out of the 15 composing the modules for severe and not-severe phenotypes, upon 1000 iterations of label shuffling.


### E) Gene Expression Atlas Supplementary Figure:
	*Before doing anything, please go to ~/CMS/data/fibroblast_expression/ and extract the zip file in the same directory; the script uses the uncrompessed file. 
	1. Open a terminal/console, and go to the directory (cd) where your CMS folder has been dowloaded.
	2. Within the CMS folder, enter Scripts folder. From this location (~/CMS/Scripts/) launch the python Script 'Suppl_Figure_Gene_Expression_Atlas.py':
		~/CMS/Scripts$ python2.7 Suppl_Figure_Gene_Expression_Atlas.py
	3. After launching the script, we will be requested to type the gene ID. Gene IDs are from each one of the genes in the severe module:
	       ID   Gene		
		0: 'AGRN'
		1: 'CHGB'
		2: 'COL13A1'
		3: 'COL15A1'
		4: 'HSPG2'
		5: 'LAMA2'
		6: 'LAMA5'
		7: 'LAMB2'
		8: 'LOXL3'
		9: 'LRP4'
		10:'PLEC'
		11:'TNC'
		12:'TNXB'
		13:'USH2A'
		14:'VCAN'

	*IMPORTANT: ID 13 (i.e. 'USH2A') will retrieve no result as there is no information in the Gene Expression Atlas output for the gene neither in GTEx nor in Ilumina Body Map (as of 2018)


### F) Fisher test of association between clinical tests and severity (@cirillodavide)
	1. Jupyter Notebook in /Scripts folder: Fisher_CMS.ipynb
	2. Output will be saved at /Plots folder: A barplot with the pvalues of the test. 
