###### README for 'Rare disease research workflow using multilayer networks elucidates the molecular determinants of severity in Congenital Myasthenic Syndromes' (https://doi.org/10.1101/2023.01.19.524736) Please, feel free to contact the author for any information and doubts at: iker.nunez@bsc.es

###### PLEASE NOTE THAT THE PAPER IS CURRENTLY UNDER REVISION, some of the scripts and the Cytoscape Session are under active development, with changes to some of the current figures of the manuscript being expected. An update of the revised version of the manuscript in bioRxiv is expected in short-time.

### If the Jupyter Notebooks presented in this repository do not render natively in Github (this is an internal issue from Github), we highly recommend the user to try *nbviewer* (https://nbviewer.org/) to render them.


###### The 3 network files used in the paper are represented here with the Gene Entrez ID. Internal processing of the files by the scripts change it to the gene name (e.g., "375790" to "AGRN"). 

This Readme file contains the instructions to launch the scripts generated for the different network analysis.

Information on the source data provided for reproducibility can be found in the file 'Source_Information_README'

##### Before launching any of the scripts, please take a moment to check if you have installed the library dependencies of these scripts:

	R libraries
		igraph
		AnnotationDbi
		biomaRt
		org.Hs.eg.db
		gplots
		parallelDist
		VennDiagram
		tidyverse
		brainGraph
		reshape2

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

R version used: 3.6.1

Python version used: 2.7.15

##### 
All Source Data files for plotting are provided within the github repository of the project: https://github.com/ikernunezca/CMS
We specifically provide the Cytoscape Session file (‘.cys’) containing all the plots used to produce Figures 3 and 4, as well as Supplementary Figures 3, 6, and 8 in this link: https://github.com/ikernunezca/CMS/blob/master/Cytoscape_Session/CMS_Session.cys

Specific input Source Data Files for creating the Cytoscape Session used to build Figure 3, 4A, 4B, 6 and 8 can be accessed from the following link as csv files: 
https://github.com/ikernunezca/CMS/tree/master/Cytoscape_Session

Additionally, the Cytoscape Session provides an extra plot with the incident interactions considered to build Figure 5, which was manually rendered.

Supplementary Figure 1 source data is provided as Supplementary Table 1. 

Input Data for reproducing Supplementary Figure 2 can be accessed from:
https://github.com/ikernunezca/CMS/tree/master/data/InputGenes
#####


#### First -MUST- steps to do before running any script:
    1. Install Jupyter Notebook. Most of the analysis are presented in Jupyter Notebook format in order to simplify reproducibility.
	2. Download the full repository and remember the destination directory where you have saved it: We will work from that very repository in most of the cases.

 ### A) Cytoscape Session (Rendering Figure 3, 4A and B, Input data for Figure 4C, Supplementary Figures 3, 6 and 8): 
	1. We provide the Cytoscape session file /Cytoscape_Session/CMS_Session.cys with the network plots of the paper (Figure 3, 4 and Supplementary Figures 3, 6 and 8). Additionally, we provide high resolution versions of the same figures within the /Plots folder.

### B) Monolayer louvain clustering analysis (Supplementary Figure 4A and B, Input data for rendering Supplementary Figure 3)
	1. Jupyter notebook in /Scripts folder: LouvainClustering_Supplementary_Figure_3_&_4.ipynb
 	2. Please note that the script must be set for working with severe-specific or non-severe specific mutations prior to be run (Cell 3 of the Jupyter Notebook)
	3. Outputs are:
		3.1 ->  3 edgelists, in .csv format, the 3 community clustering structures, from the phenotype selected, obtained from applying louvain algorithm to the reactome, metabolome and interactome networks. Saved as: /Plots/SupplementaryFigure_3NETWORK.csv
		3.2 -> Supplementary Figure 4 heatmap, depending on the phenotype selected.

### C) Identification of Figure 4 modules (Generation of the input for Cytoscape Session):
	1. Jupyter notebook in /Scripts folder: Figure_4_Plots.ipynb
 	2. Please note that the script must be set for working with severe-specific or non-severe specific mutations prior to be run (Cell 2 of the Jupyter Notebook)
	3. Outputs: 
        3.1. CSV file with the edgelist to produce the shared community pertinance graph in Cytoscape (Actual Figure 4)
        3.2. Shared community pertinance graph plotted with the jupyter notebook as output from igraph.

### D) Fisher test of association between clinical tests and severity (@cirillodavide) (Supplementary Figure 1)
	1. Jupyter Notebook in /Scripts folder: Fisher_CMS.ipynb
	2. Output will be saved at /Plots folder: A barplot with the pvalues of the test. 
	

### E) Wilcoxon-Test (Supplementary Figure 5):
	* Splicit files generated during the study's run are available in ~/CMS/data/DisGeNET/random/RandomizationsWith12to30/. Each file correspond to a randomized version of the original downloaded database from DisGeNET.
	1. Jupyter notebook in /Scripts folder: DisGeNET_Preprocessing_WilcoxonTests.ipynb
    2. Outputs: Supplementary Figures 5 A, B & C.


### F) Module significance test (Supplementary Figure 6) :
	1. Jupyter Notebook in /Scripts folder: SignificanceShufflingTest.ipynb
    2. p-values indicate the probability of finding the given number of genes, out of the 15 composing the modules for severe and not-severe phenotypes, upon 1000 iterations of label shuffling.


### G) Gene Expression Atlas (Supplementary Figure 7):
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

	*IMPORTANT: ID 13 (i.e. 'USH2A') will retrieve no result as there is no information in the Gene Expression Atlas output for the gene neither in GTEx nor in Ilumina Body Map (as of 2018). 
