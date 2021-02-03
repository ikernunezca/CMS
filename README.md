###### README for CMS Paper Multilayer network analysis scripts. Please, feel free to contact the author for any information and doubts at: iker.nunez@bsc.es

##The 3 network files used in the paper are represented here with the Gene Entrez ID. Internal processing of the files by the scripts change it to the gene name (e.g., "375790" to "AGRN"). 

This Readme file contains the instructions to launch the scripts generated for the different network analysis. As some of them can take some time to end, we highly recommend the user to launch everything from a terminal rather than other programs such as RStudio to improve the performance.

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


####First -MUST- steps to do before running any script:
	1. Download the full repository and remember the destination directory where you have saved it: We will work from that very repository in most of the cases.
	2. Go to the ~/CMS/data/Networks directory and uncompress "InteractomaSinDuplciadosJurisica.zip" file.




### A) Monolayer louvain clustering analysis (Suppl. Figure 4, Figure 3
	1. Open a terminal/console, and go to the directory (cd) where your CMS folder has been downloaded.
	2. Access CMS folder and launch R from the terminal.
	3. Open 'LouvainClustering.R' in a text editor and change line 3 as desired to set the plot generation between severe and not severe. Copy all content from the script and run it in the terminal you opened.
	4. Outputs will be saved in /Plots folder, in your CMS folder.
	5. Outputs are:
		-5.1 ->  3 edgelists, in .csv format, the 3 community clustering structures, from the phenotype selected, obtained from applying louvain algorithm to the reactome, metabolome and interactome networks.
		-5.2 -> A pdf, with the heatmap of figure 3, depending on the phenotype selected.

### B) Figure 4 modules:
	1. Open a terminal/console, and go to the directory (cd) where your CMS folder has been downloaded.
	2. Access CMS folder and launch R from the terminal.
	3. Open 'Figure_4_Plots.R' in a text editor and change line 2 as desired to set the plot generation between severe and not severe.
	3. Copy all content from Figure_4_Plots.R and paste it into the R session. 
	4. Output files (pdf graph and csv edgelist) will be in /Plots folder.


### C) Wilcoxon-Test (Suppl. Fig 2) (Suppl. Figure 3):
	* Splicit files generated during the study's run are available in ~/CMS/data/DisGeNET/random/RandomizationsWith12to30/. Each archive correspond to a randomized version of the original downloaded database from DisGeNET. To run the process with this files, just jump to the point 5 of this paragraph.

	1. Open a terminal/console, and go to the directory (cd) where your CMS folder has been dowloaded.
	2. Access CMS folder and launch R from the terminal.
	3. Launch the 'DisGeNET_Preprocessing_v2.R' script's content in this terminal. DO NOT close the R session when the run is over.
	4. Once 'DisGeNET_Preprocessing_v2.R' script has ended, without exiting the R session, copy all content from 'Randomized_DisGenET_Generator_v2.R' script and paste it into the R session. This step may take quite some time. In local, 1 module Intel® Xeon(R) E-2124 CPU, process lasted for 26 h. 
	5. Once 'Randomized_DisGenET_Generator_v2.R' script is over, copy all content from 'WilcoxonTest_DisGeNET_Randomizations.R' script and paste it into the R session. After running it, the output plot -Suppl Figure 3- will be saved in /Plots folder.
	
	*. To run the same analysis, from resolution parameter = 0 to 30, use the other versions of the scripts (the ones that have not the structure "XXXX_v2.R").


### D) Module significance test:
	1. Open a terminal/console, and go to the directory (cd) where your CMS folder has been dowloaded.
	2. Access CMS folder and launch R from the terminal.
	3. Before launching 'SignificanceShufflingTest.R', open the script in a text editor and follow the instructions from line 5 to select the desired analysis (between 'severe' and 'notsevere'). Copy and paste the script's content into te R session.
	4. When the process is finished, you can retrieve the output information from the following variables in the R session:
		-largestcomponents: A list with all the sizes of the largest components from each iteration.
		-amigosenlargestcomponents: A list with the number of genes from the module of interest that are included in the largest component of each iteration.
		-amigoscausalesenlargestcomponents: A list with the number of causal genes from the module of interest that are included in the largest component of each iteration.


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

	*IMPORTANT: ID 13 (i.e. 'USH2A') will retrieve no result as there is no information in the Gene Expression Atlas output for the gene neither in GTEx nor in Ilumina Body Map.


### F) Fisher test of association between clinical tests and severity (@cirillodavide)
	1. Run the script for fisher test: 
		~/CMS$ Rscript Scripts/Fisher_CMS.R
	2. Output will be saved at /Plots folder: A barplot with the pvalues of the test. 
