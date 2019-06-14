from __future__ import division
from IPython.display import Image
import pandas as pd
import numpy as np
import itertools as it
import PyQt4 # and next comment are only executed if we want plots.
import matplotlib
matplotlib.use('qt4agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.linear_model
from scipy.stats.stats import pearsonr
import scipy.stats as st
import random



AGRN = np.transpose(pd.read_csv("Gene_Expression_Atlas/AGRN.tsv",sep= "\t",skiprows= 4))
CHGB = np.transpose(pd.read_csv("Gene_Expression_Atlas/CHGB.tsv",sep= "\t",skiprows= 4))
COL13A1 = np.transpose(pd.read_csv("Gene_Expression_Atlas/COL13A1.tsv",sep= "\t",skiprows= 4))
COL15A1 = np.transpose(pd.read_csv("Gene_Expression_Atlas/COL15A1.tsv",sep= "\t",skiprows= 4))
HSPG2 = np.transpose(pd.read_csv("Gene_Expression_Atlas/HSPG2.tsv",sep= "\t",skiprows= 4))
LAMA2 = np.transpose(pd.read_csv("Gene_Expression_Atlas/LAMA2.tsv",sep= "\t",skiprows= 4))
LAMA5 = np.transpose(pd.read_csv("Gene_Expression_Atlas/LAMA5.tsv",sep= "\t",skiprows= 4))
LAMB2 = np.transpose(pd.read_csv("Gene_Expression_Atlas/LAMB2.tsv",sep= "\t",skiprows= 4))
LOXL3 = np.transpose(pd.read_csv("Gene_Expression_Atlas/LOXL3.tsv",sep= "\t",skiprows= 4))
LRP4 = np.transpose(pd.read_csv("Gene_Expression_Atlas/LRP4.tsv",sep= "\t",skiprows= 4))
PLEC = np.transpose(pd.read_csv("Gene_Expression_Atlas/PLEC.tsv",sep= "\t",skiprows= 4))
TNC = np.transpose(pd.read_csv("Gene_Expression_Atlas/TNC.tsv",sep= "\t",skiprows= 4))
TNXB = np.transpose(pd.read_csv("Gene_Expression_Atlas/TNXB.tsv",sep= "\t",skiprows= 4))
USH2A = np.transpose(pd.read_csv("Gene_Expression_Atlas/USH2A.tsv",sep= "\t",skiprows= 4))
VCAN = np.transpose(pd.read_csv("Gene_Expression_Atlas/VCAN.tsv",sep= "\t",skiprows= 4))

kazoo = [AGRN,CHGB,COL13A1,COL15A1,HSPG2,LAMA2,LAMA5,LAMB2,LOXL3,LRP4,PLEC,TNC,TNXB,USH2A,VCAN]

inputgene = input("Please enter gene ID > ") 
for i in range(15):
    kazoo[i].columns = kazoo[i].iloc[0]
    kazoo[i].index.name = 'newhead'
    kazoo[i].reset_index(inplace = True)
    kazoo[i] = kazoo[i].drop(labels = 0)
    kazoo[i].reset_index()

###Gene value is the index of kazoo vector.
###after creating a plot, you should type "colnames" to know the new value of x.
def getcolnames(geneid):
    gene= geneid
    colnames = list(kazoo[gene])
    return colnames

def plotGeneExpression(geneid):
    gene= geneid
    genenames= ['AGRN','CHGB','COL13A1','COL15A1','HSPG2','LAMA2','LAMA5','LAMB2','LOXL3','LRP4','PLEC','TNC','TNXB','USH2A','VCAN']
    genename = genenames[gene]
    colnames = list(kazoo[gene])
    lengu = len(list(kazoo[gene]))
    x = 'GTEx'
    x2= 'Illumina Body Map'
    kazoo[gene] = kazoo[gene].sort_values(by= x,ascending=False)
    suc = kazoo[gene].loc[:,["newhead",x]].dropna(axis = 0) #for iloc, in THIS VERY CASE, we start from 0.
    suc = suc[(suc != 0).all(1)]
    suc2 = kazoo[gene].loc[:,["newhead",x2]] #for iloc, in THIS VERY CASE, we start from 0.
    suc2 = suc2.reset_index()
    kazumba = suc.loc[:,'newhead']
    kazumba = kazumba.tolist()
    suc2 = suc2.loc[range(len(kazumba)),:]
    print suc
    print suc2
    #suc2 = suc2[(suc2 != 0).all(1)]
    #annotate axis = seaborn axis
    plt.figure(1)
    plt.subplot(1,2,1)
    ax1 = sns.barplot(x="newhead", y= x, data=suc)
    for p in ax1.patches:
                 ax1.annotate("%.2f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()),
                     ha='center', va='center', fontsize=10, color='black', xytext=(0, 20),
                     textcoords='offset points')
    
    ax1.set_ylim([0, 260])
    plt.xticks(fontsize=12,rotation=10)
    plt.yticks(fontsize=13,rotation=0)
    ax1.set_ylabel(x,fontsize=12)
    plt.subplot(1,2,2)
    ax2= sns.barplot(x= 'newhead', y= x2 ,data= suc2)
    for p in ax2.patches:
                 ax2.annotate("%.2f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()),
                     ha='center', va='center', fontsize=10, color='black', xytext=(0, 20),
                     textcoords='offset points')
    
    ax2.set_ylim([0, 260])
    title = 'Tissue Expression for ' + genename + ' (TPM)'
    plt.suptitle(title, fontsize=14)
    plt.xticks(fontsize=12,rotation=10)
    plt.yticks(fontsize=13,rotation=0)
    ax2.set_ylabel(x2,fontsize=12)
    plt.subplots_adjust(bottom=0.15, wspace=0.125)
    plt.show()

plotGeneExpression(inputgene)
    
