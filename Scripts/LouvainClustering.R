###Before using the script please get sure that your current directory is ~/CMS/
#####Louvain for each layer
rm(list=ls()) #remove current environment

phenotype <- "severe" # "severe" for severe-only analysis, "notsevere" for not-severe genes.

destfile= "Plots/Figure3_" #Directory destination and name of the csv and pdf files that will be retrieved as output.


#libraries needed:
library("AnnotationDbi")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("igraph")
library("gplots")

##First we get the genes in each patient
patients2 <- read.table(file= "data/MolTi/Community_Analysis/originalSevereDiagnosis.csv",sep= ",",header= F)
patients2[,1] <- as.character(patients2[,1])
patients2[,2] <- as.character(patients2[,2])
patients2 <- cbind(patients2,0)
patients2[,3] <- sample(patients2[,2])

#Compound Heterozygous mutations
patients <- read.table(file= "data/InputGenes/compound_all_nov2018.csv")
patients[,1] <- as.character(patients[,1])
patients[,2] <- as.character(patients[,2])
patients <- split(patients[,2],patients[,1])
#severepatients <- patients[which(patients2[,3]=="severe")] #Uncomment and comment the next one for randomizing
severepatients <- patients[patients2[,1][which(patients2[,2]=="severe")]] ##Uncomment and comment the previous one for the original plot
#notseverepatients <- patients[which(patients2[,3]=="notsevere")] #Uncomment and comment the next one for randomizing
notseverepatients <- patients[patients2[,1][which(patients2[,2]=="notsevere")]] #Uncomment and comment the previous one for the original plot
severe<- unlist(severepatients)
names(severe) <- NULL
notsevere <- unlist(notseverepatients)
names(notsevere) <- NULL
severeonly <- setdiff(severe,notsevere)
notsevereonly <- setdiff(notsevere,severe)

severepatientsnames <- names(severepatients)
notseverepatientsnames <- names(notseverepatients)

#CNVs
cnvdata <- read.csv("data/InputGenes/CNV_genes_for_iker.tsv",sep= "\t")
cnvdata[,"gene"] <- as.character(cnvdata[,"gene"])
genis <- c(cnvdata["gene"])$gene
genis <- mapIds(org.Hs.eg.db,keys = genis,column = "SYMBOL",keytype="ENSEMBL",multiVals = "first")
genis["ENSG00000005955"] <- "GGNBP2"
genis["ENSG00000006114"] <- "SYNGR"
genis["ENSG00000108264"] <- "TADA2A"
genis["ENSG00000108270"] <- "AATF"
genis["ENSG00000108272"] <- "DHRS11"
genis["ENSG00000108278"] <- "ZNHI3"
genis["ENSG00000108753"] <- "HNF1B"
genis["ENSG00000129282"] <- "MRM1"
genis["ENSG00000132130"] <- "LHX1"
genis["ENSG00000141140"] <- "MYO19"
genis["ENSG00000141141"] <- "DDX52"
genis["ENSG00000161326"] <- "DUSP14"
genis["ENSG00000167230"] <- "C17orf78"
genis["ENSG00000174093"] <- "RP11-1407O15.2"
genis["ENSG00000184886"] <- "PIGW"
genis["ENSG00000197681"] <- "TBC1D3"
genis["ENSG00000203815"] <- "FAM231D"
genis["ENSG00000219492"] <- "RP11-1396O31.13"
genis["ENSG00000229924"] <- "FAM90A26"
genis["ENSG00000250913"] <- "USP17L23"
genis["ENSG00000268172"] <- "AL590452.1"
genis <- as.matrix(genis)
rownames(genis) <- NULL
genis <- genis[,1]
cnvdata <- cbind(cnvdata,genis)
severedata <- cnvdata[,c("gene","whole_gene",severepatientsnames,"genis")]
milddata <- cnvdata[,c("gene","whole_gene",notseverepatientsnames,"genis")]


i <- 1
coli <- nrow(severedata)
i <- 1
coli <- nrow(severedata)
repeat{
  if(i > coli){
    break()
  }
  vectorinho <- table(unlist(severedata[i,-1]))["2"] #8 times 2 copies
  if(is.na(vectorinho)){
    while(is.na(vectorinho)){
      i <- i+1
      vectorinho <- table(unlist(severedata[i,-1]))["2"]
    }
  }
  if(vectorinho == 8){
    severedata <- severedata[-i,]
  }
  else if(vectorinho != 8){
    i = i+1
  }
  coli <- nrow(severedata)
  #print(i)
}
d <- 1
colo <- nrow(milddata)
repeat{
  if(d > colo){
    break()
  }
  vectorinhod <- table(unlist(milddata[d,-1]))["2"]
  if(is.na(vectorinhod)){
    while(is.na(vectorinhod)){
      d <- d+1
      vectorinhod <- table(unlist(milddata[d,-1]))["2"]
    }
  }
  if(vectorinhod == 12){
    milddata <- milddata[-d,]
  }
  else if(vectorinhod != 12){
    d = d+1
  }
  colo <- nrow(milddata)
  #print(d)
}
rownames(severedata) <- NULL
rownames(milddata) <- NULL
severecnvs <- unique(as.character(severedata[,"genis"]))
mildcnvs <- unique(as.character(milddata[,"genis"]))
cnvsevereonly <- setdiff(severecnvs,mildcnvs)
cnvsnotsevereonly <- setdiff(mildcnvs,severecnvs)

cosa <- read.csv("data/InputGenes/cmsgenes.csv",header=FALSE) #http://www.musclegenetable.fr/4DACTION/Blob_groupe2
#cosa <- read.csv("data/InputGenes/genespaper.csv",header=FALSE) #Table PMID:30552423
cosa <- cosa[,1]
cosa <- as.character(cosa)
if(phenotype == "severe"){
  algenes <- c(severeonly,cosa,cnvsevereonly) #Severe-Only Analysis 
}
if(phenotype == "notsevere"){
  algenes <- c(notsevereonly,cosa,cnvsnotsevereonly) #To get a Mild-Only  Analysis
}
if("KIAA1919" %in% algenes){
  algenes <- replace(algenes,algenes== "KIAA1919","MFSD4B")
}
if("FAM188B" %in% algenes){
  algenes <- replace(algenes,algenes== "FAM188B","MINDY4")
}
if("C14orf159" %in% algenes){
  algenes <- replace(algenes,algenes== "C14orf159","DGLUCY")
}
if("C1orf168" %in% algenes){
  algenes <- replace(algenes,algenes== "C1orf168","FYB2")
}
if("IKBKAP" %in% algenes){
  algenes <- replace(algenes,algenes== "IKBKAP","ELP1")
}
if("SYNGR" %in% algenes){
  algenes <- replace(algenes,algenes== "SYNGR","SYNGR1")
}
if("ZNHI3" %in% algenes){
  algenes <- replace(algenes,algenes== "ZNHI3","ZNHIT3")
}
if("FAM231D" %in% algenes){
  algenes <- replace(algenes,algenes== "FAM231D","LINC01145")
}
if("SELO" %in% algenes){
  algenes <- replace(algenes,algenes== "SELO","SELENOO")
}
if("MFSD7" %in% algenes){
  algenes <- replace(algenes,algenes== "MFSD7","SLC49A3")
}
algenes2 <- algenes
#this lines are commented because we only use them if we want EntrzID.
#genescausales <- mapIds(org.Hs.eg.db,keys = algenes,column = "ENTREZID",keytype="SYMBOL",multiVals = "first")
#genescausales <- as.matrix(genescausales)
#rownames(genescausales) <- NULL
#genescausales <- genescausales[,1]
#genescausales <- unlist(genescausales) #need this unlist as sometimes R does not create a c() vector with the previous line

algenes <- mapIds(org.Hs.eg.db,keys = algenes,column = "ENTREZID",keytype="SYMBOL",multiVals = "first")
names(algenes) <- NULL

### -> Interactome

interactome <- read.table(file= "data/Networks/InteractomaSinDuplciadosJurisica.csv")
interactome <- graph_from_data_frame(interactome[,1:2],directed = FALSE)
community <- cluster_louvain(interactome)
membresia <- membership(community)
### Graph where CNVs, causal genes and compound heterozygous mutations are connected if they are in the same community.
louvaini <- membresia[algenes] 
names(louvaini) <- algenes ###some of the genes are not part of the interactome, so we lose the name. 
louvaini <- split(names(louvaini),louvaini)

### -> Reactome

reactome <- read.table(file= "data/Networks/ReactomeSinDuplicados.csv")
reactome <- graph_from_data_frame(reactome,directed = FALSE)
communityr <- cluster_louvain(reactome)
membresiar <- membership(communityr)
louvainr <- membresiar[algenes] 
names(louvainr) <- algenes ###some of the genes are not part of the interactome, so we lose the name. 
louvainr <- split(names(louvainr),louvainr)

### -> Metabolome

metabolome <- read.table(file= "data/Networks/Recon3DSinDuplicados.csv")
metabolome <- graph_from_data_frame(metabolome,directed = FALSE)
communitym <- cluster_louvain(metabolome)
membresiam <- membership(communitym)
louvainm <- membresiam[algenes] 
names(louvainm) <- algenes ###some of the genes are not part of the interactome, so we lose the name. 
louvainm <- split(names(louvainm),louvainm)

ininteractome <- algenes[which(algenes %in% names(as.list(V(interactome))))]
inmetabolome <- algenes[which(algenes %in% names(as.list(V(metabolome))))]
inreactome <- algenes[which(algenes %in% names(as.list(V(reactome))))]

graphinteractome <- matrix(0,nrow = 0,ncol=2)
graphreactome <- matrix(0,nrow=0,ncol=2)
graphmetabolome <- matrix(0,nrow= 0, ncol=2)

for(i in 1:length(louvaini)){
  louvaini[[i]] <- unique(louvaini[[i]])
  kante <- mapIds(org.Hs.eg.db,keys = louvaini[[i]],column = "SYMBOL",keytype="ENTREZID",multiVals = "first")
  names(kante) <- NULL
  louvaini[[i]] <- kante
  if(length(louvaini[[i]])==1){
    graphinteractome <- rbind(graphinteractome,matrix(louvaini[[i]][1],nrow=1,ncol=2))
  }
  else if(length(louvaini[[i]])>1){
    graphinteractome <- rbind(graphinteractome,t(combn(louvaini[[i]],m=2)))
  }
  print(i)
}

for(r in 1:length(louvainr)){
  louvainr[[r]] <- unique(louvainr[[r]])
  kanter <- mapIds(org.Hs.eg.db,keys = louvainr[[r]],column = "SYMBOL",keytype="ENTREZID",multiVals = "first")
  names(kanter) <- NULL
  louvainr[[r]] <- kanter
  if(length(louvainr[[r]])==1){
    graphreactome <- rbind(graphreactome,matrix(louvainr[[r]][1],nrow=1,ncol=2))
  }
  else if(length(louvainr[[r]])>1){
    graphreactome <- rbind(graphreactome,t(combn(louvainr[[r]],m=2)))
  }
  print(r)
}

for(m in 1:length(louvainm)){
  louvainm[[m]] <- unique(louvainm[[m]])
  kantem <- mapIds(org.Hs.eg.db,keys = louvainm[[m]],column = "SYMBOL",keytype="ENTREZID",multiVals = "first")
  names(kantem) <- NULL
  louvainm[[m]] <- kantem
  if(length(louvainm[[m]])==1){
    graphmetabolome <- rbind(graphmetabolome,matrix(louvainm[[m]][1],nrow=1,ncol=2))
  }
  else if(length(louvainm[[m]])>1){
    graphmetabolome <- rbind(graphmetabolome,t(combn(louvainm[[m]],m=2)))
  }
  print(m)
}

interactomilla <- graph.edgelist(graphinteractome,directed = FALSE)
reactomilla <- graph.edgelist(graphreactome,directed= FALSE)
metabolomilla <- graph.edgelist(graphmetabolome, directed= FALSE)#We first greate a network from the first two columns, which has the list of vertices
#We filter out autointeractions and make the graph undirected
interactomilla <- simplify(interactomilla)
reactomilla <- simplify(reactomilla)
metabolomilla <- simplify(metabolomilla)

write.csv(x= graphinteractome,file= paste0(destfile,"interactome",phenotype,".csv"),row.names=FALSE)
write.csv(x= graphreactome,file= paste0(destfile,"reactome",phenotype,".csv"),row.names=FALSE)
write.csv(x= graphmetabolome,file= paste0(destfile,"metabolome",phenotype,".csv"),row.names=FALSE)


#to obtain graphs for color mapping:
# cosita <- t(combn(cosa,2))
# cnvsita <- t(combn(cnvsevereonly,2))
# cnvsnotsita <- t(combn(cnvsnotsevereonly,2))

# write.csv(x= cosita,file= paste0(destfile,"causalforcolor.csv"),row.names=FALSE)
# write.csv(x= cnvsita,file= paste0(destfile,"severeforcolour.csv"),row.names=FALSE)
# write.csv(x= cnvsnotsita,file= paste0(destfile,"notsevereforcolour.csv"),row.names=FALSE)




#Generate community data table for Edit distance analysis

##
algenes2 <- unique(algenes2)
idataframe <- matrix(data= 0, ncol= length(algenes2),nrow=length(louvaini))
colnames(idataframe) <- algenes2
namesinteractome <- paste0("interactome",names(louvaini)) 
rownames(idataframe)<- namesinteractome

for(i in 1:length(louvaini)){
  kazoo <- colnames(idataframe) %in% louvaini[[i]]
  kazoo <- replace(kazoo, kazoo == TRUE, 1)
  idataframe[i,] <- kazoo
}

mdataframe <- matrix(data= 0, ncol= length(algenes2),nrow=length(louvainm))
colnames(mdataframe) <- algenes2
namesmetabolome <- paste0("metabolome",names(louvainm)) 
rownames(mdataframe)<- namesmetabolome

for(i in 1:length(louvainm)){
  kazoo <- colnames(mdataframe) %in% louvainm[[i]]
  kazoo <- replace(kazoo, kazoo == TRUE, 1)
  mdataframe[i,] <- kazoo
}

rdataframe <- matrix(data= 0, ncol= length(algenes2),nrow=length(louvainr))
colnames(rdataframe) <- algenes2
namesreactome <- paste0("reactome",names(louvainr)) 
rownames(rdataframe)<- namesreactome

for(i in 1:length(louvainr)){
  kazoo <- colnames(rdataframe) %in% louvainr[[i]]
  kazoo <- replace(kazoo, kazoo == TRUE, 1)
  rdataframe[i,] <- kazoo
}

louvaindataframe <- rbind(idataframe,rdataframe,mdataframe)
editdistanceframe <- matrix(0, nrow=nrow(louvaindataframe),ncol=nrow(louvaindataframe))
normalizededitdistanceframe <- matrix(0, nrow=nrow(louvaindataframe),ncol=nrow(louvaindataframe))
rownames(editdistanceframe) <- rownames(louvaindataframe)
colnames(editdistanceframe) <- rownames(editdistanceframe)
rownames(normalizededitdistanceframe) <- rownames(louvaindataframe)
colnames(normalizededitdistanceframe) <- rownames(normalizededitdistanceframe)
mydist <- function(x,y){return(length(which(x+y>1))/length(unique(names(x),names(y))))}
for(i in(1:nrow(louvaindataframe))){
  for(j in (1:nrow(louvaindataframe))){
    sat <- which((louvaindataframe[i,] + louvaindataframe[j,]) != 0)
    reducted <- as.matrix(louvaindataframe[,sat])
    cus <- mydist(x=reducted[i,],y=reducted[j,])
    editdistanceframe[i,j] <- cus
    # cas <- louvaindataframe[i,] + louvaindataframe[j,]
    # cas <- unlist(table(cas != 0)["TRUE"],use.names=FALSE)
    # normalizededitdistanceframe[i,j] <- cus / cas
  }
}

#heatmap

heatmap.2(editdistanceframe,col= cm.colors(512),dendrogram= "column",trace="none")

pdf(file= paste0(destfile,"heatmap",phenotype,".pdf"))
heatmap.2(editdistanceframe,col= cm.colors(128),dendrogram= "column",trace="none",key= FALSE,density.info="none",cexRow=0.7,cexCol=0.7)
dev.off()