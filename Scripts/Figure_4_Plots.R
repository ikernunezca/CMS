rm(list=ls()) #remove current environment
phenotype <- "notsevere" # "severe" for severe-only analysis, "notsevere" for not-severe genes.

destfile= "Plots/Figure4_" #Directory destination and name of the csv and pdf files that will be retrieved as output.

 
#libraries needed:
library("AnnotationDbi")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("igraph")

#Read Community structure output files from MolTi
cero <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution0MultilayerNetwork5112018.csv")
cero <- strsplit(cero,"\t")
cerocinco <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution05MultilayerNetwork5112018.csv")
cerocinco <- strsplit(cerocinco,"\t")
uno <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution1MultilayerNetwork5112018.csv")
uno <- strsplit(uno,"\t")
unocinco <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution15MultilayerNetwork5112018.csv")
unocinco <- strsplit(unocinco,"\t")
dos <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution2MultilayerNetwork5112018.csv")
dos <- strsplit(dos,"\t")
doscinco <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution25MultilayerNetwork5112018.csv")
doscinco <- strsplit(doscinco,"\t")
tres <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution3MultilayerNetwork5112018.csv")
tres <- strsplit(tres,"\t")
trescinco <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution35MultilayerNetwork5112018.csv")
trescinco <- strsplit(trescinco,"\t")
cuatro <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution4MultilayerNetwork5112018.csv")
cuatro <- strsplit(cuatro,"\t")

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
 
  genescausales <- mapIds(org.Hs.eg.db,keys = algenes,column = "ENTREZID",keytype="SYMBOL",multiVals = "first")
  genescausales <- as.matrix(genescausales)
  rownames(genescausales) <- NULL
  genescausales <- genescausales[,1]
  genescausales <- unlist(genescausales) #need this unlist as sometimes R does not create a c() vector with the previous line
  disease <- genescausales
  res <- matrix(nrow= length(disease),ncol= 10)
  res[,1] <- disease
  res <- cbind(res,0)
  colnames(res) <- c("gene","0","0.5","1","1.5","2","2.5","3","3.5","4","disease") #Create a table with the community ID of each gene in 'genescausales' vector.
  genecommunity <- 0
  for(i in 1:length(disease)){
    gen <- disease[i]
    community <- 1
    cu <- 0
    repeat{
      su <- gen %in% cero[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(cero)){
          res[i,2] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in zero analysis is in community",community," time= ",time))
        res[i,2] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% cerocinco[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(cerocinco)){
          res[i,3] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in zerofive analysis is in community",community," time= ",time))
        res[i,3] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% uno[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(uno)){
          res[i,4] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in one analysis is in community",community," time= ",time))
        res[i,4] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% unocinco[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(unocinco)){
          res[i,5] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in onefive analysis is in community",community," time= ",time))
        res[i,5] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% dos[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(dos)){
          res[i,6] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in two analysis is in community",community," time= ",time))
        res[i,6] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% doscinco[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(doscinco)){
          res[i,7] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in twofive analysis is in community",community," time= ",time))
        res[i,7] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% tres[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(tres)){
          res[i,8] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in three analysis is in community",community," time= ",time))
        res[i,8] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% trescinco[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(trescinco)){
          res[i,9] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in threefive analysis is in community",community," time= ",time))
        res[i,9] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% cuatro[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(cuatro)){
          res[i,10] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in four analysis is in community",community," time= ",time))
        res[i,10] <- community
        break()
      }
    }
  }
  
  ###Delete genes that are not present in the multilayer network 'NA rows'
  x <- 1
  lengu <- nrow(res)
  erased <- c()
  repeat{
    if(x > lengu){
      break()
    }
    if(is.na(res[x,2])==TRUE){
      erased <- c(erased,res[x,1])
      res <- res[-x,]
      #print(paste("row",x,"eliminated"))
      lengu <- nrow(res)
    }
    else if(is.na(res[x,2])==FALSE){
      #print(paste("row",x,"is ok"))
      x <- x+1
      lengu <- nrow(res)
    }
  }
  #Generate the starplot edges: A plot where the edge weight represents the number of times that two genes are in the same community.
  rownames(res) <- res[,1]
  disease <- res[,1]
  edges <- combn(disease,m=2)
  edges <- t(edges)
  edges <- cbind(edges,0,0,0)
  kanoone <- mapIds(org.Hs.eg.db,keys = edges[,1],column = "SYMBOL",keytype="ENTREZID",multiVals = "first")
  kanotwo <- mapIds(org.Hs.eg.db,keys = edges[,2],column = "SYMBOL",keytype="ENTREZID",multiVals = "first")
  kanoone <- as.matrix(kanoone)
  rownames(kanoone) <- NULL
  kanoone <- kanoone[,1]
  kanotwo <- as.matrix(kanotwo)
  rownames(kanotwo) <- NULL
  kanotwo <- kanotwo[,1]
  edges[,4] <- kanoone
  edges[,5] <- kanotwo
  
  x <- 1
  repeat{
    if(x > nrow(edges)){
      break
    }
    a <- edges[x,1]
    b <- edges[x,2]
    sap <- res[a,3:10]
    sbp <- res[b,3:10]
    tuc <- sap==sbp
    kir <- table(tuc)
    kir <- kir["TRUE"]
    kir <- kir[[1]]
    edges[x,3] <- kir
    #print(x)
    x <- x+1
  }
  
  
  ###Get rid of edges that are not connected.
  x <- 1
  repeat{
    if(x > nrow(edges)){
      break()
    }
    if(is.na(edges[x,3])==TRUE){
      edges <- edges[-x,]
      #print(paste("row",x,"eliminated"))
      lengu <- nrow(edges)
    }
    else if(is.na(edges[x,3])==FALSE){
      #print(paste("row",x,"is ok"))
      x <- x+1
      lengu <- nrow(edges)
    }
  }
  
  ###Get only the edges with a value of 8 shared communities.
  x <- 1
  repeat{
    if(x > nrow(edges)){
      break()
    }
    if(edges[x,3] < 8){ 
      edges <- edges[-x,]
      #print(paste("row",x,"eliminated"))
      lengu <- nrow(edges)
    }
    else if(edges[x,3] >=8){
      #print(paste("row",x,"is ok"))
      x <- x+1
      lengu <- nrow(edges)
    }
  }
  g=graph.edgelist(edges[,4:5]) #We first greate a network from the first two columns, which has the list of vertices
  #We filter out autointeractions and make the graph undirected
  g <- simplify(g)
  g <- as.undirected(g) 

pdf(file= paste0(destfile,phenotype,".pdf"))
test.layout <- layout_(g,with_dh(weight.edge.lengths = edge_density(g)/1000))
plot(g, layout = test.layout)
dev.off()

write.csv(x= as_edgelist(g),file= paste0(destfile,phenotype,".csv"),row.names=FALSE) 
