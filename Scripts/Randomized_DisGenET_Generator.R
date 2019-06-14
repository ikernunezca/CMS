###Before using the script please get sure that your current directory is ~/CMS/

#### WARNING: This Script MUST be launched, in the same R session, after DisGeNET_Preprocessing.R, otherwise it will not work correctly. 
#### Here we will generate 100 randomizations of the genes associated to each DisGenNET and analyze them. 
set.seed(123)


numerodeveces <- 0
repeat{
  if(numerodeveces == 1001){
    message("1000 REPETITIONS DONE. RANDOMIZATION ENDED!")
    break()
  }
  

###First we do the same filtering as in DisGeNET_Preprocessing.R
disgenet <- readLines("data/DisGeNET/original/curated_gene_disease_associations.tsv")
disgenet <- strsplit(disgenet,"\t")
genes <- unlist(lapply(disgenet, function(x) x[1]))
diseases <- unlist(lapply(disgenet, function(x) x[3]))
kazoo <- split(genes,diseases)
diseases <- unique(diseases)
genes <- unique(genes)
length(kazoo)
cus <- c()
for(i in 1:length(kazoo)){
  cus <- c(cus, length(kazoo[[i]]))
  #print(i)
}
mean(cus)
lenkazoo <- lapply(kazoo,function(x) length(x))

kazi <- unlist(lenkazoo)
names(kazi) <- NULL
x <- 1
lengu <- length(kazoo)
sup <- 0
eliminated <- c()
repeat{
  if(x > lengu){
    break()
  }
  sup <- kazi[x]
  if(sup ==1){
    eliminated <- c(eliminated,kazoo[x])
    kazi <- kazi[-x]
    kazoo <- kazoo[-x]
    #print(paste(names(kazoo[x]),"eliminated from analysis (gene length=1)", x))
    lengu <- length(kazi)
  }
  else if(sup > 1){
    kazi <- kazi
    kazoo <- kazoo
    #print(paste(names(kazoo[x]),"mantained in analysis (gene length > 1)", x))
    x <- x+1
    lengu <- length(kazi)
  }
}
# d <- boxplot(kazi,outline=FALSE) #Uncomment to display boxplot distribution.
#Now we filter out outlier diseases using boxpot criteria (1.5* Interquartilic Range).
##This turns into filtering out diseases with more than 33 genes from the analysis. 
x <- 1
lengu <- length(kazoo)
sup <- 0
repeat{
  if(x > lengu){
    break()
  }
  sup <- kazi[x]
  if(sup > 33){
    eliminated <- c(eliminated,kazoo[x])
    kazi <- kazi[-x]
    kazoo <- kazoo[-x]
    #print(paste(names(kazoo[x]),"eliminated from analysis (gene length=1)", x))
    lengu <- length(kazi)
  }
  else if(sup <= 33){
    kazi <- kazi
    kazoo <- kazoo
    #print(paste(names(kazoo[x]),"mantained in analysis (gene length > 1)", x))
    x <- x+1
    lengu <- length(kazi)
  }
}
eliminated <- names(eliminated)
#Excluded diseases from the analysis are in vector eliminated.

###Randomized Disegenet Analysis

disgenetcsv <- read.csv("data/DisGeNET/original/curated_gene_disease_associations.tsv",sep= "\t",header = TRUE)
disgenetcsv <- cbind(disgenetcsv,0)
disgenetcsv[,3] <- as.character(disgenetcsv[,3])
prueba <- which(disgenetcsv[,3] %in% names(kazoo))
subset <- disgenetcsv[prueba,]
#### We here take 'erased' vector generated in DisGeNET_Preprocessing.R, to eliminate genes not present in the multilayer network before randomizing everything.
numrow <- 1
lengo <- nrow(subset)
x <- numrow
repeat{
  if(x > lengo){
    break()
  }
  if((subset[x,1] %in% erased)==TRUE){
    #print(paste(subset[x,1],"eliminated"))
    subset <- subset[-x,]
    lengo <- nrow(subset)
  }
  else if((subset[x,1] %in% erased)==FALSE){
    #print(paste(subset[x,1],"nothing to do"))
    x <- x+1
    lengo <- nrow(subset)
  }
}

genes <- subset[,1]
diseases <- subset[,3]
compendio <- split(genes,diseases)
cus <- c()
for(i in 1:length(compendio)){
  cus <- c(cus, length(compendio[[i]]))
  #print(i)
}
mean(cus)
lenkazoo <- lapply(compendio,function(x) length(x))
colnames(subset)[9] <- "randomgene"
subset[,9] <- sample(subset[,1])
randomgenes <- subset[,9]
randomcompendio <- split(randomgenes,diseases) # We probe that each disease has the same genes as in the original database, but randomized.
cus <- c()
for(i in 1:length(randomcompendio)){
  cus <- c(cus, length(randomcompendio[[i]]))
  #print(i)
}
lenrandom <- lapply(randomcompendio,function(x) length(x))
write.table(subset,"data/DisGeNET/random/curated_gene_disease_associations_randomized.tsv",sep= "\t",row.names=FALSE) #we generate a randomized version of this file.


### Now we generate the same table as DisGeNET_Preprocessing.R, but for each randomization it wll be different, so we will be saving a new table in each iteration.

disgenet <- read.table("data/DisGeNET/random/curated_gene_disease_associations_randomized.tsv",sep="\t",header= TRUE)
genes <- disgenet[,9] ###genes random
disgenet[,3] <- as.character(disgenet[,3])
diseases <- disgenet[,3]
kazoo <- split(genes,diseases)
diseases <- unique(diseases)
genes <- unique(genes)
length(kazoo)
cus <- c()
for(i in 1:length(kazoo)){
  cus <- c(cus, length(kazoo[[i]]))
  #print(i)
}
mean(cus)
lenkazoo <- lapply(kazoo,function(x) length(x))

hora <- Sys.time()
print(paste("start time=",hora))
diseasecode <- 1
ker <- matrix(0,nrow=0,ncol= 20)
ker[,1] <- disease
ker <- cbind(ker,0)
colnames(ker) <- c("gene","0","0.5","1","1.5","2","2.5","3","3.5","4","5","6","7","8","9","10","11","12","20","30","disease")
repeat{
  if(diseasecode==5893){ #Here the number of diseases is fewer because we eliminated the ones with no genes in the multilayer before the study, not after it, as in the previous script.
    time <- Sys.time()
    print(paste("RANDOMIZED ANALISIS ENDED FOR ALL DISEASES. TIME IS",time))
    break()
  }
  old <- Sys.time() #get start time
  disease <- kazoo[[diseasecode]]
  diseasename <- names(kazoo[diseasecode])
  #print(paste("for disease",diseasename,"creating table and get community for each gene",old))
  res <- matrix(nrow= length(disease),ncol= 20)
  res[,1] <- disease
  res <- cbind(res,0)
  colnames(res) <- c("gene","0","0.5","1","1.5","2","2.5","3","3.5","4","5","6","7","8","9","10","11","12","20","30","disease")
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
    community <- 1
    repeat{
      su <- gen %in% cinco[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(cinco)){
          res[i,11] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in five analysis is in community",community," time= ",time))
        res[i,11] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% seis[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(seis)){
          res[i,12] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in six analysis is in community",community," time= ",time))
        res[i,12] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% siete[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(siete)){
          res[i,13] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in seven analysis is in community",community," time= ",time))
        res[i,13] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% ocho[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(ocho)){
          res[i,14] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in eight analysis is in community",community," time= ",time))
        res[i,14] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% nueve[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(nueve)){
          res[i,15] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in nine analysis is in community",community," time= ",time))
        res[i,15] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% diez[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(diez)){
          res[i,16] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in ten analysis is in community",community," time= ",time))
        res[i,16] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% once[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(once)){
          res[i,17] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in eleven analysis is in community",community," time= ",time))
        res[i,17] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% doce[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(doce)){
          res[i,18] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in twelve analysis is in community",community," time= ",time))
        res[i,18] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% veinte[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(veinte)){
          res[i,19] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in twenty analysis is in community",community," time= ",time))
        res[i,19] <- community
        break()
      }
    }
    community <- 1
    repeat{
      su <- gen %in% treinta[[community]]
      if(su == FALSE){
        community <- community + 1
        if(community > length(treinta)){
          res[i,20] <- NA
          break()
        }
      }
      else if(su== TRUE){
        time <- Sys.time()
        #print(paste(gen,"in thirty analysis is in community",community," time= ",time))
        res[i,20] <- community
        break()
      }
    }
  }
  res[,21] <- diseasename
  ker <- rbind(ker,res)
  new= time - old
  print(paste("analysis for",diseasename,"elapsed for",new,"seconds, and ended at",time,"disease",diseasecode," out of 5893, randomization", numerodeveces))
  diseasecode <- diseasecode + 1
}
write.table(ker,file= paste0("data/DisGeNET/random/RandomizationsWith12to30/AllDisgenetDiseasesGeneCommunityRandomizedAssociations",numerodeveces,".csv"),sep= ",")
message(paste("repetition",numerodeveces,"done"))
numerodeveces <- numerodeveces + 1
}
