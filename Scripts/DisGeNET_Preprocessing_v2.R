###Before using the script please get sure that your current directory is ~/CMS/
### WARNING: DO NOT END THE R SESSION WHEN THIS SCRIPT HAS FINISHED OR YOU WON'T BE ABLE TO USE THE NEXT SCRIPT. 

###DisGeNET database preprocessing

#First, we obtain the number of genes associated to each disease
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

###Then, we filter out disease with just 1 gene
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
  if(sup==1){
    kazi <- kazi[-x]
    kazoo <- kazoo[-x]
    #print(paste(names(kazoo[x]),"eliminated from analysis (gene length=1)", x))
    eliminated <- c(eliminated,kazoo[x])
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
#d <- boxplot(kazi,outline=FALSE) #Uncomment to display boxplot distribution.

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
    kazi <- kazi[-x]
    kazoo <- kazoo[-x]
    #print(paste(names(kazoo[x]),"eliminated from analysis (gene length=1)", x))
    eliminated <- c(eliminated,kazoo[x])
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

##Read Community structure output files from MolTi
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
cinco <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution5MultilayerNetwork5112018.csv")
cinco <- strsplit(cinco,"\t")
seis <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution6MultilayerNetwork5112018.csv")
seis <- strsplit(seis,"\t")
siete <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution7MultilayerNetwork5112018.csv")
siete <- strsplit(siete,"\t")
ocho <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution8MultilayerNetwork5112018.csv")
ocho <- strsplit(ocho,"\t")
nueve <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution9MultilayerNetwork5112018.csv")
nueve <- strsplit(nueve,"\t")
diez <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution10MultilayerNetwork5112018.csv")
diez <- strsplit(diez,"\t")
once <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution11MultilayerNetwork5112018.csv")
once <- strsplit(once,"\t")
doce <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution12MultilayerNetwork5112018.csv")
doce <- strsplit(doce,"\t")
veinte <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution20MultilayerNetwork5112018.csv")
veinte <- strsplit(veinte,"\t")
treinta <- readLines("data/MolTi/Community_Analysis/CommunitiesResolution30MultilayerNetwork5112018.csv")
treinta <- strsplit(treinta,"\t")

edges <- combn(kazoo[[1]],m=2)
edges <- t(edges)

##Next bucle will go elapsed for some time depending on the machine. Please run it in a terminal. 
##In this bucle we obtain a table were all genes that have not been filtered out have their corresponding community in each study. 
hora <- Sys.time()
print(paste("start time=",hora))
diseasecode <- 1
ker <- matrix(0,nrow=0,ncol= 10)
ker[,1] <- disease
ker <- cbind(ker,0)
colnames(ker) <- c("gene","0","0.5","1","1.5","2","2.5","3","3.5","4","disease")
repeat{
  if(diseasecode==5901){
    time <- Sys.time()
    print(paste("ANALISIS ENDED FOR ALL DISEASES. TIME IS",time))
    break()
  }
  old <- Sys.time() #get start time
  disease <- kazoo[[diseasecode]]
  diseasename <- names(kazoo[diseasecode])
  print(paste("for disease",diseasename,"creating table and retrieving community for each gene",old))
  res <- matrix(nrow= length(disease),ncol= 10)
  res[,1] <- disease
  res <- cbind(res,0)
  colnames(res) <- c("gene","0","0.5","1","1.5","2","2.5","3","3.5","4","disease")
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
    # community <- 1
    # repeat{
    #   su <- gen %in% cinco[[community]]
    #   if(su == FALSE){
    #     community <- community + 1
    #     if(community > length(cinco)){
    #       res[i,11] <- NA
    #       break()
    #     }
    #   }
    #   else if(su== TRUE){
    #     time <- Sys.time()
    #     #print(paste(gen,"in five analysis is in community",community," time= ",time))
    #     res[i,11] <- community
    #     break()
    #   }
    # }
    # community <- 1
    # repeat{
    #   su <- gen %in% seis[[community]]
    #   if(su == FALSE){
    #     community <- community + 1
    #     if(community > length(seis)){
    #       res[i,12] <- NA
    #       break()
    #     }
    #   }
    #   else if(su== TRUE){
    #     time <- Sys.time()
    #     #print(paste(gen,"in six analysis is in community",community," time= ",time))
    #     res[i,12] <- community
    #     break()
    #   }
    # }
    # community <- 1
    # repeat{
    #   su <- gen %in% siete[[community]]
    #   if(su == FALSE){
    #     community <- community + 1
    #     if(community > length(siete)){
    #       res[i,13] <- NA
    #       break()
    #     }
    #   }
    #   else if(su== TRUE){
    #     time <- Sys.time()
    #     #print(paste(gen,"in seven analysis is in community",community," time= ",time))
    #     res[i,13] <- community
    #     break()
    #   }
    # }
    # community <- 1
    # repeat{
    #   su <- gen %in% ocho[[community]]
    #   if(su == FALSE){
    #     community <- community + 1
    #     if(community > length(ocho)){
    #       res[i,14] <- NA
    #       break()
    #     }
    #   }
    #   else if(su== TRUE){
    #     time <- Sys.time()
    #     #print(paste(gen,"in eight analysis is in community",community," time= ",time))
    #     res[i,14] <- community
    #     break()
    #   }
    # }
    # community <- 1
    # repeat{
    #   su <- gen %in% nueve[[community]]
    #   if(su == FALSE){
    #     community <- community + 1
    #     if(community > length(nueve)){
    #       res[i,15] <- NA
    #       break()
    #     }
    #   }
    #   else if(su== TRUE){
    #     time <- Sys.time()
    #     #print(paste(gen,"in nine analysis is in community",community," time= ",time))
    #     res[i,15] <- community
    #     break()
    #   }
    # }
    # community <- 1
    # repeat{
    #   su <- gen %in% diez[[community]]
    #   if(su == FALSE){
    #     community <- community + 1
    #     if(community > length(diez)){
    #       res[i,16] <- NA
    #       break()
    #     }
    #   }
    #   else if(su== TRUE){
    #     time <- Sys.time()
    #     #print(paste(gen,"in ten analysis is in community",community," time= ",time))
    #     res[i,16] <- community
    #     break()
    #   }
    # }
    # community <- 1
    # repeat{
    #   su <- gen %in% once[[community]]
    #   if(su == FALSE){
    #     community <- community + 1
    #     if(community > length(once)){
    #       res[i,17] <- NA
    #       break()
    #     }
    #   }
    #   else if(su== TRUE){
    #     time <- Sys.time()
    #     #print(paste(gen,"in eleven analysis is in community",community," time= ",time))
    #     res[i,17] <- community
    #     break()
    #   }
    # }
    # community <- 1
    # repeat{
    #   su <- gen %in% doce[[community]]
    #   if(su == FALSE){
    #     community <- community + 1
    #     if(community > length(doce)){
    #       res[i,18] <- NA
    #       break()
    #     }
    #   }
    #   else if(su== TRUE){
    #     time <- Sys.time()
    #     #print(paste(gen,"in twelve analysis is in community",community," time= ",time))
    #     res[i,18] <- community
    #     break()
    #   }
    # }
    # community <- 1
    # repeat{
    #   su <- gen %in% veinte[[community]]
    #   if(su == FALSE){
    #     community <- community + 1
    #     if(community > length(veinte)){
    #       res[i,19] <- NA
    #       break()
    #     }
    #   }
    #   else if(su== TRUE){
    #     time <- Sys.time()
    #     #print(paste(gen,"in twenty analysis is in community",community," time= ",time))
    #     res[i,19] <- community
    #     break()
    #   }
    # }
    # community <- 1
    # repeat{
    #   su <- gen %in% treinta[[community]]
    #   if(su == FALSE){
    #     community <- community + 1
    #     if(community > length(treinta)){
    #       res[i,20] <- NA
    #       break()
    #     }
    #   }
    #   else if(su== TRUE){
    #     time <- Sys.time()
    #     #print(paste(gen,"in thirty analysis is in community",community," time= ",time))
    #     res[i,20] <- community
    #     break()
    #   }
    # }
  }
  res[,11] <- diseasename
  ker <- rbind(ker,res)
  new= time - old
  print(paste("analysis for",diseasename,"elapsed for",new,"seconds, and ended at",time,"disease",diseasecode," out of 5900"))
  diseasecode <- diseasecode + 1
}
write.table(ker,file= "data/DisGeNET/original/AllDisgenetDiseasesGeneCommunityAssociations.csv",sep= ",",row.names = FALSE)

original <- read.table(file= "data/DisGeNET/original/AllDisgenetDiseasesGeneCommunityAssociations.csv",sep= ",",header=T)

###After obtaining the table we filter out genes that are not included in the multilayer network.
#VECTOR erased is critical, as we will use it in the next script, DisGeNET_Random_Generator.R
x <- 1
lengu <- nrow(original)
erased <- c()
repeat{
  if(x > lengu){
    break()
  }
  if(is.na(original[x,2])==TRUE){
    erased <- c(erased,original[x,1])
    original <- original[-x,]
    #print(paste("row",x,"eliminated"))
    lengu <- nrow(original)
    
  }
  else if(is.na(original[x,2])==FALSE){
    #print(paste("row",x,"is ok"))
    x <- x+1
    lengu <- nrow(original)
  }
}
##Finally we save the table. 
write.table(original,file= "data/DisGeNET/original/AllDisgenetDiseasesGeneCommunityAssociations.csv",sep= ",")


### WARNING: DO NOT END THE R SESSION OR YOU WON'T BE ABLE TO USE THE NEXT SCRIPT. 