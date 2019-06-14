###Before using the script please get sure that your current directory is ~/CMS/

### This final script generates Supplementary figure 3. By using Wilcoxon test we test the difference 
## between the communities associated to the genes of each disease from the DisGeNET database, against 100 randomly generated DisGeNET databases.
## In each one of the randomizations we have the same diseases with the same number of genes as the original database, 
## but the genes within each disease are random. 
## The process is performed for each of the resolution parameters separately. 
# Finally we get a boxplot where X is the resolution parameter and Y is the logarithm of each one of the Wilcoxon test p-Values. 
# No value is displayed from resolution parameter as p-Values are closer to zero than what R is capable of calculating.
# PDF with the plot is saved in : plots/boxplotCommunityAnalysisWilconxonTestN100.pdf"

x <- 0
original <- read.table(file= "data/DisGeNET/original/AllDisgenetDiseasesGeneCommunityAssociations.csv",sep= ",",header=T)
original[,"disease"] <- as.character(original[,"disease"])
zero <- c()
zerofunf <- c()
eins <- c()
einsfunf <- c()
zwei <- c()
zweifunf <- c()
drei <- c()
dreifunf <- c()
vier <- c()
funf <- c()
sechs <- c()
sieben <- c()
acht <- c()
neun <- c()
zehn <- c()
elf <- c()
zwolf <- c()
zwanzig <- c()
dreizig <- c()
repeat{
  if(x > 99){
    print(paste("Process ended at",Sys.time()))
    break()
  }
random <- read.table(file= paste0("data/DisGeNET/random/RandomizationsWith12to30/AllDisgenetDiseasesGeneCommunityRandomizedAssociations",x,".csv"),sep= ",",header= T)
random[,"disease"] <- as.character(random[,"disease"])
####Random
lst <- lapply(split(random[,"X0"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X0"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
 freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
zero <- c(zero,res$p.value)

####Random
lst <- lapply(split(random[,"X0.5"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X0.5"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
zerofunf <- c(zerofunf,res$p.value)

####Random
lst <- lapply(split(random[,"X1"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X1"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
eins <- c(eins,res$p.value)

####Random
lst <- lapply(split(random[,"X1.5"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X1.5"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
einsfunf <- c(einsfunf,res$p.value)

####Random
lst <- lapply(split(random[,"X2"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X2"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
zwei <- c(zwei,res$p.value)

####Random
lst <- lapply(split(random[,"X2.5"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X2.5"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
zweifunf <- c(zweifunf,res$p.value)

####Random
lst <- lapply(split(random[,"X3"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X3"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
drei <- c(drei,res$p.value)

####Random
lst <- lapply(split(random[,"X3.5"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X3.5"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
dreifunf <- c(dreifunf,res$p.value)

####Random
lst <- lapply(split(random[,"X4"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X4"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
vier <- c(vier,res$p.value)

####Random
lst <- lapply(split(random[,"X5"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X5"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
funf <- c(funf,res$p.value)

####Random
lst <- lapply(split(random[,"X6"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X6"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
sechs <- c(sechs,res$p.value)

####Random
lst <- lapply(split(random[,"X7"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X7"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
sieben <- c(sieben,res$p.value)

####Random
lst <- lapply(split(random[,"X8"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X8"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
acht <- c(acht,res$p.value)

####Random
lst <- lapply(split(random[,"X9"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X9"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
neun <- c(neun,res$p.value)

####Random
lst <- lapply(split(random[,"X10"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X10"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
zehn <- c(zehn,res$p.value)

####Random
lst <- lapply(split(random[,"X11"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X11"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
elf <- c(elf,res$p.value)

####Random
lst <- lapply(split(random[,"X12"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X12"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
zwolf <- c(zwolf,res$p.value)

####Random
lst <- lapply(split(random[,"X20"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X20"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
zwanzig <- c(zwanzig,res$p.value)

####Random
lst <- lapply(split(random[,"X30"],random[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(random[,"X1.5"])))
freq <- lapply(lst, function(x) x)
randomin <- unlist(freq)

###Original
lst <- lapply(split(original[,"X30"],original[,"disease"]),function(x) table(x))
# freq <- lapply(lst, function(x) length(x[x>1])/length(unique(original[,"X1.5"])))
freq <- lapply(lst, function(x) x)
origin <- unlist(freq)
res <- wilcox.test(randomin,origin)
dreizig <- c(dreizig,res$p.value)

print(paste("Archive",x,"processed.",Sys.time()))
x <- x+1
}

pdf(file= "Plots/boxplotCommunityAnalysisWilconxonTestN100.pdf")
boxplot(log(zero),log(zerofunf),log(eins),log(einsfunf),log(zwei),log(zweifunf),log(drei),log(dreifunf),log(vier),log(funf),log(sechs),log(sieben),log(acht),log(neun),log(zehn),log(elf),log(zwolf),log(zwanzig),log(dreizig),ylab= "log (Wilcoxon test p-values)",main= "Boxplot log (p-values) Original vs randomized community analysis",xaxt="n",xlab= "Resolution Parameter")
#boxplot(log(zero),log(zerofunf),log(eins),log(einsfunf),log(zwei),log(zweifunf),log(drei),log(dreifunf),log(vier),ylab= "log (Wilcoxon test p-values)",main= "Boxplot log (p-values) Original vs randomized community analysis",xaxt="n",xlab= "Resolution Parameter")
axis(side= 1,at= 1:19,labels= c(0,0.5,1,1.5,2,2.5,3,3.5,4,5,6,7,8,9,10,11,12,20,30))
#axis(side= 1,at= 1:9,labels= c(0,0.5,1,1.5,2,2.5,3,3.5,4))
dev.off()
