tab <- read.table('data/Fisher_Input/phenotyping_table.tsv',sep='\t',header=T)

tab <- tab[,-c(1:4)]
colnames(tab) <- 
c(
'Age',
'Sex',
'Severity',
'Pyridostigmine',
'Ptosis',
'Swallowing dysfunction',
'Speech dysfunction',
'Head lift',
'Shoulder lift',
'Leg lift',
'Rise from floor',
'Respiratory dysfunction',
'FCV'
)

tab$Severity <- ifelse(tab$Severity=="severe",'severe','not-severe')
tab$Age <- ifelse(tab$Age<mean(tab$Age),'Y','N')
tab$Sex <- ifelse(tab$Sex=="F",'Y','N')

fisherP <- function(col_num){
	a <- tab[tab[,col_num]!="NI",]
	cms <- matrix(c( 
		length(which(a[a$Severity=="severe",][,col_num]=='Y')),
		length(which(a[a$Severity=="severe",][,col_num]=='N')),
		length(which(a[a$Severity=="not-severe",][,col_num]=='Y')),
		length(which(a[a$Severity=="not-severe",][,col_num]=='N'))),
	nrow = 2,
	dimnames = list(feature = c("Y", "N"),
                	severity = c("Y", "N")))
	p <- fisher.test(cms, alternative = "two.sided")$p.value
	return(p)
}

c <- which(1:ncol(tab)!=which(colnames(tab) == "Severity"))
res <- sapply(c,function(x) fisherP(x))
names(res) <- colnames(tab)[c]
res <- sort(res)

pdf('Plots/Fisher_CMS.pdf',width=10)
op <- par(mar=c(15,4,4,4),ps=18)
bp <- barplot(res,ylab='p-value',names='')
abline(h=0.05,lty='dashed')
text(bp,-0.05,names(res),srt=90,xpd=T,adj=1)
par(op)
dev.off()

print(res)
