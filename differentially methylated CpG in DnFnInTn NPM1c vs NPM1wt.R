setwd("C:/Atom/TCGA")
#read the index file
readSets <- function(fname) {
  # read one line at a time
  tmp = readLines(fname)
  # split each line on tab
  tmp2 = sapply(tmp,strsplit,'\t')
  # create a list of the protein IDs by
  # removing the first member of each line
  tmp3 = sapply(tmp2,'[',-1)
  # name the list with the interpro IDs
  names(tmp3) = sapply(tmp2,'[',1)
  # remove any list items with length 0
  # These were blank lines in the original file
  tmp4 = tmp3[sapply(tmp3,length)>0]
  # return result.
  return(tmp4)
}
dat = readSets('index.txt')

#generate list of sample number
Flt3ITD<-as.numeric(dat$Flt3ITD)
Flt3ITD<-as.character(Flt3ITD[!is.na(Flt3ITD)])
NPM1c<-as.numeric(dat$NPM1c)
NPM1c<-as.character(NPM1c[!is.na(NPM1c)])
IDH1<-as.numeric(dat$IDH1)
IDH1<-as.character(IDH1[!is.na(IDH1)])
IDH2<-as.numeric(dat$IDH2)
IDH2<-as.character(IDH2[!is.na(IDH2)])
DNMT3a<-as.numeric(dat$DNMT3a)
DNMT3a<-as.character(DNMT3a[!is.na(DNMT3a)])
MLL<-as.numeric(dat$MLL)
MLL<-as.character(MLL[!is.na(MLL)])
FLT3<-as.numeric(dat$FLT3)
FLT3<-as.character(FLT3[!is.na(FLT3)])
TET2<-as.numeric(dat$TET2)
TET2<-as.character(TET2[!is.na(TET2)])
PMLRAR<-as.numeric(dat$"PML-RAR")
PMLRAR<-as.character(PMLRAR[!is.na(PMLRAR)])

NormCyto<-as.numeric(dat$NormCyto)
NormCyto<-as.character(NormCyto[!is.na(NormCyto)])

Total<-as.numeric(dat$Total)
Total<-as.character(Total[!is.na(Total)])

#jan 08, 2015, Remove Tet2 from analysis
TET2<-intersect(intersect(TET2,Total),NormCyto)
Total<-Total[-match(TET2,Total)]

#tripple subs: Flt3, IDH, DNMT3a for NPM1c patients
#note: Flt3ITD is only considered, Flt3KTD is ignored
FLT<-intersect(intersect(Flt3ITD,Total),NormCyto)
Fn<-intersect(NormCyto,Total[-match(FLT,Total)])
IDH<-intersect(NormCyto,intersect(c(IDH1,IDH2),Total))
In<-intersect(NormCyto,Total[-match(IDH,Total)])
DNMT<-intersect(NormCyto,intersect(Total,DNMT3a))
Dn<-intersect(NormCyto,Total[-match(DNMT,Total)])

DpFpIp<-intersect(intersect(FLT,IDH),DNMT)
DpFnIn<-intersect(DNMT,intersect(Fn,In))
DpFpIn<-intersect(DNMT,intersect(FLT,In))
DpFnIp<-intersect(DNMT,intersect(Fn,IDH))
DnFpIp<-intersect(Dn,intersect(FLT,IDH))
DnFpIn<-intersect(Dn,intersect(FLT,In))
DnFnIp<-intersect(Dn,intersect(Fn,IDH))
DnFnIn<-intersect(Dn,intersect(Fn,In)) # repaired bug on July 05, 2017 changed from intersect(Dn,In)-->intersect(Fn,In)
 
NPM1cList<-list(
  D1F1I1=intersect(NPM1c,DpFpIp),
  D1F1I0=intersect(NPM1c,DpFpIn),
  D1F0I1=intersect(NPM1c,DpFnIp),
  D1F0I0=intersect(NPM1c,DpFnIn),
  D0F1I1=intersect(NPM1c,DnFpIp),
  D0F1I0=intersect(NPM1c,DnFpIn),
  D0F0I1=intersect(NPM1c,DnFnIp),
  D0F0I0=intersect(NPM1c,DnFnIn)
)
NPM1wtList<-list(
  D1F1I1=DpFpIp,
  D1F1I0=DpFpIn[-match(NPM1cList$D1F1I0,DpFpIn)],
  D1F0I1=DpFnIp[-match(NPM1cList$D1F0I1,DpFnIp)],
  D1F0I0=DpFnIn[-match(NPM1cList$D1F0I0,DpFnIn)],
  D0F1I1=DnFpIp[-match(NPM1cList$D0F1I1,DnFpIp)],
  D0F1I0=DnFpIn[-match(NPM1cList$D0F1I0,DnFpIn)],
  D0F0I1=DnFnIp[-match(NPM1cList$D0F0I1,DnFnIp)],
  D0F0I0=DnFnIn[-match(NPM1cList$D0F0I0,DnFnIn)]
)

NPM1c<-NPM1cList$D0F0I0
NPM1wt<-NPM1wtList$D0F0I0

NPM1ctxt<-paste(NPM1cList$D0F0I0,".txt",sep="")
NPM1wttxt<-paste(NPM1wtList$D0F0I0,".txt",sep="")
setwd("C:/Atom/TCGA/DM")


NPM1cdf<-NA
for (i in 1:length(NPM1c)){
  a<-read.table(NPM1ctxt[i])
  b<-a$beta
  NPM1cdf<-data.frame(NPM1cdf,b)
}
rownames(NPM1cdf)<-rownames(a)
NPM1cdf<-NPM1cdf[,-1]
colnames(NPM1cdf)<-NPM1c

NPM1wtdf<-NA
for (i in 1:length(NPM1wt)){
  a<-read.table(NPM1wttxt[i])
  b<-a$beta
  NPM1wtdf<-data.frame(NPM1wtdf,b)
}
rownames(NPM1wtdf)<-rownames(a)
NPM1wtdf<-NPM1wtdf[,-1]
colnames(NPM1wtdf)<-NPM1wt

dfexp<-na.omit(as.matrix(data.frame(NPM1cdf,NPM1wtdf)))
dbeta<-na.omit(rowMeans(NPM1cdf)-rowMeans(NPM1wtdf))

library(minfi)
dmpexp<-dmpFinder(dfexp,pheno=c(rep("NPM1c",length(NPM1c)),rep("NPM1wt",length(NPM1wt))),
                  type="categorical")

dmpexp<-data.frame(dmpexp,dbeta[rownames(dmpexp)])
dmpexp<-na.omit(dmpexp)
colnames(dmpexp)[5]<-"dbeta"
#dbeta<-NA
#for (i in rownames(dmpexp)){
#  dbeta<-c(dbeta,rowMeans(NPM1cdf[i,])-rowMeans(NPM1wtdf[i,]))
#}
#dbeta<-dbeta[-1]

dmpexp.sig<-subset(dmpexp,dmpexp$qval<0.05 & abs(dmpexp$dbeta)>0.25)
dfexp.sig<-dfexp[rownames(dmpexp.sig),]

distance <- dist(as.matrix(t(dfexp.sig)), method="euclidean")
clust <- hclust(distance, method = "average")

distance2 <- dist(as.matrix(dfexp.sig), method="euclidean")
clust2 <- hclust(distance2, method = "average")

png("C:/Atom/TCGA/Results/July 05, 2017/DnFnIpNPM1cvsNPM1wtsig_heatmap_dbeta greater than 0.25 qval smaller 0.05.png",    # create PNG for the heat map        
    width = 200*100,        # 5 x 300 pixels
    height = 150*100,
    res = 1000,            # 300 pixels per inch
    pointsize = 2)        # smaller font size
heatmap(as.matrix(dfexp.sig),col=colorRampPalette(c("darkblue", "yellow"))(256),labRow=NULL,
        margins=c(2.5,1), scale="none", Colv=as.dendrogram(clust), Rowv= as.dendrogram(clust2))
dev.off()
library(Heatplus)
cluster.colname<-colnames(dfexp.sig)[order.dendrogram(as.dendrogram(clust))] # order.dendro gets the exact order of patient smaple on the heatmap - but it's the index of the colnames of the matrix where the heatmap is created from
annNPM1c<-rep("NPM1c",length(NPM1c))
annNPM1wt<-rep("NPM1wt",length(NPM1wt))
genotype<-cbind(c(NPM1c,NPM1wt),c(annNPM1c,annNPM1wt))
rownames(genotype)<-paste("X",c(NPM1c,NPM1wt),sep="")
picketPlot(genotype[cluster.colname,2])
png("C:/Atom/TCGA/Results/July 05, 2017/DnFnIpNPM1cvsNPM1wtsig_heatmap_dbeta greater than 0.25 picketplot.png",    # create PNG for the heat map        
    width = 400*10,        # 5 x 300 pixels
    height = 1000,
    res = 100,            # 300 pixels per inch
    pointsize = 20)        # smaller font size
picketPlot(genotype[cluster.colname,2])
dev.off()

#count how many dmp are found
sum(dmpexp$qval<0.05, na.rm=TRUE)	
#Jan 05, 2015 omit all the NA value

x<-cbind(dfexp[rownames(dmpexp.sig),],dmpexp.sig)
x<-cbind(x,dbeta[rownames(dmpexp.sig)])
cat("",
    colnames(x),
    file="C:/Atom/TCGA/Results/July 05, 2017/SigMethpoints.dat",
    sep="\t",
    append=TRUE)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/SigMethpoints.dat", 
    append=TRUE)
for (i in 1:nrow(x)){
  cat(rownames(x)[i],
      as.numeric(x[i,]),
      file="C:/Atom/TCGA/Results/July 05, 2017/SigMethpoints.dat",
      sep="\t",
      append=TRUE)
  cat("\n", 
      file="C:/Atom/TCGA/Results/July 05, 2017/SigMethpoints.dat", 
      append=TRUE)
}

png("C:/Atom/TCGA/Results/July 05, 2017/volcanoplotofallmethProbes.png",    # create PNG for the heat map        
    width = 250*10,        # 5 x 300 pixels
    height = 200*10,
    res = 500,            # 300 pixels per inch
    pointsize = 5)        # smaller font size
plot(x=dbeta[rownames(dmpexp.sig)],y=(-log10(dmpexp.sig$qval)),xlab="dbeta",ylab="-log10qval")
title(main="DMP in (DNMT3a-,FLT3ITD-,IDH-,TET2 NPM1c) vs (DNMT3a-,FLT3ITD-,IDH-,TET2 NPM1wt), qval<0.05, dbeta>0.25")
dev.off()

cpggrp<-cutree(clust2,k=4)
mpgrp1.hyper<-names(cpggrp[which(cpggrp==2)])
mpgrp2.hypo<-names(cpggrp[which(cpggrp==1)])
mpgrp3.hypo<-names(cpggrp[which(cpggrp==3)])

