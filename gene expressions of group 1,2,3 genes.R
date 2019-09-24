#observation of CTCF and nCTCF influencced genes (including regulatory regions) in group 1,2,3
#libraries
library(Heatplus)


#extracting gene names
grp1CBS.hypr<-unique(read.table("C:/Atom/TCGA/Results/July 05, 2017/Group1CBSgenes.bed",sep = "\t",stringsAsFactors=F)[,8])
grp2CBS.hypo<-unique(read.table("C:/Atom/TCGA/Results/July 05, 2017/Group2CBSgenes.bed",sep = "\t",stringsAsFactors=F)[,8])
grp3CBS.hypo<-unique(read.table("C:/Atom/TCGA/Results/July 05, 2017/Group3CBSgenes.bed",sep = "\t",stringsAsFactors=F)[,8])
grp1nCBS.hypr<-unique(read.table("C:/Atom/TCGA/Results/July 05, 2017/Group1nCBSgenes.bed",sep = "\t",stringsAsFactors=F)[,8])
grp2nCBS.hypo<-unique(read.table("C:/Atom/TCGA/Results/July 05, 2017/Group2nCBSgenes.bed",sep = "\t",stringsAsFactors=F)[,8])
grp3nCBS.hypo<-unique(read.table("C:/Atom/TCGA/Results/July 05, 2017/Group3nCBSgenes.bed",sep = "\t",stringsAsFactors=F)[,8])


####generating expression profile####
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
DnFnIn<-intersect(Dn,intersect(Fn,In))

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

NPM1ctxt<-paste("Xprsn",paste(NPM1cList$D0F0I0,".txt",sep=""),sep="")
NPM1wttxt<-paste("Xprsn",paste(NPM1wtList$D0F0I0,".txt",sep=""),sep="")

setwd("C:/Atom/TCGA/Xprsn")
NPM1cdf<-NA
for (i in 1:length(NPM1c)){
  a<-read.table(NPM1ctxt[i],header = F, stringsAsFactors = FALSE, skip = 2)
  b<-a[,2]
  NPM1cdf<-data.frame(NPM1cdf,b)
}
rownames(NPM1cdf)<-a[,1]
NPM1cdf<-NPM1cdf[,-1]
colnames(NPM1cdf)<-NPM1c

NPM1wtdf<-NA
NPM1wtdf<-NA
for (i in 1:length(NPM1wt)){
  a<-read.table(NPM1wttxt[i],header = F, stringsAsFactors = FALSE, skip = 2)
  b<-a[,2]
  NPM1wtdf<-data.frame(NPM1wtdf,b)
}
rownames(NPM1wtdf)<-a[,1]
NPM1wtdf<-NPM1wtdf[,-1]
colnames(NPM1wtdf)<-NPM1wt
#########################
#Total gene expression###
#########################
dfexp<-as.matrix(log2(data.frame(NPM1cdf,NPM1wtdf)))### log2 is the official way to do this
sig.pvalue<-NA
fold.change<-NA
for (i in 1:nrow(dfexp)){
  mut<-dfexp[i,1:length(NPM1c)]
  wt<-dfexp[i,(length(NPM1c)+1):ncol(dfexp)]
  sig.pvalue<-c(sig.pvalue,t.test(mut,wt)$p.value)
  fold.change<-c(fold.change,mean(mut)-mean(wt))
}
sig.pvalue<-sig.pvalue[-1]
fold.change<-fold.change[-1]
table.ttest<-cbind(fold.change,sig.pvalue)
rownames(table.ttest)<-rownames(dfexp)
table.sig<-subset(table.ttest,sig.pvalue<0.05)
table.sig<-table.sig[order(table.sig[,1]),]
#table.sig<-subset(table.sig,table.sig[,1]>1)
dfexp.sig<-dfexp[rownames(table.sig),]
pheno.data<-c(rep.int("NPM1c",length(NPM1c)),rep.int("NPM1wt",length(NPM1wt)))

dfexp.sig.grp1CBS<-dfexp.sig[intersect(grp1CBS.hypr,rownames(dfexp.sig)),]
ann1<-annHeatmap(dfexp.sig.grp1CBS,ann=(pheno.data))
png("C:/Atom/TCGA/Results/July 05, 2017/Group1CBS.png",    # create PNG for the heat map        
    width = 9*900,        # 5 x 300 pixels
    height = 9*900,
    res = 1000,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
plot(ann1)
dev.off()

dfexp.sig.grp2CBS<-dfexp.sig[intersect(grp2CBS.hypo,rownames(dfexp.sig)),]
ann1<-annHeatmap(dfexp.sig.grp2CBS,ann=(pheno.data))
png("C:/Atom/TCGA/Results/July 05, 2017/Group2CBS.png",    # create PNG for the heat map        
    width = 9*900,        # 5 x 300 pixels
    height = 9*900,
    res = 1000,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
plot(ann1)
dev.off()

dfexp.sig.grp3CBS<-dfexp.sig[intersect(grp3CBS.hypo,rownames(dfexp.sig)),]
ann1<-annHeatmap(dfexp.sig.grp3CBS,ann=(pheno.data))
png("C:/Atom/TCGA/Results/July 05, 2017/Group3CBS.png",    # create PNG for the heat map        
    width = 9*900,        # 5 x 300 pixels
    height = 9*900,
    res = 1000,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
plot(ann1)
dev.off()

dfexp.sig.grp1nCBS<-dfexp.sig[intersect(grp1nCBS.hypr,rownames(dfexp.sig)),]
ann1<-annHeatmap(dfexp.sig.grp1nCBS,ann=(pheno.data))
png("C:/Atom/TCGA/Results/July 05, 2017/Group1nCBS.png",    # create PNG for the heat map        
    width = 9*900,        # 5 x 300 pixels
    height = 9*900,
    res = 1000,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
plot(ann1)
dev.off()
dfexp.sig.grp2nCBS<-dfexp.sig[intersect(grp2nCBS.hypo,rownames(dfexp.sig)),]
ann1<-annHeatmap(dfexp.sig.grp2nCBS,ann=(pheno.data))
png("C:/Atom/TCGA/Results/July 05, 2017/Group2nCBS.png",    # create PNG for the heat map        
    width = 9*900,        # 5 x 300 pixels
    height = 9*900,
    res = 1000,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
plot(ann1)
dev.off()
dfexp.sig.grp3nCBS<-dfexp.sig[intersect(grp3nCBS.hypo,rownames(dfexp.sig)),]
ann1<-annHeatmap(dfexp.sig.grp3nCBS,ann=(pheno.data))
png("C:/Atom/TCGA/Results/July 05, 2017/Group3nCBS.png",    # create PNG for the heat map        
    width = 9*900,        # 5 x 300 pixels
    height = 9*900,
    res = 1000,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
plot(ann1)
dev.off()

######Whole group1,2,3
dfexp.sig.grp1<-rbind(dfexp.sig.grp1CBS,dfexp.sig.grp1nCBS)
ann1<-annHeatmap(dfexp.sig.grp1,ann=(pheno.data))
png("C:/Atom/TCGA/Results/July 05, 2017/Group1.png",    # create PNG for the heat map        
    width = 9*900,        # 5 x 300 pixels
    height = 9*900,
    res = 1000,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
plot(ann1)
dev.off()

dfexp.sig.grp2<-rbind(dfexp.sig.grp2CBS,dfexp.sig.grp2nCBS)
ann1<-annHeatmap(dfexp.sig.grp2,ann=(pheno.data))
png("C:/Atom/TCGA/Results/July 05, 2017/Group2.png",    # create PNG for the heat map        
    width = 9*900,        # 5 x 300 pixels
    height = 9*900,
    res = 1000,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
plot(ann1)
dev.off()

dfexp.sig.grp3<-rbind(dfexp.sig.grp3CBS,dfexp.sig.grp3nCBS)
ann1<-annHeatmap(dfexp.sig.grp3,ann=(pheno.data))
png("C:/Atom/TCGA/Results/July 05, 2017/Group3.png",    # create PNG for the heat map        
    width = 9*900,        # 5 x 300 pixels
    height = 9*900,
    res = 1000,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
plot(ann1)
dev.off()


cat("Group1.CBS",
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat",
    append=TRUE)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)
cat(grp1CBS.hypr,
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    sep = "\t",
    append=TRUE)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)
cat("Group1.nCBS",
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat",
    append=TRUE)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)
cat(grp1nCBS.hypr,
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    sep = "\t",
    append=TRUE
)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)
cat("Group2.CBS",
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat",
    append=TRUE)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)
cat(grp2CBS.hypo,
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    sep = "\t",
    append=TRUE
)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)
cat("Group2.nCBS",
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat",
    append=TRUE)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)
cat(grp2nCBS.hypo,
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    sep = "\t",
    append=TRUE
)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)

cat("Group3.CBS",
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat",
    append=TRUE)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)
cat(grp3CBS.hypo,
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    sep = "\t",
    append=TRUE
)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)
cat("Group3.nCBS",
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat",
    append=TRUE)
cat("\n", 
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    append=TRUE)
cat(grp3nCBS.hypo,
    file="C:/Atom/TCGA/Results/July 05, 2017/All genes affected by DMP.dat", 
    sep = "\t",
    append=TRUE
)
