###To compare differentially methylated probes near/within CTCF binding sites to nonCTCF binding site probes
###in DNMT3a-,FLT3-,IDH-TET2- AML in newer coding script
mpCBS<-as.vector(read.table("C:/Atom/mpCBS/mpCBS.txt")[,1])
dmpexp<-read.table("C:/Atom/TCGA/Results/July 05, 2017/SigMethpoints.dat",header = T, sep="\t",stringsAsFactors = FALSE)
rownames(dmpexp)<-dmpexp$X
dmpexp.sig<-subset(dmpexp,dmpexp$qval<0.05 & abs(dmpexp$dbeta)>0.25)
#plot(x=dmpexp$dbeta,y=(-log10(dmpexp$qval)),xlab="beta",ylab="-log10qval",xlim=c(-1,1)) 

dmpCBS.sig<-dmpexp.sig[intersect(dmpexp.sig$X,mpCBS),]
#write.table(dmpCBS.sig, file = "C:/Atom/TCGA/Results/July 05, 2017/dmpCBSDnFnIn.txt", sep = "\t")
png("C:/Atom/TCGA/Results/July 05, 2017/volcanoplotofsigDMPCBS.png",    # create PNG for the heat map        
    width = 250*10,        # 5 x 300 pixels
    height = 200*10,
    res = 500,            # 300 pixels per inch
    pointsize = 5)        # smaller font size
plot(x=dmpCBS.sig$dbeta,y=(-log10(dmpCBS.sig$qval)),xlab="dbeta",ylab="-log10qval",xlim=c(-1,1),ylim = c(0,12))
title(main="DMP near or within CBS (DNMT3a-,FLT3ITD-,IDH-,TET2-) NPM1c vs NPM1wt, qval<0.05, dbeta>0.25")
dev.off()


x<-dmpexp.sig
for (i in rownames(dmpCBS.sig)){
  x<-x[-match(i,rownames(x)),]
}
dmpnCBS.sig<-x
mpnCBS.sig<-rownames(dmpnCBS.sig)

png("C:/Atom/TCGA/Results/July 05, 2017/volcanoplotofsigDMPnCBS.png",    # create PNG for the heat map        
    width = 250*10,        # 5 x 300 pixels
    height = 200*10,
    res = 500,            # 300 pixels per inch
    pointsize = 5)        # smaller font size
plot(x=dmpnCBS.sig$dbeta,y=(-log10(dmpnCBS.sig$qval)),xlab="dbeta",ylab="-log10qval",xlim=c(-1,1),ylim = c(0,12))
title(main="DMP nonCBS (DNMT3a-,FLT3ITD-,IDH-,TET2) NPM1c vs NPM1wt, qval<0.05, dbeta>0.25")
dev.off()
#chisq tst/Fisher
chisq.test(rbind(c(sum(dmpnCBS.sig$dbeta<0, na.rm=TRUE),sum(dmpnCBS.sig$dbeta>0, na.rm=TRUE)),c(sum(dmpCBS.sig$dbeta<0, na.rm=TRUE),sum(dmpCBS.sig$dbeta>0, na.rm=TRUE))))
fisher.test(rbind(c(sum(dmpnCBS.sig$dbeta<0, na.rm=TRUE),sum(dmpnCBS.sig$dbeta>0, na.rm=TRUE)),c(sum(dmpCBS.sig$dbeta<0, na.rm=TRUE),sum(dmpCBS.sig$dbeta>0, na.rm=TRUE))))
