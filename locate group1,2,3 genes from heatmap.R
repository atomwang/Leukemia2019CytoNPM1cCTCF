#grouping DMP into 3 groups
###BEDTools was used to intersect mpCBS_bed.bed and hg_19.bed 
#bedtools intersect -a mpCBS_bed.bed -b hg_19.bed -wa -wb > mpCBS_genes.bed
#mpCBS_genes.bed was stored in "C:/Atom/mpCBS_genes.bed"
mpCBS<-as.vector(read.table("C:/Atom/mpCBS/mpCBS.txt")[,1])
mpgrp1.hyper.CBS<-intersect(mpgrp1.hyper,mpCBS)
mpgrp2.hypo.CBS<-intersect(mpgrp2.hypo,mpCBS)
mpgrp3.hypo.CBS<-intersect(mpgrp3.hypo,mpCBS)
mpnCBS.sig<-rownames(dmpnCBS.sig)# from CBS DMP vs nCBS DMP in DnFnInTn
mpgrp1.hyper.nCBS<-intersect(mpgrp1.hyper,mpnCBS.sig)
mpgrp2.hypo.nCBS<-intersect(mpgrp2.hypo,mpnCBS.sig)
mpgrp3.hypo.nCBS<-intersect(mpgrp3.hypo,mpnCBS.sig)
#bed files
location<-read.table("C:/Atom/TCGA/450K_CpGs_hg19_from rajat.bed",stringsAsFactors = FALSE)
rownames(location)<-location$V4
location.Group1<-location[mpgrp1.hyper,]
location.Group1.CBS<-location[mpgrp1.hyper.CBS,]
location.Group1.nCBS<-location[mpgrp1.hyper.nCBS,]
location.Group2<-location[mpgrp2.hypo,]
location.Group2.CBS<-location[mpgrp2.hypo.CBS,]
location.Group2.nCBS<-location[mpgrp2.hypo.nCBS,]
location.Group3<-location[mpgrp3.hypo,]
location.Group3.CBS<-location[mpgrp3.hypo.CBS,]
location.Group3.nCBS<-location[mpgrp3.hypo.nCBS,]

for (i in 1:nrow(location.Group1.CBS)){
  cat(location.Group1.CBS[i,"V1"],location.Group1.CBS[i,"V2"],location.Group1.CBS[i,"V3"],location.Group1.CBS[i,"V4"],file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup1.CBS_bed.bed",sep="\t",append=TRUE)
  cat("\n", file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup1.CBS_bed.bed", append=TRUE)
}
for (i in 1:nrow(location.Group1.nCBS)){
  cat(location.Group1.nCBS[i,"V1"],location.Group1.nCBS[i,"V2"],location.Group1.nCBS[i,"V3"],location.Group1.nCBS[i,"V4"],file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup1.nCBS_bed.bed",sep="\t",append=TRUE)
  cat("\n", file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup1.nCBS_bed.bed", append=TRUE)
}	
for (i in 1:nrow(location.Group2.CBS)){
  cat(location.Group2.CBS[i,"V1"],location.Group2.CBS[i,"V2"],location.Group2.CBS[i,"V3"],location.Group2.CBS[i,"V4"],file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup2.CBS_bed.bed",sep="\t",append=TRUE)
  cat("\n", file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup2.CBS_bed.bed", append=TRUE)
}
for (i in 1:nrow(location.Group2.nCBS)){
  cat(location.Group2.nCBS[i,"V1"],location.Group2.nCBS[i,"V2"],location.Group2.nCBS[i,"V3"],location.Group2.nCBS[i,"V4"],file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup2.nCBS_bed.bed",sep="\t",append=TRUE)
  cat("\n", file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup2.nCBS_bed.bed", append=TRUE)
}	
for (i in 1:nrow(location.Group3.CBS)){
  cat(location.Group3.CBS[i,"V1"],location.Group3.CBS[i,"V2"],location.Group3.CBS[i,"V3"],location.Group3.CBS[i,"V4"],file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup3.CBS_bed.bed",sep="\t",append=TRUE)
  cat("\n", file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup3.CBS_bed.bed", append=TRUE)
}
for (i in 1:nrow(location.Group3.nCBS)){
  cat(location.Group3.nCBS[i,"V1"],location.Group3.nCBS[i,"V2"],location.Group3.nCBS[i,"V3"],location.Group3.nCBS[i,"V4"],file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup3.nCBS_bed.bed",sep="\t",append=TRUE)
  cat("\n", file="C:/Atom/TCGA/Results/July 05, 2017/dmpGroup3.nCBS_bed.bed", append=TRUE)
}	
