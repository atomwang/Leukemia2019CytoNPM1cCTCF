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
DnFnIn<-intersect(Dn,intersect(Fn,In))# repaired bug on July 05, 2017 changed from intersect(Dn,In)-->intersect(Fn,In)

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
