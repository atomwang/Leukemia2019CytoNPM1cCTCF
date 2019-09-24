library(GenomicFeatures)
library(ASpli)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)

## extract GRanges
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
#### if look only at PTPRC set to Chr1 - it's on Chr1####



## subset to get only ranges from chr 1 to 22
auto <- extractSeqlevelsByGroup(species="Homo sapiens", style="UCSC",
                                group="auto")
txby<-txdb
txby <- keepSeqlevels(txby,auto)#, pruning.mode = "coarse")
#seqlevels(txby)

## renaming styles.
newStyle <- mapSeqlevels(seqlevels(txby),"NCBI")
ncbi_gr <- renameSeqlevels(txby, newStyle)
#seqlevels(ncbi_gr)<-"1"
# extract features from annotation
TxDb<-ncbi_gr
features <- binGenome( TxDb )
junctionCoord <- featuresj( features )
geneCoord <- featuresg( features )
binCoord <- featuresb( features )
binMetadata <- mcols( binCoord )


symbols1 <- data.frame( row.names = unique(genes( TxDb )),
                       symbol = paste( 'This is symbol of gene:',
                                       unique(genes( TxDb ))$gene_id )
                      )

features <- binGenome( TxDb, md = symbols1 )


setwd("C:/Atom/TCGA/RNAseqLvl1")
bamFiles<-list.files(pattern = "\\.bam$")
Exon5ass<-c("5788:J003","5788:J035","5788:J018","5788:J019","5788:J004","5788:J020")
Exon5<-NA
gc()
memory.limit(size = 1500000)
for (i in bamFiles){
  bam<-NA
  targets <- data.frame(row.names = c("Reads"
                                       #,"NPM1c_2","NPM1c_3", "NPM1wt_1","NPM1wt_2","NPM1wt_3"
                                        ),
                        bam = i,
                        genotype = c("NPM1"
               #,"Mut", "Mut", "WT","WT","WT"
                                      ) ,
                        stringsAsFactors = FALSE )
  bam<-loadBAM(targets)
  counts <- readCounts (
    features,
    bam,
    targets,
    cores = 1,
    l = 100L,
    maxISize = 50000,
    minAnchor = NULL )
  PTPRC<-countsj(counts)[countsj(counts)$gene=="5788",]
  rownames(PTPRC)<-PTPRC$junction
  Exon5<-cbind(Exon5,PTPRC[Exon5ass,"Reads"])
  bam<-NA
  counts<-NA
  gc()
}

Exon5<-Exon5[,-1]
rownames(Exon5)<-Exon5ass
colnames(Exon5)<-gsub("-.*","\\1",gsub("TCGA-AB-","",bamFiles))


Ex5XR<-Exon5["5788:J035",]/(Exon5["5788:J035",]+0.25*(Exon5["5788:J003",]+Exon5["5788:J019",]+Exon5["5788:J004",]+Exon5["5788:J020",]))

