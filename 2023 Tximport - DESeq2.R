#### TxImport and DESeq2 Analysis #####
## Make sure quant.sf files are moved from scratch to another directory like the one below
setwd("/Users/zac/Box Sync/Thesis/RNAseq_Salmon/")
library(tximport)
library(DESeq2)

#####  Download gtf from ensembl ####
gtf <- rtracklayer::import('./Mus_musculus.GRCm38.101.gtf')
gtf_df=as.data.frame(gtf)

tx2gene <- gtf_df[,c("transcript_id","gene_name")]
tx2gene <- unique(tx2gene)

tx2gene <- tx2gene[!is.na(tx2gene$transcript_id),]
tx2gene <- tx2gene[!is.na(tx2gene$gene_name),]
tx2gene$gene_id <- gsub("\\.[0-9]","",tx2gene$gene_name)
## Saving, we do not need to keep making this after the first time
save(tx2gene,file = "./Common_Mousegenecode2gene.rda")

#### If tx2gene file was previously made, start from here ####
load("Common_Mousegenecode2gene.rda")
files<-file.path("GSE119441", list.files("GSE119441"),"quant.sf")
names(files)<-list.files("GSE119441")
txi <- tximport(files, type="salmon", tx2gene=tx2gene,ignoreAfterBar = T,ignoreTxVersion = T)
data<-txi$counts%>%round()%>%tibble() 

### save expression data ###
save(txi,file="GSE119441_expression.rda")
## In order to make use of metadata, we had to download meta via GEO, then load this 
load("TBT_GSE119441.rda")

####  PLEASE READ ALL TEXT ######
## Utilizing metadata to establish levels

coldata<- as.data.frame(colnames(txi$counts))
coldata[,2]<-as.factor(result$original$`treatment:ch1`)
## If you cant pull because of design from meta, we can manually annotate in excel ##

### Since this example is not neat, we will manually annotate in excel by creating file "conditions.csv" and include sample and treatment as the two column names ###

##  Manual annotation using excel
coldata<-read.csv("conditions.csv") 
colnames(coldata) <- c("sample", "treatment")
## running deseq2 ##
dds <- DESeqDataSetFromTximport(txi,colData=coldata,design = ~treatment )
keep <- rowMeans(counts(dds) >= 1) >= 0.5
dds <- dds[keep,]
dds<-DESeq(dds)
res <- as.data.frame(results(dds, contrast=c("treatment","chemical","control")))
summary(res)
# Order by p-value
res <- res[order(res$pvalue),]
# Saving as csv
write.csv(res, file="GSE119441_deseq2results.csv")

