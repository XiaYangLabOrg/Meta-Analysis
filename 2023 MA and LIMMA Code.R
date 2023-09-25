##### Microarrary and LIMMA #####
setwd("/Users/zac/Box Sync/Thesis/DEHP/")
## Make sure you have the following in this directory:
#1 old_HUGO_new_HUGO.rda
#2 organ_search_codeZ.r
#3 HUGO_symbols.txt
#4 human_rat_hcop_fifteen_column.txt
#5 human_mouse_hcop_fifteen_column.txt
source("./organ_search_codeZ.r")
load("./old_HUGO_new_HUGO.rda")
#Human genome organzition nomenclature for symbols
EE <- new.env(hash = TRUE)  # equivalent to new.env()
set.seed(123)
list2env(
  setNames(
    as.list(final$newgene), 
    final$oldgene
  ),
  envir = EE
)
EE[["A2MP"]] #hashtable for searching data
### All packages required should be here ###
library(readr)
library(clusterProfiler)
library(Biobase)

library(GEOquery)
library(limma)
library(plyr)
library(dplyr)
library(impute)
library(dplyr)
##############  Adding symbols to genes 
HUGO_symbols <- read_delim("./HUGO_symbols.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)
HUGO_symbols2 <- HUGO_symbols[which(HUGO_symbols$Status == "Approved"),2, drop = FALSE]
RAT_symbols <- read_delim("./human_rat_hcop_fifteen_column.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
RAT_symbols2 <- RAT_symbols[,c(5,12)]
Mouse_symbols <- read_delim("./human_mouse_hcop_fifteen_column.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
Mouse_symbols2 <- Mouse_symbols[,c(5,12)]

################

#### Now we implement the Matrix created below, we are switching it to FA for ease
load("/Users/zac/Box Sync/ToxSci/DEHP/DEHP_Studies_GSE_CtrlGSM_TestGSM.RData")
DEHP <- as.data.frame(DEHP)
### indicating we are storing as characters 
data <- matrix(ncol = 4,nrow = 0)
storage.mode(data) <- "character"
####Pulling our GSE list from matrix and loading the data given the GSE number and RDA
gse_number <- unlist(DEHP$GSE)
for(i in 1:length(gse_number)){
  current_gds <- NA
  try(load(paste0(gse_number[i],".rda")))
  if(is.na(current_gds)){
    dir.create("./test")
    Sys.sleep("1")
    gset <- NULL
    repind = 0
    while(!inherits(gset,'list')) {
      repind <- repind + 1
      if(repind == 25){break}
      tryCatch(
        {
          gset = getGEO(gse_number[i], destdir = paste0(getwd(),"/test"), AnnotGPL = T)
          current_gds <- gset[[1]]
          save(current_gds, file = paste0(gse_number[i],".rda"))
        },
        error=function(e) {Sys.sleep(5); message('retrying'); return(e)}
      )
    }
    unlink("./test",recursive = T)
  }
  
  if(is.na(current_gds)){
    cat("no sample detected",gse_number[i])
    next
  }
  gset <- current_gds
  features <- pData(gset)
  #identify labels
  for(j in 1:ncol(features)){
    features[,j] <- as.character(features[,j])
  }
  #adding sample names for  reference
  Samples <- unlist(DEHP[i,3])
  Controls <- unlist(DEHP[i,2])
  
  ex <- as.data.frame(exprs(gset))
  
  annotation <- fData(gset)
  #checking for RNASeq data, optional
  if(nrow(annotation) == 0){
    if(length(grep("sequencing",as.character(pmid$type[which(pmid$gse %in% gse_number[i])]))) > 0){
      cat("potentially NGS file in ",gse_number[i])
      NGSfilelist <- rbind(NGSfilelist,c(input,gse_number[i]))
      #try(suppfiles <- getGEOSuppFiles(gse_number[i]))
      #filepath <- rownames(suppfiles)
    }
    next
  }
  #checking for NA 
  na_ind <- which(rowMeans(is.na(ex)) > 0.5) 
  if(length(na_ind)/nrow(ex) >= 0.7){
    cat("too much NA in original dataframe ",gse_number[i])
    next
  }
  
  for(column in 1:ncol(annotation)){
    annotation[,column] <- as.character(annotation[,column])
  }
  ENTREZ_SYMBOL <- FALSE
  symbolcolumn <- agrep("Symbol",colnames(annotation),ignore.case = T)

  if(length(symbolcolumn) == 0){
    symbolcolumn <- grep("ENTREZ",colnames(annotation),ignore.case = T)
    if(length(symbolcolumn) > 0){
      ENTREZ_SYMBOL <- TRUE
      cat("ENTREZ ID used in ",gse_number[i])
    }
  }
  if(length(symbolcolumn) == 0){
    orfind <- grep("ORF|SPOT_ID",colnames(annotation),ignore.case = T)
    if(length(orfind) > 0){
      if(length(intersect(annotation[,orfind[1]],HUGO_symbols2$ENTREZID)) > 5000 ||
         length(intersect(annotation[,orfind[1]],Mouse_symbols2$ENTREZID)) > 5000 ||
         length(intersect(annotation[,orfind[1]],RAT_symbols2$ENTREZID)) > 5000
      ){
        symbolcolumn <- orfind
        ENTREZ_SYMBOL <- TRUE
      }else if(length(intersect(annotation[,orfind[1]],HUGO_symbols2$`Approved Symbol`)) > 5000 ||
               length(intersect(annotation[,orfind[1]],Mouse_symbols2$mouse_symbol)) > 5000 ||
               length(intersect(annotation[,orfind[1]],RAT_symbols2$rat_symbol)) > 5000
      ){symbolcolumn <- orfind}
    }
  }
  #more annotations 
  if(length(symbolcolumn) == 0){
    gene_assignmentcol <- grep("gene.assignment", colnames(annotation),ignore.case = T)
    if(length(gene_assignmentcol) > 0){
      cat("complex gene symbol in ",gse_number[i])
      genelist <- gsub(" ","",annotation[,gene_assignmentcol[1]])
      cat("total gene number = ",length(genelist))
      #if(length(genelist) > 500000){
      #  cat("too many genes")
      #  next
      #}
      genelist <- strsplit(genelist,"/|//|///")
      if(length(grep("Homo sapiens",unique(features$organism_ch1),ignore.case = T)) > 0){
        symbols_for_test <- unique(HUGO_symbols2$`Approved Symbol`)
      }else if(length(grep("mus musculus",unique(features$organism_ch1),ignore.case = T)) > 0){
        symbols_for_test <- unique(Mouse_symbols2$mouse_symbol)
      }else{
        symbols_for_test <- unique(RAT_symbols2$rat_symbol)
      }
      currentsymbol <- NULL
      for(gene in 1:length(genelist)){
        currentsymbol[gene] <- paste(symbols_for_test[which(symbols_for_test %in% genelist[[gene]])],collapse = ",")
        if(currentsymbol[gene] == ""){currentsymbol[gene] <- NA}
        if(gene %% 1000 == 0){cat(gene," ")}
      }
      annotation <- data.frame(ID = annotation[,1],Genesymbol = currentsymbol)
      symbolcolumn <- 2
    }
  }
  
  #convert old labels 
  if(length(symbolcolumn) == 0){
    cat("unable to detect symbol column in ",gse_number[i])
    next
  }
  annotation <- annotation[,c(1,symbolcolumn[1])] #here we will ignore multiple mappings
  indold <- which(annotation[,2] %in% final$oldgene)
  annotation[,2] <- as.character(annotation[,2])
  if(length(indold) != 0){
    for( k in 1:length(indold)){
      annotation[indold[k],2] <- EE[[as.character(annotation[indold[k],2])]]
    }
  } 
  #convert old labels to new HUGO symbol via hashtable
  colnames(annotation)[2] <- "Gene symbol"
  annotation$`Gene symbol` <- as.character(annotation$`Gene symbol`)
  ex$ID <- rownames(ex)
  ex <- inner_join(ex,annotation)
  if(length(which(nchar(ex$`Gene symbol`)==0 )) !=0){
    ex <- ex[-which(nchar(ex$`Gene symbol`) == 0),]
  }
  
  
  data_column_ind <- which(colnames(ex) %in% c(Controls,Samples))
  if(length(data_column_ind) == 0){
    cat("unable to match sample and control in ",gse_number[i])
  }
  #choose only sample columns, combine same gene mappings using median
  ex <- aggregate(ex[,data_column_ind], list(Symbol=ex$`Gene symbol`), median, na.rm=T)
  
  if(nrow(ex) <= 1500){
    cat(paste0("no enough gene in ", gse_number[i]))
    next
  }
  #move IDENTIFIER colum
  rownames(ex) <- ex[,"Symbol"]
  ex <- ex[,-1] #delete Identifier
  
  #log2 transformation
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { 
    if(length(which(ex <= 0)) >0){ex[ex <= 0] <- NaN}
    result <- log2(ex) 
  }else{result <- ex
  }
  
  #dealing if there are other unrelated samples to control and drug of interest #something wrong with here
  if(length(c(Samples,Controls)) != nrow(features)){
    print(paste("Some samples were excluded in", gse_number[i]))
    features <- features[which(features$geo_accession %in% c(Samples,Controls)),]
    features <- features[order(features$geo_accession),]
    result <- result[,c(Samples,Controls)]
    result <- result[,order(colnames(result))]
    #if(!identical(as.character(features$geo_accession),colnames(result))){
    #stop("false pruning at", soft_names, input)
    #}
    indsample <- which(features$geo_accession %in% Samples)
    indcontrol <- which(features$geo_accession %in% Controls)
  }

  #To make samples/controls different levels for LIMMA
  indsample <- which(features$geo_accession %in% Samples)
  indcontrol <- which(features$geo_accession %in% Controls)
  agent_factor <- NULL
  agent_factor[indcontrol] <- 1
  agent_factor[indsample] <- 2
  agent_factor <- as.factor(agent_factor)
  if(sum(table(agent_factor) > 1) != 2){
    cat("no enough samples in each group", gse_number[i])
    next
  }
  
  # data processing results
  result$genename <- rownames(result)
  result <- result[,c(ncol(result),1:(ncol(result)-1))]
  
  #result[is.na(result)] <- 0
  NAlabel <- NA
  na_ind <- which(rowSums(is.na(result[,-1])) > 0) 
  if(length(na_ind) > 0){result <- result[-na_ind,]}
  if(length(na_ind)/nrow(result) >= 0.3){
    NAlabel <- "Many NAs here"
  }
  
  #pre-filter with HUGO genes to save time
  ## Make sure to change corresponding chemical where "FA" is 
  annotation <- c("DEHP",gse_number[i],as.character(unique(features$organism_ch1)),as.character(unique(features$characteristics_ch1)),as.character(unique(features$characteristics_ch1.1)),as.character(unique(features$source_name_ch1)),as.character(unique(features$title)))
  annotation2 <- c("DEHP",gse_number[i],Brenda_onto_annot(annotation,processed_BrendaOnto2,processed_BrendaOnto3,finalitems),as.character(unique(features$organism_ch1)))
  
  #trying to match processed data gene names to symbols gene names  (object is called results)
  if(length(grep("Mus musculus",annotation,ignore.case = T)) > 0 ){
    if(ENTREZ_SYMBOL){
      result <- result[which(result$genename %in% Mouse_symbols2$ENTREZID),]
    }else{
      result <- result[which(result$genename %in% Mouse_symbols2$mouse_symbol),]
    }
  }else if(length(grep("Homo sapiens",annotation,ignore.case = T)) > 0){
    if(ENTREZ_SYMBOL){
      result <- result[which(result$genename %in% HUGO_symbols2$ENTREZID),]
    }else{
      result <- result[which(result$genename %in% HUGO_symbols2$`Approved Symbol`),]
    }
  }else if(length(grep("rattus norvegicus",annotation,ignore.case = T)) > 0){
    if(ENTREZ_SYMBOL){
      result <- result[which(result$genename %in% RAT_symbols2$ENTREZID),]
    }else{
      result <- result[which(result$genename %in% RAT_symbols2$rat_symbol),]
    }
  }else{
    cat("No species found for ",gse_number[i])
    next
  }
  #Running ChDir analysis
 #chdir_analysis_example <- chdirAnalysis(result,agent_factor,1,CalculateSig = TRUE)
  #chdir_analysis_example
  #dev.off(dev.list()["RStudioGD"]) #off the graphics 
  
  #test <- as.data.frame(chdir_analysis_example$chdirprops$chdir[[1]])
  #test$`Approved Symbol` <- rownames(test)

  

  #Limma here, not used lets try
  design <- model.matrix(~ agent_factor  + 0)
  colnames(design) <-c("CONT","Test")
  results<-subset(result, select=-c(genename))
  cont.matrix <-  makeContrasts(Diff=Test-CONT, levels = design)
  # Just checking
  cont.matrix
  fit <- lmFit(results, design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit3 <- eBayes(fit2)
  #deg<-topTable(fit3, coef = "Diff", p.value = .05, adjust.method = "fdr", lfc = log2(1.5), number = nrow(results))
  tT <- topTable(fit3, coef = "Diff", adjust.method = "fdr",sort.by = "logFC",number = nrow(results))
  tT <- tT[order(tT$logFC),]

  if(sum(is.na(tT$logFC)) > 0){
    tT <- tT[-which(is.na(tT$logFC)),]
  }
  test<-as.data.frame(tT)
  test[1]<-tT$logFC
  rownames(test)<-rownames(tT)
  test
  ### Correcting alphabetization error from previous mistakes (YW helped)
  tT$genename <- rownames(tT)
  test <- tT[,c("genename","logFC")]
  if(sum(is.na(tT$logFC)) > 0){
    tT <- tT[!is.na(tT$logFC),]
  }  
  #tried changing line down here lets see 
  test <- test[,c(2,1)]
  test
  SUM <- decideTests(fit3)
  summary(SUM)
  setwd("/Users/zac/Box Sync/ToxSci/")
  
  save(test, tT, fit3, file= paste0(GSENAMES[i],"_LIMMA.rda"))
  
  
  #finish final annotation
  colnames(test)[1] <- gse_number[i]
  repind <- grep(gse_number[i],c(colnames(HUGO_symbols2),colnames(RAT_symbols2),colnames(Mouse_symbols2)))
  if(length(repind) > 0){
    colnames(test)[1] <- paste0(gse_number[i],"_",length(repind))
    annotation2[2] <- paste0(gse_number[i],"_",length(repind))
  }
  if(length(grep("Homo sapiens",annotation,ignore.case = T)) > 0 ){
    if(ENTREZ_SYMBOL){colnames(test)[2] <- "ENTREZID"}
    else{colnames(test)[2] <- "Approved Symbol"}
    HUGO_symbols2 <- left_join(HUGO_symbols2,test)
    annotation2[4] <- "Homo sapiens"
  }else if(length(grep("rattus",annotation,ignore.case = T)) > 0 ){
    if(ENTREZ_SYMBOL){colnames(test)[2] <- "ENTREZID"
    }else{colnames(test)[2] <- "rat_symbol"}
    RAT_symbols2 <- left_join(RAT_symbols2, test)
    annotation2[4] <- "rattus norvegicus"
  }else if(length(grep("Mus musculus",annotation,ignore.case = T)) > 0 ){
    if(ENTREZ_SYMBOL){colnames(test)[2] <- "ENTREZID"
    }else{colnames(test)[2] <- "mouse_symbol"}
    Mouse_symbols2 <- left_join(Mouse_symbols2,test)
    annotation2[4] <- "Mus musculus"
  }
  data <- rbind(data,annotation2)
  save(data, HUGO_symbols2, Mouse_symbols2, RAT_symbols2, file = "/Users/zac/Box Sync/ToxSci/DEHP/DEHP_database.rda")
}
#####################################################################################################
load("/Users/zac/Box Sync/ToxSci/DEHP/DEHP_database.rda")
HUGO_symbols3 <- HUGO_symbols2[!(rowMeans(is.na(HUGO_symbols2[,-c(1)])) == 1),]
#IN the case of one for Species!!! must adjust!!!! see how i deleted the ",2" in the above line
HUGO_symbols3 <- HUGO_symbols2[!duplicated(HUGO_symbols2[,1]),]
Mouse_symbols3 <- Mouse_symbols2[!(rowMeans(is.na(Mouse_symbols2[,-c(1,2)])) == 1),]
Mouse_symbols3 <- Mouse_symbols3[!duplicated(Mouse_symbols3),]
Mouse_symbols3 <- Mouse_symbols3[!duplicated(Mouse_symbols3[,1]),]
RAT_symbols3 <- RAT_symbols2[!(rowMeans(is.na(RAT_symbols2[,-c(1,2)])) == 1),]
RAT_symbols3 <- RAT_symbols3[!duplicated(RAT_symbols3),]
RAT_symbols3 <- RAT_symbols3[!duplicated(RAT_symbols3[,1]),]

data <- data[data[,4] %in% c("Homo sapiens","Mus musculus","rattus norvegicus","Rattus norvegicus"),]
#colnames(HUGO_symbols3)[names(HUGO_symbols3)=="Approved Symbol"]<-"human_symbol"   #changing approved to human symbol
data[data[,4] %in% "rattus norvegicus",4] <- "Rattus norvegicus"
Mouse_symbols3<-Mouse_symbols3[,-3]
save(HUGO_symbols3,Mouse_symbols3,RAT_symbols3, data, file = "/Users/zac/Box Sync/DEHP_database.rda_V1.rda")
##############

load("/Users/zac/Box Sync/ToxSci/DEHP/DEHP_database.rda_V1.rda")
library(clusterProfiler)  
library(RobustRankAggreg)
library(gplots)
library(readr)
library(dplyr)
library(enrichR)

names(HUGO_symbols3)[names(HUGO_symbols3) == 'Approved Symbol'] <- 'human_symbol'
HUGO_symbols3<-HUGO_symbols3[,1]

animal_list_up <- list()
animal_list_down <- list()
common_up_frame <- data.frame(Name = character(0),Score = numeric(0),Species = character(0),Organ = character(0))
common_down_frame <- data.frame(Name = character(0),Score = numeric(0),Species = character(0),Organ = character(0))
subdata<-data.frame()
subdata <- data[data[,1] %in% "DEHP2",] 
subdata<- data
View(subdata)
colnames(subdata)<-c("DEHP","currentGsets","currentorgan","species") 
length(subdata[,4])
for(i in 1:length(subdata[,4])){
  currentGsets <- subdata[subdata[,4] %in% species[i],2]
  currentGsets <- unname(currentGsets)
  
  if(length(currentGsets) == 0){
    cat("No ",species[i], "in ",currentGsets)
    next
  }
  if(species[i] %in% "Homo sapiens"){
    subframe <- HUGO_symbols3[,c(T,T,colnames(HUGO_symbols3[,-c(1,2)]) %in% currentGsets)]
  }else if(species[i] %in% "Mus musculus"){
    subframe <- Mouse_symbols3[,c(T,T,colnames(Mouse_symbols3[,-c(1,2)]) %in% currentGsets)]
  }else {
    subframe <- RAT_symbols3[,c(T,T,colnames(RAT_symbols3[,-c(1,2)]) %in% currentGsets)]
  }
  
  library(tibble)
  library(enrichR)
  #rownames(subframe)=subframe[,1]
  names(subframe)[1]<-"human_symbol" #need to ranme to maintain unifrm column names
  subframe<-column_to_rownames(subframe,var = "human_symbol")
  #subframe <- subframe[,-c("human_symbol","mouse_symbol","rat_symbol"),drop = F]
  #cat(currentGsets %in% colnames(subframe))
  
  currentGsets <- colnames(subframe)
  
  if(length(currentGsets) == 0){
    cat("No ",species[i], "in ",currentGsets)
    next
  }
  
  drug_organs <- unique(subdata[subdata[,2] %in% currentGsets,3])
  
  subdata_glist_up <- list()
  subdata_glist_down <- list()
  
  for(j in 1:ncol(subframe)){
    test <- subframe[,j]
    names(test) <- rownames(subframe)
    test <- test[order(test,na.last = NA,decreasing = T)]
    subdata_glist_up[[colnames(subframe)[j]]] <- names(test)
    test <- test[order(test,na.last = NA)]
    subdata_glist_down[[colnames(subframe)[j]]] <- names(test)
  }
  
  organresult_up <- list()
  organresult_down <- list()
  list_for_common_up <- list()
  list_for_common_down <- list()
  
  
  for(j in 1:length(drug_organs)){
    currentorgan <- drug_organs[j]
    Gsetorgans <- intersect(currentGsets,subdata[subdata[,3] %in% currentorgan,2])
    
    aggregateframe  <- aggregateRanks(subdata_glist_up[Gsetorgans], full = T)
    organresult_up[[currentorgan]] <- rownames(aggregateframe)[aggregateframe$Score < 0.05]
    list_for_common_up[[currentorgan]] <- rownames(aggregateframe)
    aggregateframe$Species <- species[i]
    aggregateframe$Organ <- drug_organs[j]
    common_up_frame <- rbind.data.frame(common_up_frame,aggregateframe)
    
    aggregateframe  <- aggregateRanks(subdata_glist_down[Gsetorgans], full = T)
    organresult_down[[currentorgan]] <- rownames(aggregateframe)[aggregateframe$Score < 0.05]
    list_for_common_down[[currentorgan]] <- rownames(aggregateframe)
    aggregateframe$Species <- species[i]
    aggregateframe$Organ <- drug_organs[j]
    common_down_frame <- rbind.data.frame(common_down_frame,aggregateframe)
  }
  
  aggregateframe  <- aggregateRanks(list_for_common_up, full = T)
  organresult_up[["common"]] <- rownames(aggregateframe)[aggregateframe$Score < 0.05]
  list_for_common_up[["common"]] <- rownames(aggregateframe)
  aggregateframe$Species <- species[i]
  aggregateframe$Organ <- "common"
  common_up_frame <- rbind.data.frame(common_up_frame,aggregateframe)
  
  aggregateframe  <- aggregateRanks(list_for_common_down, full = T)
  organresult_down[["common"]] <- rownames(aggregateframe)[aggregateframe$Score < 0.05]
  list_for_common_down[["common"]] <- rownames(aggregateframe)
  aggregateframe$Species <- species[i]
  aggregateframe$Organ <- "common"
  common_down_frame <- rbind.data.frame(common_down_frame,aggregateframe)
  
  animal_list_up[[species[i]]] <- organresult_up
  animal_list_down[[species[i]]] <- organresult_down
  #animal_full_list_up[[species[i]]] <- list_for_common_up
  #animal_full_list_down[[species[i]]] <- list_for_common_down
  
}
allspecies <- names(animal_list_up)

for(j in 1:length(allspecies)){
  currentspecies <-  allspecies[j]
  allorgans <- names(animal_list_up[[currentspecies]])
  for(k in 1:length(allorgans)){
    currentorgan <- allorgans[k]
    overlap <- intersect(animal_list_up[[currentspecies]][[currentorgan]],animal_list_down[[currentspecies]][[currentorgan]])
    if(length(overlap) > 0){
      animal_list_up[[currentspecies]][[currentorgan]] <- setdiff(animal_list_up[[currentspecies]][[currentorgan]],overlap)
      animal_list_down[[currentspecies]][[currentorgan]] <- setdiff(animal_list_down[[currentspecies]][[currentorgan]],overlap)
      common_up_frame <- common_up_frame[-which(common_up_frame$Name %in% overlap & common_up_frame$Species %in% currentspecies & common_up_frame$Organ %in% currentorgan),]
      common_down_frame <- common_down_frame[-which(common_down_frame$Name %in% overlap & common_down_frame$Species %in% currentspecies & common_down_frame$Organ %in% currentorgan),]
      cat(length(overlap)," ")
    }
  }
}

save(animal_list_up,animal_list_down, file = "List of DEHP genes_0.05 cutoff.rda")
#save(animal_list_up,animal_list_down, file = "List of DEHP upregulated genes_0.01 cutoff.rda")



###
