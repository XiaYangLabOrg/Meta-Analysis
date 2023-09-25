library(readr)
library(clusterProfiler)
library(Biobase)
library(GEOquery)
library(limma)
library(GeoDE)
library(plyr)
library(dplyr)
library(tibble)
library(RobustRankAggreg)

setwd("/Users/zac/Box Sync/GEO Meta Processed Data/Completed_datasets/Genistein/")

## after loading up database for MA change names to join easily.
load("/Users/zac/Box Sync/GEO Meta Processed Data/Completed_datasets/Genistein/Genistein_July27_database.rda_V1.rda")
#change name of symbols for easy binding... dont change this ever
names(HUGO_symbols3)[names(HUGO_symbols3) == 'Approved Symbol'] <- 'human_symbol'

###  with RNAseq thrown in do the following...###
## LOAD Deseq2 results in ##
# change GSE/dose/time acordingly #
load("/Users/zac/Box Sync/GEO Meta Processed Data/Completed_datasets/Genistein/GSE38234_genistein_ECC-1_deseq2results.rda")
rnaseq<-res
rnaseq<-rownames_to_column(rnaseq)
rnaseq$rowname<-toupper(rnaseq$rowname)
rnaseq<-rnaseq[,c(1,3)]
### Dont change human_symbol below (this can stay the same cause mouse and rat frames have human_symbol in them)
colnames(rnaseq) <- c("human_symbol", "GSE38234_ECC1")
HUGO_symbols3<- join(HUGO_symbols3,rnaseq, type="left")

### with 2-color arrays do the following... ###
load("/Users/zac/Box Sync/GEO Meta Processed Data/Completed_datasets/Genistein_TwoColor/GSE19477_LIMMA.rda")

colnames(test) <- c("GSE19477","human_symbol")
### change both HUGO_symbols3 below to whatever species (ie RAT_symbols3 if a rat study)
HUGO_symbols3<-join(HUGO_symbols3,test, type="left")

##if you need to add row of parameters if not dont run these lines (like RNAseq and two-color)
data<-rbind(data,c("Genistein","GSE38234_ECC1","urogenital system","Homo sapiens"))

data<-rbind(data,c("Genistein","GSE38234_T47D","breast","Homo sapiens"))

data<-rbind(data,c("Genistein","GSE19477","urogenital system","Homo sapiens"))


#If you need to edit specific columns in data IF they are wrong! 
#data[6,4]<-c("Homo sapiens")





species <- c("Homo sapiens","Mus musculus","Rattus norvegicus")
animal_list_up <- list()
animal_list_down <- list()
common_up_frame <- data.frame(Name = character(0),Score = numeric(0),Species = character(0),Organ = character(0))
common_down_frame <- data.frame(Name = character(0),Score = numeric(0),Species = character(0),Organ = character(0))
subdata<-data.frame()
subdata <- data[data[,1] %in% "Genistein",] 
subdata<- data
View(subdata)
colnames(subdata)<-c("Genistein","currentGsets","currentorgan","species") 
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
## can do 0.01 by modifying above block with aggregateframe$Score < 0.01. For now just use 0.05 as standard
save(animal_list_up,animal_list_down, file = "List of Genistein genes_0.05 cutoff.rda")
save(data,HUGO_symbols3,Mouse_symbols3,RAT_symbols3, file= "Genistein Meta-signature.rda")
#save(animal_list_up,animal_list_down, file = "List of DEHP upregulated genes_0.01 cutoff.rda")

##### For single studies with LIMMA only
load("/Users/zac/Box Sync/GEO Meta Processed Data/Completed_datasets/HydroxyMethylOctanone/GSE15658_LIMMA.rda")

setwd("/Users/zac/Box Sync/GEO Meta Processed Data/Completed_datasets/Hydroxynonanol/")

results<-decideTests(fit3)
### tells you how many DEGs there are
summary(results)
results<-as.data.frame(results)
UP<- subset(results, Diff==1)
UP<-rownames_to_column(UP)
UP<-UP[,1]
DOWN<- subset(results, Diff==-1)
DOWN<-rownames_to_column(DOWN)
DOWN<-DOWN[,1]

save(UP,DOWN, file="List of Hydroxynonanol genes_0.05 cutoff.rda")
### enrichment
complete <- enrichr(c(UP,DOWN),
                   databases = c("KEGG_2016","GO_Biological_Process_2017","ChEA_2016","ENCODE_TF_ChIP-seq_2015","DisGeNET"))
down <- enrichr(DOWN,
                databases = c("KEGG_2016","GO_Biological_Process_2017","ChEA_2016","ENCODE_TF_ChIP-seq_2015","DisGeNET"))
up <- enrichr(UP,
              databases = c("KEGG_2016","GO_Biological_Process_2017","ChEA_2016","ENCODE_TF_ChIP-seq_2015","DisGeNET"))         

for(j in 1:length(up)){
  currentframe <- up[[j]]
  currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
  if(nrow(currentframe)  == 0){next
  }else{
    if(!exists("drugframe")){
      drugframe <- currentframe
    }else{
      drugframe <- rbind.data.frame(drugframe,currentframe)
    }
  }
}
up_0.05<-drugframe
rm("drugframe")
for(j in 1:length(down)){
  currentframe <- down[[j]]
  currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
  if(nrow(currentframe)  == 0){next
  }else{
    if(!exists("drugframe")){
      drugframe <- currentframe
    }else{
      drugframe <- rbind.data.frame(drugframe,currentframe)
    }
  }
}
down_0.05<-drugframe
rm("drugframe")
for(j in 1:length(complete)){
  currentframe <- complete[[j]]
  currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
  if(nrow(currentframe)  == 0){next
  }else{
    if(!exists("drugframe")){
      drugframe <- currentframe
    }else{
      drugframe <- rbind.data.frame(drugframe,currentframe)
    }
  }
}
complete_0.05<-drugframe

save(complete, down, up, file="Hydroxynonanol_Complete_Pathways.rda")
save(complete_0.05, down_0.05, up_0.05, file="Hydroxynonanol_0.05_Pathways_Combined.rda")


load("List of Genistein genes_0.05 cutoff.rda")
load("Genistein Meta-signature.rda")

### Enrichment analysis

  allspecies <- names(animal_list_up)
 
  for(i in 1:length(allspecies)){
    currentspecies <- allspecies[i]
    allorgans <- names(animal_list_up[[currentspecies]])
    
    for(k in 1:length(allorgans)){
      currentorgan <- allorgans[k]
      res <- enrichr(c(animal_list_up[[currentspecies]][[currentorgan]],animal_list_down[[currentspecies]][[currentorgan]]),
                     databases = c("KEGG_2016","GO_Biological_Process_2017","ChEA_2016","ENCODE_TF_ChIP-seq_2015","DisGeNET"))
      complete <- enrichr(c(animal_list_up[[currentspecies]][[currentorgan]],animal_list_down[[currentspecies]][[currentorgan]]),
                     databases = c("KEGG_2016","GO_Biological_Process_2017","ChEA_2016","ENCODE_TF_ChIP-seq_2015","DisGeNET"))
      down <- enrichr(animal_list_down[[currentspecies]][[currentorgan]],
                        databases = c("KEGG_2016","GO_Biological_Process_2017","ChEA_2016","ENCODE_TF_ChIP-seq_2015","DisGeNET"))
      up <- enrichr(animal_list_up[[currentspecies]][[currentorgan]],
                      databases = c("KEGG_2016","GO_Biological_Process_2017","ChEA_2016","ENCODE_TF_ChIP-seq_2015","DisGeNET"))         
      rm("drugframe")
      for(j in 1:length(res)){
        currentframe <- res[[j]]
        currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
        if(nrow(currentframe)  == 0){next
        }else{
          if(!exists("drugframe")){
            drugframe <- currentframe
          }else{
            drugframe <- rbind.data.frame(drugframe,currentframe)
          }
        }
      }
      
      if(!exists("drugframe")){next}
      drugframe$species <- currentspecies
      drugframe$organ <- currentorgan
      
      if(!exists("finalframe")){
        finalframe <- drugframe
      }else{
        finalframe <- rbind.data.frame(finalframe,drugframe)
      save(complete,down,up, file= paste0("Complete ",currentspecies ," ", currentorgan, " Pathways.rda"))
    
        }}}
  
      save(finalframe, file="0.05_Pathways_Combined.rda")
     
    
  