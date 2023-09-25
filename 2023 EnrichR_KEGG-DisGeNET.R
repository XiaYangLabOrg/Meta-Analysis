### Enrichment Analysis  ###
setwd(dir = "/Users/zac/Box Sync/ToxSci/")

library(enrichR)  
library(clusterProfiler)
dbs <- c("KEGG_2016", "DisGeNET")
library(UpSetR)
################# TBT ###############
load("List of clustered TBT genes_0.01 cutoff.rda")


## Human In vitro ##
siggenes_down_human<-unlist(animal_list_down$`Homo sapiens`$`liver`)
siggenes_up_human<-unlist(animal_list_up$`Homo sapiens`$`liver`)

humandown<-enrichr(siggenes_down_human,dbs)
humanup<-enrichr(siggenes_up_human,dbs)

# Save
save(humandown,humanup,mixdown,mixup,file = "TBT_human_mix_ENRICHR_results.rda")
## as csv
TBT_human_liver<-c(siggenes_down_human,siggenes_up_human)
#write.csv(TBT_human_liver,file = "TBT_human_liver.csv")

########### PFOA #############
load("List of clustered PFOA genes_0.01 cutoff.rda")
## Mouse is combined of all studies ##
siggenes_down_mouse<-unlist(animal_list_down$`Mus musculus`$`liver`)
siggenes_up_mouse<-unlist(animal_list_up$`Mus musculus`$`liver`)
alllll<-c(siggenes_up_mouse,siggenes_down_mouse)
newpaths<-enrichr(alllll,dbs)
alldown<-enrichr(siggenes_down_mouse,dbs)
allup<-enrichr(siggenes_up_mouse,dbs)

# Save
save(alldown,allup,file = "PFOA_all_ENRICHR_results.rda")

# change according to cluster
PFOA_all_liver<-c(siggenes_down_mouse,siggenes_up_mouse)
#write.csv(PFOA_all_liver,file = "PFOA_all_liver.csv")

############ DEHP ###########
load("List of clustered DEHP genes_0.01 cutoff.rda")
## Rat is combined of all studies ##
siggenes_down_rat<-animal_list_down$`Rattus norvegicus`$liver
siggenes_up_rat<-animal_list_up$`Rattus norvegicus`$liver

alldown<-enrichr(siggenes_down_rat,dbs)
allup<-enrichr(siggenes_up_rat,dbs)

# Save
save(alldown,allup,file = "DEHP_all_ENRICHR_results.rda")
DEHP_all_liver<-c(siggenes_down_rat,siggenes_up_rat)
write.csv(DEHP_all_liver,file = "DEHP_all_liver.csv")

############ BPA ###########
load("List of clustered BPA genes_0.01 cutoff.rda")
##### SPLIT BPA #####
##  Rat only! exlcude GSE130434 and GSE19662 ##
siggenes_down_rat<-animal_list_down$`Rattus norvegicus`$liver
siggenes_up_rat<-animal_list_up$`Rattus norvegicus`$liver

ratdown<-enrichr(siggenes_down_rat,dbs)
ratup<-enrichr(siggenes_up_rat,dbs)

## Human only! exclude GSE69844 ##
siggenes_down_human<-unlist(animal_list_down$`Homo sapiens`$`liver`)
siggenes_up_human<-unlist(animal_list_up$`Homo sapiens`$`liver`)

humandown<-enrichr(siggenes_down_human,dbs)
humanup<-enrichr(siggenes_up_human,dbs)

# Save
save(ratdown,ratup,humandown,humanup, file = "BPA_SPLIT_ENRICHR_results.rda")

BPA_rat_liver<-c(siggenes_down_rat,siggenes_up_rat)
BPA_human_liver<-c(siggenes_down_human,siggenes_up_human)

write.csv(BPA_rat_liver,file = "BPA_rat_liver.csv")
write.csv(BPA_human_liver,file = "BPA_human_liver.csv")


## change according to cluster
BPA_rar_liver<-c(siggenes_down_rat,siggenes_up_rat)
BPA_human_liver<-c(siggenes_down_human,siggenes_up_human)
write.csv(BPA_mix_liver,file = "BPA_mix_liver.csv")

save(ratdown,ratup,humandown,humanup, file ="BPA_enrichr_results.rda")


######################### Not needed, visualization of significant genes overlap UPSET PLOTTTT  #########################  #########################
setwd(dir = "/Users/zac/Box Sync/ToxSci/")

################# TBT ###############
load("List of clustered TBT genes_0.01 cutoff.rda")
## Human In vitro ##
TBT_human_down<-unlist(animal_list_down$`Homo sapiens`$`liver`)
TBT_human_up<-unlist(animal_list_up$`Homo sapiens`$`liver`)

## RAT is mix of two studies (not rat) ##
#TBT_mix_down<-animal_list_down$`Rattus norvegicus`$liver
#TBT_mix_up<-animal_list_up$`Rattus norvegicus`$liver

########### PFOA #############
load("List of clustered PFOA genes_0.01 cutoff.rda")
## Mouse is combined of all studies ##
PFOA_all_down<-unlist(animal_list_down$`Mus musculus`$`liver`)
PFOA_all_up<-unlist(animal_list_up$`Mus musculus`$`liver`)

############ DEHP ###########
load("List of clustered DEHP genes_0.01 cutoff.rda")
## Rat is combined of all studies ##
DEHP_all_down<-animal_list_down$`Rattus norvegicus`$liver
DEHP_all_up<-animal_list_up$`Rattus norvegicus`$liver


############ BPA ###########
load("List of clustered BPA genes_0.01 cutoff.rda")
##### SPLIT BPA #####
##  Rat only! exlcude GSE130434 and GSE19662 ##
BPA_rat_down<-animal_list_down$`Rattus norvegicus`$liver
BPA_rat_up<-animal_list_up$`Rattus norvegicus`$liver

## Human only! exclude GSE69844 ##
BPA_human_down<-unlist(animal_list_down$`Homo sapiens`$`liver`)
BPA_human_up<-unlist(animal_list_up$`Homo sapiens`$`liver`)
library(UpSetR)


### rerutn list of unions

fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  dat<- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  dat[is.na(dat)] <- as.integer(0)
  dat[dat != 0] <- as.integer(1)
  dat <- data.frame(matrix(dat, ncol = length(input), byrow = F))
  dat <- dat[which(rowSums(dat) != 0), ]
  names(dat) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(dat) <- elements
  return(dat)
}

venn.diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
display_venn(ListInput)
##########3
ListInput<-list("BPA_human"= BPA_hum_down,"BPA_rat"=BPA_rat_down, "DEHP_all"= DEHP_all_down,"PFOA_all"= PFOA_all_down, "TBT_human"=TBT_human_down)
upset(fromList(ListInput),order.by = "freq",nsets=6, mainbar.y.label = "Number of Genes", set_size.show = TRUE, keep.order = TRUE,text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.5),set_size.scale_max = 1300)

D<-fromList(ListInput)

Down<-rownames_to_column(D)


ListInput<-list("BPA_human"= BPA_hum_up,"BPA_rat"=BPA_rat_up, "DEHP_all"= DEHP_all_up,"PFOA_all"= PFOA_all_up, "TBT_human"=TBT_human_up)
upset(fromList(ListInput),order.by = "freq",nsets=6, mainbar.y.label = "Number of Genes", set_size.show = TRUE, keep.order = TRUE,text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.5),set_size.scale_max = 1200)


U<-fromList(ListInput)

Up<-rownames_to_column(U)

BPA_human<-c(BPA_hum_up,BPA_hum_down)
BPA_rat<-c(BPA_rat_up,BPA_rat_down)
DEHP_all<-c(DEHP_all_up,DEHP_all_down)
TBT_human<-c(TBT_human_up,TBT_human_down)
PFOA_all<-c(PFOA_all_up,PFOA_all_down)

BPA<-c(BPA_human,BPA_rat)


ListInput<-list("BPA"= BPA, "DEHP_all"= DEHP_all,"PFOA_all"= PFOA_all, "TBT_human"=TBT_human)
upset(fromList(ListInput),order.by = "freq",nsets=6, mainbar.y.label = "Number of Genes", set_size.show = TRUE, keep.order = TRUE,text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.5),set_size.scale_max = 2500)
A<-fromList(ListInput)

