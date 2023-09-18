## Download the necessary packages if you do not already have them!!!
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(GEOquery)

#sample identification funcion
sampleidentification <- function(input, gse_number){
  current_gds <- NA
  try(load(paste0(gse_number,".rda")))
  
  if(is.na(current_gds)){
    dir.create("./test")
    Sys.sleep("1")
    gset <- NULL
    repind = 0
    #repeat 25 times on connection
    while(!inherits(gset,'list')) {
      repind <- repind + 1
      if(repind == 25){break}
      tryCatch(
        {
          gset = getGEO(gse_number, destdir = paste0(getwd(),"/test"), AnnotGPL = T)
          current_gds <- gset[[1]]
          save(current_gds, file = paste0(gse_number,".rda"))
        },
        error=function(e) {Sys.sleep(5); message('retrying'); return(e)}
      )
    }
    unlink("./test",recursive = T)
  }
  if(is.na(current_gds)){
    cat("no sample detected",gse_number)
    return(0)
  }
  features <- pData(current_gds)
  if(length(grep("Homo sapiens|rattus|mus musculus",unique(features$organism_ch1),ignore.case = T)) == 0 ){
    cat("no species for",gse_number)
    return(0)
  }
  
  #identify labels
  for(j in 1:ncol(features)){
    features[,j] <- as.character(features[,j])
  }
  #agilent detection
  Agilent_label <- F
  if(length(grep("ch2",colnames(features))) > 0){
    cat("two color array detected",gse_number)
    Agilent_label <- T
    sub <- features[,grep("ch2",colnames(features))]
    if(sum(duplicated(sub)) == (nrow(sub)-1)){
      cat("Agilent array with common referance")
      Agilent_label <- F
    }
  }
  
  if(Agilent_label){
    cat("Agilent complex array in ",gse_number)
    return(0)
  }
  
  #adding sample names for  referance
  
  #for acolumns could contain treatment and control, we use source and character plus title column to search for treatment vs control
  indagent <- unique(c(grep("source",colnames(features)), grep("characteristic",colnames(features)),grep("treatment",colnames(features)),grep("description",colnames(features))))
  indagent <- c(indagent,1)
  indagent <- setdiff(indagent, grep("protocol",colnames(features)))
  
  treatind <- NULL
  controlind <- NULL
  indnon_drug <- NULL
  controlnames2 <- c("No drug","HFD","mDEN","SHAM","VEH","air","none treated","ctrl","normal","nodrug","unexposed","control"," none","untreated","mock", "uninduced","DMSO","saline","unexposed","vehicle","no treatment","contol","PBS","baseline","DEN","oil","placebo","no induction","EtOH","Control","Mannitol","Vehicle","Ctr","water","non.{0,1}treated","high fat", paste0("No.",input), "chow","before first","un-treated","naive",paste0("non.",input),"pre-treatment", paste0(input,".resistant"),paste0(input," sensitive"),paste0("without.",input),paste0("control.*",input),paste0("pre.",input),paste0(input,".*: N"),paste0(input," {0,1}\\-"),paste0("no.",input), paste0("before .*",input),paste0("prior .*",input),paste0(input,".* baseline"),"spontaneous differentiation")
  controlnames <- paste(controlnames2,collapse = "|")
  
  #paste0("control.*",input), paste0(input,"Veh")
  
  for(j in 1:length(indagent)){
    currentcolumn <- as.character(features[,indagent[j]])
    currenttreat <- grep(paste(as.vector(t(outer(paste0(c("^"," ","[[:punct:]]"),input), c("[[:punct:]]"," ","$"),paste, sep=""))),collapse = "|"),currentcolumn,ignore.case = T)
    currentcontrol <- grep(controlnames,currentcolumn,ignore.case = T)
    indnull <- which(nchar(currentcolumn) == 0)
    indnon_drug <- union(indnon_drug,grep(paste(c("No drug",paste0("no.",input),paste0("non.",input),paste0("pre.",input),"pre-treatment",paste0("without.",input),paste0("pre.",input),paste0(input,".*: N"),paste0("before .*",input),paste0(input,".deprivation"),paste0(input,".deficient"),paste0("prior .*",input),paste0(input,".* baseline")),collapse = "|"),currentcolumn,ignore.case = T))
    if(length(currenttreat) > 0 && length(currentcontrol) > 0){
      treatind <- union(treatind,currenttreat)
      controlind <- union(controlind,currentcontrol)
    }else if(length(currenttreat) > 0 && length(indnull) > 0){
      treatind <- union(treatind,currenttreat)
      controlind <- union(controlind,indnull)
    }
    if(length(controlind) != 0 && length(treatind) != 0 && !Agilent_label){
      break
    }#To stop if data is found in title, but special for agilent, check for other channels
  }
  
  if(length(controlind) == 0){
    controlnames2 <- c("CON","Veh","NT","Ctl")
    controlnames2 <- paste(controlnames2,collapse = "|")
    currentcolumn <- as.character(features[,indagent[1]])
    controlind <- grep(controlnames2,currentcolumn)
  }
  
  
  
  if(length(treatind) == 0){
    #for potential abbreviations based on first three words
    input2 <- substr(input,1,3)
    #input2 <- paste0(substr(input,1,3),"|",substr(input,nchar(input)-5,nchar(input)))
    #input2 <- toupper(input2)
    #input2 <- paste0(input2,"|",paste0(substr(input2,1,1),tolower(substr(input2,2,3))))
    #only run abbreviation in first column
    currentcolumn <- as.character(features[,indagent[1]])
    treatind <- grep(paste(as.vector(t(outer(paste0(c("^"," ","[[:punct:]]"),input2,"[[:alpha:]]{0,2}"), c("[[:punct:]]"," ","$"),paste, sep=""))),collapse = "|"),currentcolumn,ignore.case = T)
    if(length(controlind) == 0){
      controlind <- grep(controlnames,currentcolumn,ignore.case = T)
    }
    indnon_drug <- union(indnon_drug,grep(paste(c("No drug",paste0("no.",input2),paste0("non.",input2),paste0("pre.",input2),"pre-treatment",paste0("without.",input2),paste0("control.*",input2), paste0(input2,".* Veh"),paste0("before .*",input2),paste0(input2,".deprivation"),paste0(input2,".deficient"),paste0("prior .*",input2),paste0(input2,".* baseline")),collapse = "|"),currentcolumn,ignore.case = T))
  }
  
  controlind <- union(controlind,indnon_drug)
  if(length(indnon_drug) == nrow(features)){
    cat("all non-drug treatment in", gse_number)
    return(features)
  }
  
  #third phase acronym matching with extracting acronym from drug (acronym) style
  #if(length(treatind) == 0){
  #  currenttext <- paste0(pmid$title[which(pmid$gse %in% gse_number)]," ",pmid$summary[which(pmid$gse %in% gse_number)])
  #  acronym <- regmatches(currenttext, gregexpr(paste(c(tolower(paste(paste0(c(" ","^"),input," [[:alpha:]]* ?","\\(.*?\\)"),collapse = "|")),paste(paste0(c(" ","^"),input," [[:alpha:]]* ?","\\(.*?\\)"),collapse = "|")),collapse = "|"), currenttext))[[1]]
  #  if(length(acronym) != 0){
  #    acronym <- gsub("[\\(\\)]", "", regmatches(acronym , gregexpr("\\(.*?\\)", acronym ))[[1]])
  #    for(j in 1:length(indagent)){
  #      currentcolumn <- as.character(features[,indagent[j]])
  #      currenttreat <- grep(paste(as.vector(t(outer(paste0(c("^"," ","[[:punct:]]"),acronym), c("[[:punct:]]"," ","$"),paste, sep=""))),collapse = "|"),currentcolumn,ignore.case = T)
  #      if(length(controlind) > 0){
  #        currentcontrol <- controlind
  #      }else{
  #        currentcontrol <- grep(controlnames,currentcolumn,ignore.case = T)
  #      }
  #      if(length(currenttreat) > 0 && length(currentcontrol) > 0){
  #        treatind <- union(treatind,currenttreat)
  #        controlind <- union(controlind,currentcontrol)
  #        indnon_drug <- grep(paste(c("No drug",paste0("no.",acronym),paste0("non.",acronym),paste0("pre.",acronym),"pre-treatment",paste0("without.",acronym),paste0("before .*",acronym),paste0(acronym,".deprivation"),paste0(acronym,".deficient"),paste0("prior .*",acronym),paste0(acronym,".* baseline")),collapse = "|"),currentcolumn,ignore.case = T)
  #      }
  #      if(length(treatind) != 0){
  #        break
  #      }
  #    }
  #  }
  #}
  
  #this section was designed to detect non-drug, to exclude it, not only as control
  if(length(treatind) > 0){
    treatind <- setdiff(treatind,indnon_drug)
  }
  
  if(length(treatind) == nrow(features)){
    cat("all treated samples in", gse_number)
    return(features)
  }
  
  if(length(intersect(controlind,treatind)) > 0){
    if(Agilent_label){
      cat("duplicated label in Agilent array, it could be complex design ",gse_number)
      return(features)
    }
    cat("duplicated treatment control label", gse_number)
    if(length(controlind) > length(treatind)){
      controlind <- setdiff(controlind,treatind)
    }else if(length(controlind) < length(treatind)){
      treatind <- setdiff(treatind,controlind)
    }else if(mean(controlind %in% treatind) == 1){
      
      cat("same labels across the two sets")
      return(features)
    }else{
      cat("same labels across the two sets")
      Samples <- as.character(features$geo_accession)[treatind]
      Controls <- as.character(features$geo_accession)[controlind]
      return(features)
      
    }
  }
  
  if(length(controlind) == 0){
    cat("No control found", gse_number)
    Samples <- as.character(features$geo_accession)[treatind]
    return(list(Samples = Samples, Controls = NULL, features = features[c(Samples),], original = features))
    
  }
  if(length(treatind) == 0){
    cat("No treat found", gse_number)
    Controls <- as.character(features$geo_accession)[controlind]
    return(list(Samples = NULL, Controls = Controls, features = features[c(Controls),], original = features))
    
  }
  Samples <- as.character(features$geo_accession)[treatind]
  Controls <- as.character(features$geo_accession)[controlind]
  
  return(list(Samples = Samples, Controls = Controls, features = features[c(Samples,Controls),], original = features))
}

#testing codes, change input and gse number based on the frame

###conditions to look into
#if function returns 0 -> indicates poor quality data, you can just skip this
#if function returns only the table without sample/control, it indicates some complex design and cannot be identified easily, you can still look into this or skip it
#if function returns NULL in either sample or control, you should look into "original dataframe" to see if you can identify lost elements, if not, you should skip this
#Even if function works normally, you should check if Control and Sample are making sense (under the metadata dataframe)

############################################################################################################
#filtration code, please modify keywords for control or sample if needed
#####.   DONT run if you think it's looking fine!!!!  ########
features <- result[[4]]
#control filtration
new_controlword <- "0hr"
index <- grep(new_controlword,test$source_name_ch1,ignore.case = T)
control2 <- rownames(test[index,])
result[[2]] <- control2

#sample filtration
new_sampleword <- "BPA"
index <- grep(new_sampleword,test$source_name_ch1,ignore.case = T)
sample2 <- rownames(test[index,])
result[[1]] <- sample2
result[[3]] <- features[c(sample2,control2),]


control2 <- features$geo_accession[1:6]
result[[2]] <- control2
sample2 <- features$geo_accession[7:12]
result[[1]] <- sample2
result[[3]] <- features[c(sample2,control2),]

########################################################################################################################################

#Result summary:
#If it gave 0 -> go check next one, ignore this
#Full -> go 1 -> 2 -> 3
#If it gave non empty sample and control, check if both names are matching the meaning of filteredtable (3rd one)
#not full -> go 1 -> 2 -> 4
#If any of sample/control is empty, go check 4th element to see if you can find missing part.
###If you can find missing part, manually add stuffs to substitute the list
###Otherwiese make a note on a textfile, I will check it later

########################################################################################################################################
#######    EXAMPLE  formal running through all data  EXAMPLE  ########


dir.create("./Completed_datasets/Pentanedione/")
setwd("./Completed_datasets/Pentanedione/")

#input <- "Pentanedione"
#preview <- alldatasets[[currentdataset]]
#gse_number <- preview[1]
#result <- sampleidentification(input, gse_number,all_gses)

input <- "Pentanedione"
gse_number <- "GSE155845"
result <- sampleidentification(input, gse_number)

result[1]
result[2]
features_selected <- result[[3]]
features <- result[[4]]


#1, 2 (both contained things) -> 3
#1, 2 (empty) -> 4, note the datasets

### Save ### 
save(result, file = paste0(input,"_",gse_number,".rda"))


# DONE #
############################################################################################################
### Another way to grab treated and controls is to use the following code... this is not necessary## 

result$Samples <- rownames(features)[grep("Exposed", features$title, ignore.case = T)]
result$Controls <- rownames(features)[grep("Control|Vehicle|DMSO|Corn Oil|Naive",features$title, ignore.case = T)]
result$features  <- features[c(result$Samples,result$Controls),]

save(result, file = paste0(input,"_",gse_number,".rda"))