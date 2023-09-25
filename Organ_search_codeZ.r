#tissue annotation function
load("/Users/zac/Box Sync/Thesis/BTO_tissue_converter_new_1202.rda")
load("/Users/zac/Box Sync/Thesis/BTO_include_synonym_pruned_1202.rda")
load("/Users/zac/Box Sync/Thesis/all_gse.rda")
processed_BrendaOnto2$X2 <- gsub("[(]","",processed_BrendaOnto2$X2)
processed_BrendaOnto2$X2 <- gsub("[)]","",processed_BrendaOnto2$X2)
processed_BrendaOnto3 <- processed_BrendaOnto2
indcell <- union(grep("^[A-Z].+ cell$",processed_BrendaOnto3$X2),grep("^[0-9].+ cell$",processed_BrendaOnto3$X2))
processed_BrendaOnto3$X2[indcell] <- gsub(" cell$","",processed_BrendaOnto3$X2[indcell])

indadd <- grep("^[0-9]+$",processed_BrendaOnto3$X2)
processed_BrendaOnto3$X2[indadd] <- paste0(processed_BrendaOnto3$X2[indadd]," cell")


#Brenda 3 was to substitute cell into non, because in annotation usually "cell" is omitted

finalitems <- c("kidney", "liver" ,"pancreas", "breast", "ovary", "adipose tissue" ,
                "cardiovascular system","nervous system","respiratory system","urogenital system",
                "immune system","hematopoietic system","skeletal system",
                "integument","connective tissue","muscular system","gland",
                "soft body part","viscus",
                "primary cell","pacemaker cell","adult stem cell","embryonic structure","tissues")


Brenda_onto_annot <- function(example,processed_BrendaOnto2,processed_BrendaOnto3,finalitems,minus_exclusion = TRUE){
  
  if(minus_exclusion){example <- gsub("_|:|-|_|,|\\."," ",example)
  }else{example <- gsub("_|:|_|,|\\."," ",example)
  }
  
  allind <- 1:nrow(processed_BrendaOnto2)
  max <- 25
  d1 <- split(allind, ceiling(allind/max))
  
  
  final <- NULL
  matchedword <- NULL
  finalind <- NULL
    ##chunksize grep and identification
    for(k in 1:length(d1)){
      
      
      testing_ind <- grep(
        paste(as.vector(t(outer(
          as.vector(t(outer(c("^"," "), processed_BrendaOnto3$X2[d1[[k]]], paste, sep=""))) 
          , c("s{0,1} ","s{0,1}$"), paste, sep=""))) ,
              collapse = "|"), example,ignore.case = T)
      
      if(length(testing_ind) == 0){
        next
      }else{
        finalind <- c(finalind,d1[[k]])
      }
    }
  
  #matching of mapped words
    if(length(finalind) > 0){
      for(i in finalind){
        currentword <- processed_BrendaOnto3[i,]
        ind <- NULL
        if(i %in% indcell){
        try(
          {ind <- grep(paste(as.vector(t(outer(paste0(c("^"," "),currentword$X2), c("s{0,1} ","s{0,1}$"),paste, sep=""))),collapse = "|"), example)
          })
          }else{
            try(
              {ind <- grep(paste(as.vector(t(outer(paste0(c("^"," "),currentword$X2), c("s{0,1} ","s{0,1}$"),paste, sep=""))),collapse = "|"), example,ignore.case = T)
              })
            }
        
        if(length(ind) > 0){
          final <- c(final,as.character(BTO_tissue_converter[currentword$X1,]))
          matchedword <- c(matchedword,currentword$X2)
        }
      }
    }
    
    #if not found in first path, continue searching title and summary
    if(length(final) == 0){
      text_ind <- TRUE
      GSE <- example[2]
      if(is.na(GSE)){
        return(NA)
      }
      GSE <- paste(all_gses[which(all_gses$gse %in% GSE),c(2,8)],collapse = " ")
      GSE <- gsub("_|:|-|_|,|\\."," ",GSE)
      #chunksize identification first, with more stringent version 2 Bernda
      for(k in 1:length(d1)){
        testing_ind <- grep(paste(paste0(" ",processed_BrendaOnto2$X2[d1[[k]]],"s{0,1} "),collapse = "|"), GSE, ignore.case = T)
        if(length(testing_ind) == 0){
          next
        }else{
          finalind <- c(finalind,d1[[k]])
        }
      }
      
      for(i in finalind){
        currentword <- processed_BrendaOnto2[i,]
        ind <- NULL
        try(
          {#allow some modifications here
            ind <- grep(paste0(" ",currentword$X2,"s{0,1} "), GSE, ignore.case = T)
          }
        )
        if(length(ind) > 0){
          final <- c(final,as.character(BTO_tissue_converter[currentword$X1,]))
          matchedword <- c(matchedword,currentword$X2)
        }
      }
    }
    
    if(length(final) == 0){
      final <- NA
    }else{
      final <- unique(final)
      final <- final[!is.na(final)]
      if(length(final) > 1){
        #make sorting to have single representation
        final <- factor(final,levels = finalitems)
        final <- sort(final)
        final <- as.character(final)[1]
      }
      if(final %in% "tissues"){final <- "unclassified tissues"}
    

    return(final)
    }
}

