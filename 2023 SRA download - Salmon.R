######### Download and processing of RNAseq Datasets is done through the Hoffman2 Cluster at UCLA ########

#### batch processing codes for downloading GEO datasets and searching for SRRs/SRXs ####
setwd("/Users/zac/Box Sync/Thesis/RNAseq_Salmon/")
input <- "GSE148809" 
species <- "mouse_index"  #Index here 
homefolder <- "/u/home/z/zaczamor/"
fastqfolder <- "/u/scratch/z/zaczamor/fastq/"
command1 = paste0("bash $HOME/gse2srr.sh ",input)
system(command1)
SRAfile <- readLines(paste0(homefolder,input,"_SRR.txt")) 

#make sure your gse2SRR.txt is in home directory
## Open the SRRtext if you want to modify specific SRRs/SRXs to be downloaded in cyberduck (not necessary) ##
#file.remove(paste0(homefolder,input,"_SRR.txt"))
allSRRs <- NULL
allSRXs <- NULL
for(SRA in SRAfile){
  if(startsWith(SRA,"SRR")){allSRRs <- c(allSRRs, SRA)}
  if(startsWith(SRA,"SRX")){allSRXs <- c(allSRXs, SRA)}
}
dir.create(paste0(fastqfolder,input))
length(allSRRs) == length(allSRXs)
## if length doesn't match make sure to check if dataset is RNASeq dataset/ or if there are errors

#To run on interactive mode qrsh -l h_rt=5:55:55,h_data=12G, this may need to be modified based on size of dataset 
#### Downloading, trimming, and salmon quasi-mapping ####
for(i in 1:length(allSRXs)){
  correspondingSRRs <- unlist(strsplit(allSRRs[i],"\t",fixed = T))
  if(length(correspondingSRRs) > 1){
    for(j in 1:length(correspondingSRRs)){
      command = paste0("parallel-fastq-dump -s ",correspondingSRRs[j]," --outdir $SCRATCH/fastq --split-files --threads 4 --gzip --tmpdir $SCRATCH/fastq")
      system(command)
    }
    allfiles <- list.files(path = fastqfolder,pattern = paste0(correspondingSRRs,collapse = "|"))
    ## Merge the reads from different runs
    if(length(grep("_2",allfiles)) > 0){
      system(paste0("cat ",paste0("$SCRATCH/fastq/",allfiles[grep("_1",allfiles)],collapse = " ",sep = "")," > ",paste0("$SCRATCH/fastq/",correspondingSRRs[1]),"final_1.fastq.gz"))
      system(paste0("cat ",paste0("$SCRATCH/fastq/",allfiles[grep("_2",allfiles)],collapse = " ",sep = "")," > ",paste0("$SCRATCH/fastq/",correspondingSRRs[1]),"final_2.fastq.gz"))
    }else{
      system(paste0("cat ",paste0("$SCRATCH/fastq/",allfiles,collapse = " ",sep = "")," > ",paste0("$SCRATCH/fastq/",correspondingSRRs[1]),"final.fastq.gz"))
    }
    system(paste0("rm ",paste0("$SCRATCH/fastq/",allfiles,collapse = " ",sep = "")))
  }else{
    command = paste0("parallel-fastq-dump -s ",correspondingSRRs," --outdir $SCRATCH/fastq --split-files --threads 4 --gzip --tmpdir $SCRATCH/fastq")
    system(command)
  }
  allfiles <- list.files(path = fastqfolder,pattern = correspondingSRRs[1])
  if(length(allfiles) == 0){
    system(command)
    allfiles <- list.files(path = fastqfolder,pattern = paste0(correspondingSRRs[1],"_[0-9].fastq.gz"))
    if(length(allfiles) == 0){stop("SRR is not matching with SRX")}
  }
  status <- ifelse(length(allfiles) == 2, "paired", "single")
  prefix <- gsub("_.*","",allfiles)
  dir.create(paste0(fastqfolder,input,"/",allSRXs[i]))
  if(status == "paired"){
    system(paste0("trim_galore --paired $SCRATCH/fastq/",allfiles[1]," $SCRATCH/fastq/",allfiles[2]," -o $SCRATCH/fastq -j 4"))
    allfiles2 <- paste0(prefix,c("_1_val_1.fq.gz","_2_val_2.fq.gz"))
    system(paste0("/u/home/z/zaczamor/salmon-1.7.0_linux_x86_64/bin/salmon quant -i /u/home/z/zaczamor/", species, " -l A -1 $SCRATCH/fastq/",allfiles2[1], " -2 $SCRATCH/fastq/",allfiles2[2]," -p 4 -o $SCRATCH/fastq/",input,"/",allSRXs[i]))
  }else{
    system(paste0("trim_galore $SCRATCH/fastq/",allfiles[1]," -o $SCRATCH/fastq -j 4"))
    if(length(list.files(path = "/u/scratch/z/zaczamor/fastq",pattern = "_1_trimmed.fq.gz")) > 0){
      allfiles2 <- paste0(prefix,"_1_trimmed.fq.gz")
    }else{
      allfiles2 <- list.files(path = fastqfolder,pattern = ".*trimmed.fq.gz")
    }
    system(paste0("/u/home/z/zaczamor/salmon-1.7.0_linux_x86_64/bin/salmon quant -i /u/home/z/zaczamor/", species, " -l A -r $SCRATCH/fastq/",allfiles2[1]," -p 4 -o $SCRATCH/fastq/",input,"/",allSRXs[i]))
 }
  file.remove(paste0("/u/scratch/z/zaczamor/fastq/",c(allfiles,allfiles2)))}

## End result should have "quant.sf" file in scratch under SRX directory ##