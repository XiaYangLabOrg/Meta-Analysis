# Meta-Analysis

The workflow for meta-analyis of bulk RNAseq data as established in manuscript "Systematic Transcriptome-wide Meta-Analysis across Endocrine Disrupting Chemicals Reveals Shared and Unique Liver Pathways, Gene Networks, and Disease Associations".


__1. Download of microarray data from GEO and sample finding__  

    2023_Microarray download-sample_finder code.R

__2. Normalization check and LIMMA__

    2023 MA and LIMMA Code.R

   _Files necessary for processing:_

           #1 old_HUGO_new_HUGO.rda

           #2 organ_search_codeZ.r

           #3 HUGO_symbols.txt

           #4 human_rat_hcop_fifteen_column.txt

           #5 human_mouse_hcop_fifteen_column.txt
           
__4. Download of RNAseq data from GEO (SRA) and pseudomapping using Salmon__

     2023 SRA download - Salmon.R

_Files necessary for processing:_

         #1 gse2srr.sh

__5. Importing quantification and DESeq2__

     2023 Tximport - DESeq2.R
     
_Depending on species one of the following files is necessary for processing:_

        #1 Common_Ratgenecode2gene.rda
        #2 Common_Mousegenecode2gene.rda
        #3 Common_Humangenecode2gene.rda
        
__6. Robust rank aggregation and enrichR (KEGG + DisGeNET)__

     2023_Meta-signature_enrichr_script.R

__7. wKDA analysis__

     2023_KDA.R
     
_Files necessary for processing:_

        #1 Mergeomics_Version_1.99.0.R 
   
