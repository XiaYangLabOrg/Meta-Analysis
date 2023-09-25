# Meta-Analysis

The workflow as established in manuscript "Systematic Transcriptome-wide Meta-Analysis across Endocrine Disrupting Chemicals Reveals Shared and Unique Liver Pathways, Gene Networks, and Disease Associations".


1. Download of microarray data from GEO and sample finding  

    2023_Microarray download-sample_finder code.R

2. Normalization check and LIMMA

    2023 MA and LIMMA Code.R

   Files necessary for processing:

           #1 old_HUGO_new_HUGO.rda

           #2 organ_search_codeZ.r

           #3 HUGO_symbols.txt

           #4 human_rat_hcop_fifteen_column.txt

           #5 human_mouse_hcop_fifteen_column.txt
           
4. Download of RNAseq data from GEO (SRA) and pseudomapping using Salmon

     2023 SRA download - Salmon.R

5. Importing quantification and DESeq2

     2023 Tximport - DESeq2.R

6. Robust rank aggregation and enrichR (KEGG + DisGeNET)

     2023_Meta-signature_enrichr_script.R

7. wKDA analysis

     2023_KDA.R
   
