
GSE=$1
echo "Finding all SRX associated with ${GSE}..."

mapfile -t SRX_ARRAY < <(esearch -db gds -query "${GSE}[ACCN] AND GSM[ETYP]" |\
efetch -format docsum | xtract -pattern ExtRelation -element TargetObject)

echo "Finding all SRR associated with ${GSE}..."

rm -f ${GSE}_SRR.txt

for i in "${SRX_ARRAY[@]}"
do
   echo "$i" >> ${GSE}_SRR.txt 
   esearch -db sra -query $i | efetch -format docsum | \
   xtract -pattern DocumentSummary -element Run@acc >> ${GSE}_SRR.txt
done

