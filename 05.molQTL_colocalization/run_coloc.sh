realpath */ | grep "QTL"|grep -v "coloc_" | awk '{folders[NR]=$0} END {for (i=1; i<=NR; i++) for (j=i+1; j<=NR; j++) print folders[i], folders[j]}'|grep -E "eQTL|sQTL" |grep "eQTL"|grep -v "sex"|while read qtl1  qtl2
do
sbatch -p normal -J coloc -n 4 --mem 100g -t 24:00:00 --wrap="Rscript coloc_XQTL.R ${qtl1} ${qtl2}"
done

realpath */ | grep -v "coloc_" | grep -v "sex" | awk '{folders[NR]=$0} END {for (i=1; i<=NR; i++) for (j=i+1; j<=NR; j++) print folders[i], folders[j]}'|grep  "eQTL" | grep "enhancer" |while read qtl1  qtl2
do
sbatch -p normal -J coloc -n 4 --mem 100g -t 24:00:00 --wrap="Rscript coloc_XQTL.R ${qtl1} ${qtl2}"
done
 
realpath */ | grep -v "coloc_"  | grep -v "sex"| awk '{folders[NR]=$0} END {for (i=1; i<=NR; i++) for (j=i+1; j<=NR; j++) print folders[i], folders[j]}'|grep "enhancer" |grep "sQTL" |while read qtl1  qtl2
do
sbatch -p normal -J coloc -n 4 --mem 100g -t 24:00:00 --wrap="Rscript coloc_XQTL.R ${qtl1} ${qtl2}"
done