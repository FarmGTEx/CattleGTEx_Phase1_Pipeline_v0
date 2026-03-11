#!/bin/bash

cd /faststorage/project/cattle_gtexs/script/QTLenrich/QTLEnrich-master/src
main_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/state_ieQTL/result/"
for trait_dir in "$main_dir"/*/; do
  trait_name="$(basename "$trait_dir")"
  python3 QTLEnrichV2.py --gwas_file /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/GWAS/${trait_name}.txt \
                             --trait_name ${trait_name} \
                             --qtl_directory /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/tissue_shared/Group6/ \
                             --file_name .cis_qtl_pairs.significant.txt \
                             --qtl_type best_eqtl \
                             --confounders_table /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/tissue_shared/Group1_output_confounder.txt \
                             --null_table /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/tissue_shared/Group1_null_variants_table.txt \
                             --gencode_file /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/Bos_taurus.ARS-UCD1.2.gtf \
                             --gwas_p_value 0.05 \
                             --exp_label eQTL_${trait_name} \
                             --output_directory /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/tissue_shared/result/Group6/${trait_name}/ \
                             --GeneEnrich_input
done
