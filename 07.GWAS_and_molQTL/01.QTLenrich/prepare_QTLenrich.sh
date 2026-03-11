#!/bin/bash

cd /faststorage/project/cattle_gtexs/script/QTLenrich/

python3 extract_unique_variants.py --directory /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2/ --file_name .cis_qtl_pairs.significant.txt --output_file /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_significant_variants.txt

python3 create_null_variants.py --all_variants_file /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/all_variant.txt --significant_variants_file /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_significant_variants.txt --output_file /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_null_variants.txt

python3 Generate_Null_Table.py --QTL_Directory /faststorage/project/farmgtex/Figure6/QTLenrich/eQTL/all_pairs/ --File_Extension .eqtl_allpairs.txt --Null_Variants_List /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_null_variants.txt --Output_File /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_null_variants_table.txt.gz

python3 Generate_Confounders_Table.py --QTL_Directory /faststorage/project/farmgtex/Figure6/QTLenrich/eQTL/all_pairs/ --File_Extension .eqtl_allpairs.txt --Variants_List /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/all_variant.txt --LD_Proxy /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/number_ld_per_variant.txt --GENCODE_File /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/Bos_taurus.ARS-UCD1.2.gtf --Output_File /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_output_confounder.txt.gz
