#!/bin/bash

cd /faststorage/project/cattle_gtexs/CattleGTEx/Ancient_DNA
vcftools --gzvcf ancient_cattle_20.vcf.gz --weir-fst-pop modern_euro.txt --weir-fst-pop ancient_euro.txt --out euro_fst.txt --fst-window-size 30000 --fst-window-step 10000
