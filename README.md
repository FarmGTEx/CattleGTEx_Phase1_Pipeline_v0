# Analysis pipelines for the Phase1 of CattleGTEx
## 1. Introduction

Cattle are globally distributed across diverse environments and are fundamental to human societal development and global food security. Despite a pilot cattle multi-tissue regulatory atlas has been established, the genetic and molecular architecture underlying ecologically and economically im-portant traits in cattle remain largely unexplained, primarily due to limited sample size and tissue coverage, as well as a primary focus of gene expression alone. Here, we present the Cattle Geno-type–Tissue Expression (**CattleGTEx**) Phase 1 resource (https://cattlegtex.farmgtex.org/) as an essential component of the ongoing Farm animal Genotype-Tissue Expression (**FarmGTEx**, https://www.farmgtex.org/), which substantially expands the scale of the pilot study, increasing RNA-seq sample size (from 4,889 to 12,422), tissue coverage (from 23 to 43), and molecular phenotypes (from two to seven). Leveraging this expanded dataset, we identified 433,972 primary and 161,428 non-primary regulatory effects across seven molecular phenotypes and systematically characterized their sharing patterns across tis-sues and breeds. By integrating these regulatory data with genetic associations for 44 complex traits and selection signatures (Bos taurus vs. Bos indicus and dairy vs. beef breeds), we show that 76% of GWAS signals can be explained by this expanded resource, illustrating how CattleGTEx helps elucidate the regulatory mechanisms underlying adaptive evolution and selective breeding. Finally, by examining evolutionary constraints on these regulatory effects, we highlight the translational val-ue of this resource for understanding the shared genetic basis of complex traits in humans
![CattleGTEx-Phase1](https://github.com/FarmGTEx/CattleGTEx_phase1_Pipeline_v0/blob/main/CattleGTEx_Phase1.png)

## 2. Analysis pipelines

This repository contains codes used in data analysis of the SheepGTEx pilot phase, including:

- 01.RNA-seq: RNA-seq raw data processing, quantification and normalization of gene/exon/isoform/enhancer expression, splicing, RNA stability and 3'UTR APA, as well as ASE analysis.
- 02.metadata_prediction: breed prediction for RNA-seq samples.
- 03.gene_expression_analyses: gene expression analyses including sample clustering, pathway enrichment, tissue specificity and co-expression.
- 04.molQTL_mapping: molecular quantitative trait locus (molQTL) discovery, annotation and validation.
- 05.molQTL_colocalization: comparison between eQTL and other types of molQTL with the same gene.
- 06.molQTL_tissue_specificity: analyses of molQTL tissue sharing and effect size specificity. 
- 07.GWAS_and_molQTL: Enrichment analysis, Colocalzation, TWAS and SMR between molQTL and GWAS loci from 44 complex traits.
- 08.population_genetics: analyses of population structure (ADMIXTURE), signatures of selection (FST) and enrichment analysis between molQTL and selected signals.
- 09.Evolutionary_constraint_between_cattle_and_humans: alphagenome, LDSC, functional enrichment and colocalization based on human diseases and cattle eQTL.

# 3. Citation

#### **An enhanced multi-tissue atlas of regulatory effects in cattle**
Houcheng Li1†, Huicong Zhang1†, Di Zhu1†, Pengju Zhao2†, Zhenyu Wei3†, Jingsheng Lu4,5,6†, Mian Gong1,7, Qi Zhang1,8, Weijie Zheng1,8, Xinfeng Liu1,9, Dailu Guan10, Jinyan Teng11, Qin Lin11, Yongjie Tang8, Zhe Zhang11, Junting Du12, Chao Fang12, Bingxing An1,7, Bingjin Lin1,13, Min Tian1,14, Jingjing Tian1,7, Siqian Chen8, Wansheng Liu15, Yanan Wang4,5,6, Mingshan Wang4, Eveline M. Ibeagha-Awemu16, Richard P. M. A. Crooijmans17, Martijn F. L. Derks17, Marta Gòdia17, Ole Madsen17, Hubert Pausch18, Alexander S. Leonard18, Greger Larson19, Laurent Frantz20,21, David E. MacHugh22,23,24, John F. O’Grady18, Iuliana Ionita-Laza25, Xin Zhao26, Leluo Guan27, Huaijun Zhou10, Emilio Mármol-Sánchez28, Monique G. P. van der Wijst31,32, Xubin Lu33, Hui Jiang33, Zhangping Yang33, Qien Yang34, Qinyou Liu11, Chuang Xu35, Moli Li35, Yali Hou7, Yan Chen7, Ruidong Xiang36,37,38, Amanda J. Chamberlain36,39, Mathew Littlejohn40,41, Emily L. Clark42, Huiz-eng Sun43, Bo Han8, Yi Zhang8, Ying Yu8, Dongxiao Sun8, Yaokun Li44, Dewu Liu44, Goutam Sa-hana1, Zexi Cai1, Mogens Sandø Lund1, John B. Cole45,46,47, Li Ma48, Jicai Jiang49, Wenjian Li50, Yang Wu50, Qin Zhang51, Xiao Wang52*, Xuemei Lu4,5,6*, Yu Jiang14*, Yang Zhou53,54*, Yu Wang14*, Bingjie Li1,55*, Peter Sørensen1*, Lingzhao Fang1*
        
        
        
        
