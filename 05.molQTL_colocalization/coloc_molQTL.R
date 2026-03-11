library(data.table)
library(coloc)
options(stringsAsFactors = FALSE)
ARGS <- commandArgs(trailingOnly = TRUE)
QTL1 <- ARGS[1]
QTL2 <- ARGS[2]
QTL1_file <- list.files(QTL1,full.name = FALSE)
QTL2_file <- list.files(QTL2,full.name = FALSE)
QTL1_tissue <- gsub("_coloc.bed","",QTL1_file)
QTL2_tissue <- gsub("_coloc.bed","",QTL2_file)
select_tissue <- intersect(QTL1_tissue,QTL2_tissue)
 if (!dir.exists(paste0("coloc_",basename(QTL1),"_",basename(QTL2)))) {
      dir.create(paste0("coloc_",basename(QTL1),"_",basename(QTL2)))
    }
#select_tissue <-"Gut"
if (!grepl("enhancer", QTL1) & !grepl("enhancer", QTL2)){    
 for (i in select_tissue){
  # for (i in select_tissue[(length(select_tissue)-1) :length(select_tissue)]){
    print(paste0("processing ",i," in ",basename(QTL1)," and ",basename(QTL2)))
    QTL1_result <- fread(paste0(QTL1,"/",i,"_coloc.bed"))
    QTL1_result$varbeta <- QTL1_result$beta_se^2
    QTL2_result <- fread(paste0(QTL2,"/",i,"_coloc.bed"))
    QTL2_result$varbeta <- QTL2_result$beta_se^2
    
    gene1 <- unique(QTL1_result$gene_id)
    gene2 <- unique(QTL2_result$gene_id) 
    #choose shared genes
    gene_list <- intersect(gene1,gene2)
    all_summary <- data.frame()
    all_result <- data.frame()
    if (length(gene_list) > 0){
        for (k in 1:length(gene_list)) {
          gene_qtl1 <- QTL1_result[QTL1_result$gene_id == gene_list[[k]], ]
          gene_qtl2 <- QTL2_result[QTL2_result$gene_id == gene_list[[k]], ]
            for( a in 1:length(unique(gene_qtl1$pheno_id))){
              for (b in 1:length(unique(gene_qtl2$pheno_id))){
          pheno_qtl1 <- gene_qtl1[gene_qtl1$pheno_id == unique(gene_qtl1$pheno_id)[a], ]
          pheno_qtl2 <- gene_qtl2[gene_qtl2$pheno_id == unique(gene_qtl2$pheno_id)[b], ]
          input <- merge(pheno_qtl1, pheno_qtl2, by="rs_id", all=FALSE, suffixes=c(paste0("_",basename(QTL1)),paste0("_",basename(QTL2)))) #change to your xQTL type
          input <- input[!duplicated(input$rs_id),]
          res <- coloc.abf(
          dataset1 = list(snp = input$rs_id,pvalues = input[[paste0("pval_nominal_", basename(QTL1))]], beta = input[[paste0("beta_", basename(QTL1))]], varbeta = input[[paste0("varbeta_", basename(QTL1))]],N = input[[paste0("N_sample_", basename(QTL1))]], type = "quant", MAF = input[[paste0("maf_", basename(QTL1))]]),dataset2 = list(snp = input$rs_id,pvalues = input[[paste0("pval_nominal_", basename(QTL2))]], beta = input[[paste0("beta_", basename(QTL2))]], varbeta = input[[paste0("varbeta_", basename(QTL2))]], N = input[[paste0("N_sample_", basename(QTL2))]], type = "quant", MAF = input[[paste0("maf_", basename(QTL2))]]))
          summary <- data.frame(res$summary)
          summary$gene_id <- gene_list[[k]]
          summary$QTL1 <- unique(gene_qtl1$pheno_id)[a]
          summary$QTL2 <- unique(gene_qtl2$pheno_id)[b]
          all_summary <- rbind(all_summary, summary)
          if (summary[6,1] >=0.5){
            result <- data.frame(res$result)
            result$gene_id <- gene_list[[k]]
            result$QTL1 <- unique(gene_qtl1$pheno_id)[a]
            result$QTL2 <- unique(gene_qtl2$pheno_id)[b]
            result <- result[result$SNP.PP.H4 >=0.5,]
            all_result <- rbind(all_result, result)
              }
            }
          }
        }
    }
          write.csv(all_summary,paste0("/faststorage/project/farmgtex/pipeline/coloc/coloc_",basename(QTL1),"_",basename(QTL2),"/",i,"_coloc2.csv"))
          write.csv(all_result,paste0("/faststorage/project/farmgtex/pipeline/coloc/coloc_",basename(QTL1),"_",basename(QTL2),"/",i,"_coloc_result2.csv"))
  }
} else {
    if (grepl("enhancer", QTL1)) {
    enhancer <- QTL1
    }else {
    enhancer <- QTL2 
    QTL2 <- QTL1 
    }
    for (i in select_tissue){
    #for (i in select_tissue[(length(select_tissue)-1) :length(select_tissue)]){
    print(paste0("processing ",i," in ",basename(QTL1)," and ",basename(QTL2)))
    enhancer_result <- fread(paste0(enhancer,"/",i,"_coloc.bed"))
    enhancer_result$varbeta <- enhancer_result$beta_se^2
    QTL2_result <- fread(paste0(QTL2,"/",i,"_coloc.bed"))
    QTL2_result$varbeta <- QTL2_result$beta_se^2
    gene_list <- unique(QTL2_result$gene_id) 
    all_summary <- data.frame()
    all_result <- data.frame()
    for (k in 1:length(gene_list)) {
          gene_qtl2 <- QTL2_result[QTL2_result$gene_id == gene_list[[k]], ]
          gene_enhancer <- enhancer_result[enhancer_result$rs_id %in% intersect(enhancer_result$rs_id,gene_qtl2$rs_id),]
          if (nrow(gene_enhancer) > 0){
            for (b in 1:length(unique(gene_qtl2$pheno_id))){
              for( a in 1:length(unique(gene_enhancer$pheno_id))){    
                pheno_qtl1 <- gene_enhancer[gene_enhancer$pheno_id == unique(gene_enhancer$pheno_id)[a], ]
                pheno_qtl2 <- gene_qtl2[gene_qtl2$pheno_id == unique(gene_qtl2$pheno_id)[b], ]
                input <- merge(pheno_qtl1, pheno_qtl2, by="rs_id", all=FALSE, suffixes=c(paste0("_",basename(enhancer)),paste0("_",basename(QTL2)))) #change to your xQTL type
                input <- input[!duplicated(input$rs_id),]
                res <- coloc.abf(
          dataset1 = list(snp = input$rs_id,pvalues = input[[paste0("pval_nominal_", basename(enhancer))]], beta = input[[paste0("beta_", basename(enhancer))]], varbeta = input[[paste0("varbeta_", basename(enhancer))]],N = input[[paste0("N_sample_", basename(enhancer))]], type = "quant", MAF = input[[paste0("maf_", basename(enhancer))]]),dataset2 = list(snp = input$rs_id,pvalues = input[[paste0("pval_nominal_", basename(QTL2))]], beta = input[[paste0("beta_", basename(QTL2))]], varbeta = input[[paste0("varbeta_", basename(QTL2))]], N = input[[paste0("N_sample_", basename(QTL2))]], type = "quant", MAF = input[[paste0("maf_", basename(QTL2))]]))
          summary <- data.frame(res$summary)
          summary$gene_id <- gene_list[[k]]
          summary$QTL1 <- unique(gene_enhancer$pheno_id)[a]
          summary$QTL2 <- unique(gene_qtl2$pheno_id)[b]
          all_summary <- rbind(all_summary, summary)
          if (summary[6,1] >=0.5){
          result <- data.frame(res$result)
          result$gene_id <- gene_list[[k]]
          result$QTL1 <- unique(gene_enhancer$pheno_id)[a]
          result$QTL2 <- unique(gene_qtl2$pheno_id)[b]
          result <- result[result$SNP.PP.H4 >=0.5,]
          all_result <- rbind(all_result, result)
              }
            }
          }
        }
    
    }
    write.csv(all_summary,paste0("/faststorage/project/farmgtex/pipeline/coloc/coloc_",basename(enhancer),"_",basename(QTL2),"/",i,"_coloc2.csv"))
    write.csv(all_result,paste0("/faststorage/project/farmgtex/pipeline/coloc/coloc_",basename(enhancer),"_",basename(QTL2),"/",i,"_coloc_result2.csv"))
  }
}
