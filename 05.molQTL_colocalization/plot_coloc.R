library(data.table)
library(qqman)
ARGS <- commandArgs(trailingOnly = TRUE)
QTL1 <- ARGS[1]
QTL2 <- ARGS[2]
QTL1 <- basename(QTL1)
QTL2 <- basename(QTL2)
#QTL1 <- "eQTL"
#QTL2 <- "sQTL"
path <- paste0("/faststorage/project/farmgtex/pipeline/coloc/coloc_",QTL1,"_",QTL2,"/")

tissues <- gsub("_coloc.csv","",list.files(path))
gtf <- fread("/faststorage/project/farmgtex/gtex/cattle_genome/Bos_taurus.ARS-UCD1.2.108.chr.gtf",data.table = F)
gtf <- gtf[gtf[,3] == "gene",]
for (i in tissues){
    if (!dir.exists(paste0(path,"/plot/",i,"/"))) {
      dir.create(paste0(path,"/plot/",i,"/"),recursive = TRUE)
    }
    coloc_sum <- fread(paste0(path,i,"_coloc.csv"),data.table  =F )
    result <- coloc_sum[grep("PP.H4.abf", coloc_sum[,1]), ]
    result <- result[order(result$res.summary,decreasing = T),]
    N <- 20
      for (n in 1:N){
      temp_result <- result[n,]
      chr <- gtf[grep(paste(temp_result$gene_id, collapse="|"), gtf[, 9]), 1]
      QTL1_pair_file <- list.files(path = paste0("/faststorage/project/farmgtex/QTL_result/", QTL1, "/", i, "/"),pattern = paste0("*.cis_qtl_pairs.", chr, ".txt.gz"),full.names = TRUE)
      QTL2_pair_file <- list.files(path = paste0("/faststorage/project/farmgtex/QTL_result/", QTL2, "/", i, "/"),pattern = paste0("*.cis_qtl_pairs.", chr, ".txt.gz"),full.names = TRUE)
      QTL1_pair <- fread(QTL1_pair_file,data.table = F)
      QTL2_pair <- fread(QTL2_pair_file,data.table = F)
      sub_QTL1 <- QTL1_pair[QTL1_pair$pheno_id == temp_result$QTL1,]
      sub_QTL1$CHR <- as.numeric(chr)
      sub_QTL2 <- QTL2_pair[QTL2_pair$pheno_id == temp_result$QTL2,]
      sub_QTL2$CHR <- as.numeric(chr)
      data1 <- sub_QTL1[,c("variant_id","CHR","start_distance","pval_g1")]
      colnames(data1) <- c("SNP", "CHR", "BP", "P")
      data2 <- sub_QTL2[,c("variant_id","CHR","start_distance","pval_g1")]
      colnames(data2) <- c("SNP", "CHR", "BP", "P")
      data1$BP <- 1:nrow(data1)
      data2$BP <- 1:nrow(data2)
      #data1$BP <- as.numeric(sub(".*_(.*)", "\\1", data1$SNP))
      #data2$BP <- as.numeric(sub(".*_(.*)", "\\1", data2$SNP))
      common_BP <- intersect(data1$BP, data2$BP)
      x_min <- min(common_BP, na.rm = TRUE)
      x_max <- max(common_BP, na.rm = TRUE)
      png(filename = paste0(path,"/plot/",i,"/manhattan_", QTL1, "_", QTL2, "_", i,"_",temp_result$gene_id,".png"), width = 3500, height = 3000, res = 300)
      par(mfrow = c(2, 1))

      manhattan(data1,
          col = c('#4682B4'),
          suggestiveline = -log10(1e-05),
          genomewideline = -log10(5e-08),
          annotatePval = 0.05,#标记p值小于0.05的点
          annotateTop = T, 
          #xlim = c(x_min, x_max),
          main = paste0(QTL1,"_",i,"_",temp_result$gene_id)#标题
          )
      manhattan(data2,
          col = c('#4682B4'),
          suggestiveline = -log10(1e-05),
          genomewideline = -log10(5e-08),
          annotatePval = 0.05,#标记p值小于0.05的点
          annotateTop = T,
         # xlim = c(x_min, x_max),
          main = paste0(QTL2,"_",i,"_",temp_result$gene_id)#标题 
          )  
          dev.off()            
      }
}