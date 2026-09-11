options(stringsAsFactors=FALSE)
.libPaths(c("C:/Users/ZCL/AppData/Local/R/win-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})
root <- Sys.getenv("SEPSIS_PROJECT_ROOT", unset = normalizePath("..", winslash = "/", mustWork = TRUE))
wdir <- file.path(root,"work","final_revision_20260908","wgcna_allgenes_rerun")
overlap <- readLines(file.path(wdir,"primary_sepsis_module_corrected_DEG_overlap.txt"),warn=FALSE)
assignments <- read.csv(file.path(wdir,"module_assignments_all_genes.csv"),check.names=FALSE)
mapped_overlap <- suppressMessages(bitr(overlap,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db))
mapped_universe <- suppressMessages(bitr(assignments$gene,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db))
ids <- unique(mapped_overlap$ENTREZID); universe <- unique(mapped_universe$ENTREZID)
go_obj <- enrichGO(gene=ids,universe=universe,OrgDb=org.Hs.eg.db,ont="ALL",pAdjustMethod="BH",
                   pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,readable=TRUE)
kegg_obj <- enrichKEGG(gene=ids,universe=universe,organism="hsa",keyType="kegg",pvalueCutoff=1,
                       pAdjustMethod="BH",qvalueCutoff=1,minGSSize=10)
if (nrow(as.data.frame(kegg_obj))) kegg_obj <- setReadable(kegg_obj,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
go <- as.data.frame(go_obj); kegg <- as.data.frame(kegg_obj)
write.csv(go,file.path(wdir,"overlap_GO_all_results.csv"),row.names=FALSE)
write.csv(kegg,file.path(wdir,"overlap_KEGG_all_results.csv"),row.names=FALSE)
write.csv(data.frame(input_overlap_n=length(overlap),mapped_overlap_n=length(ids),universe_n=length(assignments$gene),
                     mapped_universe_n=length(universe),GO_tested_n=nrow(go),GO_FDR_sig_n=sum(go$p.adjust<.05,na.rm=TRUE),
                     KEGG_tested_n=nrow(kegg),KEGG_FDR_sig_n=sum(kegg$p.adjust<.05,na.rm=TRUE)),
          file.path(wdir,"enrichment_summary.csv"),row.names=FALSE)
cat("GO",nrow(go),sum(go$p.adjust<.05,na.rm=TRUE),"KEGG",nrow(kegg),sum(kegg$p.adjust<.05,na.rm=TRUE),"\n")
