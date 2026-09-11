options(stringsAsFactors=FALSE)
old_file <- "C:/Users/ZCL/sepsis_workspace_link/WGCNA/wgcna_analysis_complete.RData"
out_dir <- "C:/Users/ZCL/sepsis_workspace_link/work/final_revision_20260908/old_wgcna_extract"
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
load(old_file)
corv <- as.numeric(as.matrix(module_group_cor))
pv <- as.numeric(as.matrix(module_group_p))
nms <- rownames(module_group_cor)
if (is.null(nms) || length(nms)!=length(corv)) nms <- rownames(as.matrix(module_group_cor))
nms <- sub("^ME","",nms)
sizes <- as.data.frame(table(module=as.character(moduleColors)),stringsAsFactors=FALSE)
names(sizes)[2] <- "module_size"
tab <- data.frame(module=nms,correlation_sepsis=corv,p_value=pv,FDR=p.adjust(pv,method="BH"),stringsAsFactors=FALSE)
tab <- merge(tab,sizes,by="module",all.x=TRUE,sort=FALSE)
tab <- tab[order(-abs(tab$correlation_sepsis)),]
write.csv(tab,file.path(out_dir,"old_module_trait_all_results.csv"),row.names=FALSE)
write.csv(data.frame(parameter=c("sample_n","sepsis_n","normal_n","gene_n","soft_threshold_power","networkType","TOMType","minModuleSize","deepSplit","mergeCutHeight","correlation","multiple_testing","random_seed","sample_exclusion"),
 value=c(ncol(expr_data),sum(group_info$GroupNum==1),sum(group_info$GroupNum==0),nrow(expr_data),softPower,"signed","signed",minModuleSize,2,MEDissThres,"Pearson","BH reconstructed across all module tests","not set",paste(outliers,collapse=";"))),file.path(out_dir,"old_parameters.csv"),row.names=FALSE)
print(tab)
