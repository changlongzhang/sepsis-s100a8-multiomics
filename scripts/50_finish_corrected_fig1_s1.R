options(stringsAsFactors = FALSE)
.libPaths(c("C:/Users/ZCL/AppData/Local/R/win-library/4.5", .libPaths()))
suppressPackageStartupMessages({
  library(limma); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggrepel); library(clusterProfiler); library(org.Hs.eg.db)
  library(glmnet); library(randomForest); library(mboost); library(xgboost)
  library(patchwork); library(svglite); library(ragg); library(grid)
})

root <- "C:/Users/ZCL/sepsis_workspace_link"
panel_dir <- file.path(root, "review_revision", "10_figure_sources", "Fig1_S1_corrected_panels")
formal_dir <- "C:/Users/ZCL/sepsis_final_figures"
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

font <- "Arial"; ink <- "#20252A"; mid <- "#7A858D"; light <- "#D9DEE2"; pale <- "#F2F4F5"
navy <- "#315B78"; teal <- "#5B958C"; coral <- "#D56A54"
red <- "#D66A5E"; blue <- "#4C78A8"; green <- "#59A14F"

theme_pub <- function(base_size=8) theme_classic(base_size=base_size, base_family=font) + theme(
  text=element_text(colour=ink), axis.text=element_text(size=8, colour=ink),
  axis.title=element_text(size=8.5, colour=ink), axis.line=element_line(linewidth=.35, colour=ink),
  axis.ticks=element_line(linewidth=.3, colour=ink), legend.title=element_text(size=8),
  legend.text=element_text(size=8), plot.tag=element_text(size=12, face="bold", family=font),
  plot.tag.position=c(0,1), plot.margin=ggplot2::margin(4,5,4,5)
)
save_panel <- function(p, stem, w, h) {
  ggsave(file.path(panel_dir,paste0(stem,".pdf")),p,width=w,height=h,device=cairo_pdf,family=font,bg="white")
  ggsave(file.path(panel_dir,paste0(stem,".svg")),p,width=w,height=h,device=svglite,bg="white")
}

cat("STAGE\tLOAD_CORRECTED_RESULTS\n")
load(file.path(panel_dir,"Fig1_corrected_analysis.RData"))
expr_linear <- read.csv(file.path(root,"sepsis-s100a8-multiomics","data","02_DEG","Input","GSE65682_gene.csv"),row.names=1,check.names=FALSE)
expr <- log2(as.matrix(expr_linear))
meta <- read.csv(file.path(root,"sepsis-s100a8-multiomics","data","02_DEG","Input","GSE65682_Groups.csv"),check.names=FALSE)
colnames(meta)[1:2] <- c("sample_id","group")
group <- factor(meta$group[match(colnames(expr),meta$sample_id)],levels=c("Normal","Disease"))
black <- unique(trimws(readLines(file.path(root,"sepsis-s100a8-multiomics","data","03_WGCNA_analysis","output","DIFF_MODULE_black_GENES.txt"))))
black <- black[nzchar(black)]
stopifnot(length(overlap)==200,"S100A8"%in%overlap)

cat("STAGE\tREBUILD_KEGG\n")
mapped <- suppressMessages(bitr(overlap,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db))
kegg_result <- enrichKEGG(gene=unique(mapped$ENTREZID),organism="hsa",pvalueCutoff=1,pAdjustMethod="BH",minGSSize=10)
kegg_result <- setReadable(kegg_result,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
kegg_all <- as.data.frame(kegg_result)
write.csv(kegg_all,file.path(panel_dir,"Fig1F_KEGG_all_results.csv"),row.names=FALSE)
kegg_top <- kegg_all %>% arrange(pvalue) %>% slice_head(n=5) %>% mutate(
  Category=case_when(
    grepl("Complement|Fc epsilon|TNF|IL-17|immune|phagocytosis",Description,ignore.case=TRUE)~"Immune",
    grepl("metabolism|biosynthesis|steroidogenesis",Description,ignore.case=TRUE)~"Metabolism",
    TRUE~"Related"),
  Description_plot=str_wrap(Description,width=22)) %>% arrange(Category,Count) %>%
  mutate(Description_plot=factor(Description_plot,levels=unique(Description_plot)))
pF1 <- ggplot(kegg_top,aes(Count,Description_plot,fill=Category))+geom_col(width=.62)+
  geom_text(aes(label=Count),hjust=-.25,size=8/ggplot2::.pt,family=font)+
  scale_fill_manual(values=c(Immune=red,Metabolism=green,Related=blue),drop=FALSE)+
  scale_x_continuous(expand=expansion(mult=c(0,.2)))+
  labs(title="KEGG Enrichment Analysis",x="Gene count",y=NULL,fill=NULL,caption="No pathway passed FDR < 0.05")+theme_pub()+
  theme(axis.text.y=element_text(size=8),axis.ticks.y=element_blank(),legend.position="none",
        panel.grid.major.x=element_line(linewidth=.25,colour="#E7EBEE"),
        plot.title=element_text(size=8.5,face="bold",hjust=.5,margin=ggplot2::margin(b=5)),
        plot.caption=element_text(size=8,hjust=0,colour=mid),plot.margin=ggplot2::margin(4,4,3,4))
save_panel(pF1,"Fig1F_KEGG_corrected",2.72,2.42)
cat("KEGG_ROWS\t",nrow(kegg_all),"\n",sep="")
cat("KEGG_FDR_SIG_N\t",sum(kegg_all$p.adjust<.05,na.rm=TRUE),"\n",sep="")

go_top <- go_all %>% group_by(ONTOLOGY) %>% arrange(pvalue,.by_group=TRUE) %>% slice_head(n=5) %>% ungroup() %>%
  mutate(ONTOLOGY=factor(ONTOLOGY,levels=c("MF","CC","BP")),Description_plot=str_wrap(Description,width=45)) %>%
  arrange(ONTOLOGY,Count) %>% mutate(Description_plot=factor(Description_plot,levels=unique(Description_plot)))
pE1 <- ggplot(go_top,aes(Count,Description_plot,fill=ONTOLOGY))+geom_col(width=.68)+
  geom_text(aes(label=Count),hjust=-.25,size=8/ggplot2::.pt,family=font)+
  scale_fill_manual(values=c(BP=red,CC=blue,MF=green))+
  scale_x_continuous(expand=expansion(mult=c(0,.17)))+
  labs(title="GO Enrichment Analysis",x="Gene count",y=NULL,fill=NULL)+theme_pub()+
  theme(axis.text.y=element_text(size=8,lineheight=.75),axis.ticks.y=element_blank(),
        panel.grid.major.x=element_line(linewidth=.25,colour="#E7EBEE"),
        legend.position=c(.99,.02),legend.justification=c(1,0),legend.background=element_blank(),
        legend.key.height=unit(2.6,"mm"),legend.spacing.y=unit(0,"mm"),
        plot.title=element_text(size=8.5,face="bold",hjust=.5,margin=ggplot2::margin(b=5)),
        plot.margin=ggplot2::margin(4,4,4,4))
save_panel(pE1,"Fig1E_GO_corrected",4.17,2.42)

cat("STAGE\tFIT_S1_FINAL_MODELS\n")
genes <- intersect(overlap,rownames(expr))
x <- as.data.frame(t(expr[genes,,drop=FALSE]),check.names=FALSE)
x <- x[,apply(x,2,sd,na.rm=TRUE)>0,drop=FALSE]
for(j in seq_len(ncol(x))) if(anyNA(x[[j]])) x[[j]][is.na(x[[j]])] <- median(x[[j]],na.rm=TRUE)
xs <- scale(x); y <- as.numeric(group=="Disease")

set.seed(123); lf <- glmnet(xs,y,family="binomial",alpha=1,maxit=10000)
set.seed(123); lcv <- cv.glmnet(xs,y,family="binomial",alpha=1,nfolds=10,maxit=10000)
lc <- as.matrix(coef(lf,s=lcv$lambda.min))
lasso_imp <- data.frame(Gene=rownames(lc),Value=abs(lc[,1])) %>% filter(Gene!="(Intercept)",Value>0) %>% arrange(desc(Value))

set.seed(123); rf <- randomForest(x=xs,y=group,ntree=500,mtry=min(140,ncol(xs)),importance=TRUE)
ri <- randomForest::importance(rf)
metric <- if("MeanDecreaseGini"%in%colnames(ri)) "MeanDecreaseGini" else colnames(ri)[1]
rf_imp <- data.frame(Gene=rownames(ri),Value=ri[,metric]) %>% arrange(desc(Value))

bd <- data.frame(Group=group,as.data.frame(xs,check.names=FALSE),check.names=FALSE)
set.seed(123); bf <- glmboost(Group~.,data=bd,family=Binomial(),control=boost_control(mstop=250,nu=.1))
bc <- coef(bf,off2int=TRUE)
boost_imp <- data.frame(Gene=names(bc),Value=abs(as.numeric(bc))) %>% filter(Gene!="(Intercept)",Value>0) %>% arrange(desc(Value))

dtrain <- xgb.DMatrix(as.matrix(xs),label=y)
set.seed(123); xf <- xgb.train(data=dtrain,nrounds=50,verbose=0,params=list(
  objective="binary:logistic",eval_metric="auc",max_depth=1,eta=.4,gamma=0,
  colsample_bytree=.6,min_child_weight=1,subsample=.5,nthread=6))
xi <- xgb.importance(feature_names=colnames(xs),model=xf)
xgb_imp <- data.frame(Gene=xi$Feature,Value=xi$Gain) %>% arrange(desc(Value))

contrib <- predict(xf,dtrain,predcontrib=TRUE)
feature_cols <- setdiff(colnames(contrib),"BIAS")
shap_all <- as.data.frame(contrib[,feature_cols,drop=FALSE],check.names=FALSE)
baseline <- mean(contrib[,"BIAS"])
top6 <- head(names(sort(colMeans(abs(as.matrix(shap_all))),decreasing=TRUE)),6)
shap_top <- shap_all[,top6,drop=FALSE]
x_top <- as.data.frame(xs[,top6,drop=FALSE],check.names=FALSE)

rank_df <- bind_rows(
  head(lasso_imp,5)%>%mutate(Model="LASSO"), head(boost_imp,5)%>%mutate(Model="glmBoost"),
  head(xgb_imp,5)%>%mutate(Model="XGBoost"), head(rf_imp,5)%>%mutate(Model="Ranger")) %>% dplyr::select(Gene,Model,Value)
write.csv(rank_df,file.path(panel_dir,"S1_model_rankings_corrected.csv"),row.names=FALSE)
write.csv(data.frame(Sample=rownames(shap_top),shap_top),file.path(panel_dir,"S1_shap_values_corrected.csv"),row.names=FALSE)
write.csv(data.frame(Sample=rownames(x_top),x_top),file.path(panel_dir,"S1_feature_values_corrected.csv"),row.names=FALSE)
write.csv(data.frame(Baseline_log_odds=baseline),file.path(panel_dir,"S1_shap_baseline_corrected.csv"),row.names=FALSE)
write.csv(data.frame(lambda=lcv$lambda,ll=log(lcv$lambda),cvm=lcv$cvm,cvsd=lcv$cvsd,cvup=lcv$cvup,cvlo=lcv$cvlo,nzero=lcv$nzero),
          file.path(panel_dir,"S1_lasso_cv_corrected.csv"),row.names=FALSE)
path_df <- as.data.frame(as.matrix(lf$beta)) %>% tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene,names_to="Step",values_to="Coefficient") %>%
  mutate(StepIndex=as.integer(sub("s","",Step)),Lambda=lf$lambda[StepIndex+1L],LogLambda=log(Lambda))
write.csv(path_df,file.path(panel_dir,"S1_lasso_path_corrected.csv"),row.names=FALSE)

cat("S1_INPUT_N\t",ncol(xs),"\n",sep="")
cat("S1_TOP6_SHAP\t",paste(top6,collapse=","),"\n",sep="")
cat("S1_ALL4\t",paste(Reduce(intersect,list(lasso_imp$Gene,head(rf_imp$Gene,15),boost_imp$Gene,head(xgb_imp$Gene,50))),collapse=","),"\n",sep="")

cat("STAGE\tDRAW_S1\n")
cvd <- data.frame(LogLambda=log(lcv$lambda),Mean=lcv$cvm,Lower=lcv$cvlo,Upper=lcv$cvup)
pA <- ggplot(cvd,aes(LogLambda,Mean))+geom_linerange(aes(ymin=Lower,ymax=Upper),colour=light,linewidth=.35)+
  geom_point(size=1.35,shape=21,fill=navy,colour="white",stroke=.15)+
  geom_vline(xintercept=log(lcv$lambda.min),colour=coral,linetype="dashed",linewidth=.5)+
  geom_vline(xintercept=log(lcv$lambda.1se),colour=teal,linetype="dashed",linewidth=.5)+
  labs(x=expression(log(lambda)),y="CV deviance",tag="A")+theme_pub()
pB <- ggplot(path_df,aes(LogLambda,Coefficient,group=Gene))+
  geom_line(data=filter(path_df,!Gene%in%c("S100A8","MAFG")),colour=light,linewidth=.32,alpha=.85)+
  geom_line(data=filter(path_df,Gene%in%c("S100A8","MAFG")),aes(colour=Gene),linewidth=.75)+
  geom_vline(xintercept=log(lcv$lambda.min),colour=coral,linetype="dashed",linewidth=.5)+
  geom_vline(xintercept=log(lcv$lambda.1se),colour=teal,linetype="dashed",linewidth=.5)+
  scale_colour_manual(values=c(S100A8=navy,MAFG=coral))+labs(x=expression(log(lambda)),y="Regularized coefficient",colour=NULL,tag="B")+
  theme_pub()+theme(legend.position=c(.82,.82),legend.background=element_blank())

rplot <- rank_df %>% group_by(Model) %>% arrange(desc(Value),.by_group=TRUE) %>% mutate(Rank=row_number()) %>% ungroup()
models <- c("LASSO","glmBoost","XGBoost","Ranger")
glev <- rplot %>% group_by(Gene) %>% summarise(MinRank=min(Rank),MeanRank=mean(Rank),.groups="drop") %>% arrange(MinRank,MeanRank) %>% pull(Gene)
tiles <- expand_grid(Gene=glev,Model=models) %>% left_join(rplot,by=c("Gene","Model")) %>%
  mutate(Gene=factor(Gene,levels=rev(glev)),Model=factor(Model,levels=models),Fill=ifelse(is.na(Rank),0,6-Rank))
pC <- ggplot(tiles,aes(Model,Gene))+geom_tile(aes(fill=Fill),colour="white",linewidth=.6,width=.86,height=.86)+
  geom_text(aes(label=ifelse(is.na(Rank),"",Rank),colour=ifelse(!is.na(Rank)&Rank<=2,"white",ink)),size=8/ggplot2::.pt,fontface="bold",family=font,show.legend=FALSE)+
  scale_fill_gradient(low=pale,high=teal,limits=c(0,5),guide="none")+scale_colour_identity()+
  scale_x_discrete(labels=c(Ranger="Random\nforest"))+labs(x="Within-model rank (1 = highest)",y=NULL,tag="C")+
  theme_minimal(base_size=8,base_family=font)+theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8,face="italic",colour=ink),
  panel.grid=element_blank(),plot.tag=element_text(size=12,face="bold"),plot.tag.position=c(0,1),plot.margin=ggplot2::margin(4,5,4,5))

sl <- shap_top %>% mutate(Row=row_number()) %>% pivot_longer(-Row,names_to="Gene",values_to="SHAP") %>%
  left_join(x_top%>%mutate(Row=row_number())%>%pivot_longer(-Row,names_to="Gene",values_to="FeatureValue"),by=c("Row","Gene")) %>%
  mutate(Gene=factor(Gene,levels=rev(top6)))
set.seed(123); pD <- ggplot(sl,aes(SHAP,Gene,colour=FeatureValue))+geom_vline(xintercept=0,colour=mid,linewidth=.35)+
  geom_jitter(height=.15,width=0,size=1,alpha=.82)+scale_colour_gradient2(low="#3569A8",mid="#F3F0E8",high="#CB4C64",midpoint=0)+
  labs(x="SHAP value (log-odds)",y="Feature",colour="Feature value",tag="D")+theme_pub()+theme(axis.text.y=element_text(face="italic"))

row_id <- min(2L,nrow(shap_top)); wf <- data.frame(Gene=top6,Value=as.numeric(shap_top[row_id,top6]),Feature=as.numeric(x_top[row_id,top6])) %>%
  arrange(abs(Value)) %>% mutate(Start=baseline+lag(cumsum(Value),default=0),End=Start+Value,Y=row_number(),Label=sprintf("%+.3f",Value),GeneLabel=sprintf("%s  (%.2f)",Gene,Feature))
fm <- baseline+sum(wf$Value)
pE <- ggplot(wf)+geom_rect(aes(xmin=pmin(Start,End),xmax=pmax(Start,End),ymin=Y-.31,ymax=Y+.31,fill=Value>=0),colour="white",linewidth=.3)+
  geom_text(aes(x=(Start+End)/2,y=Y,label=Label),colour="white",size=8/ggplot2::.pt,family=font)+
  geom_vline(xintercept=baseline,colour=mid,linetype="dashed",linewidth=.4)+geom_vline(xintercept=fm,colour=ink,linetype="dashed",linewidth=.4)+
  scale_fill_manual(values=c(`TRUE`=coral,`FALSE`=navy),guide="none")+scale_y_continuous(breaks=wf$Y,labels=wf$GeneLabel,expand=expansion(mult=c(.03,.12)))+
  labs(x="SHAP contribution (log-odds)",y=NULL,tag="E")+theme_pub()+theme(axis.text.y=element_text(face="italic"))

dps <- lapply(seq_along(top6),function(i){g<-top6[i]; partner<-if(g=="S100A8")setdiff(top6,g)[1] else "S100A8"; if(!partner%in%top6)partner<-setdiff(top6,g)[1]
  dd<-data.frame(X=x_top[[g]],SHAP=shap_top[[g]],Color=x_top[[partner]])
  ggplot(dd,aes(X,SHAP,colour=Color))+geom_hline(yintercept=0,colour=light,linewidth=.3)+geom_point(size=1,alpha=.82)+
    scale_colour_gradient2(low="#3569A8",mid="#F3F0E8",high="#CB4C64",midpoint=0,guide="none")+
    labs(x=g,y="SHAP value",tag=if(i==1)"F"else NULL)+annotate("text",x=Inf,y=-Inf,label=paste0("Color: ",partner),hjust=1.05,vjust=-.5,size=8/ggplot2::.pt,colour=mid,family=font)+
    theme_pub()+theme(axis.title.x=element_text(face="italic"),plot.margin=ggplot2::margin(3,4,4,4))})
pF <- wrap_plots(dps,ncol=3,nrow=2)
save_panel(pA,"S1A_lasso_cv_corrected",3.3,1.85);save_panel(pB,"S1B_lasso_path_corrected",3.3,1.85)
save_panel(pC,"S1C_model_ranks_corrected",3.3,1.9);save_panel(pD,"S1D_shap_beeswarm_corrected",3.3,1.9)
save_panel(pE,"S1E_shap_waterfall_corrected",6.7,1.35);save_panel(pF,"S1F_shap_dependence_corrected",6.7,2.3)
s1 <- ((pA|pB)/(pC|pD)/pE/pF)+plot_layout(heights=c(1.9,1.95,1.4,2.4)) & theme(plot.margin=ggplot2::margin(4,5,4,5))
ggsave(file.path(formal_dir,"S1_fig_R_master.pdf"),s1,width=7.5,height=8.75,device=cairo_pdf,family=font,bg="white")
ggsave(file.path(panel_dir,"S1_fig_corrected_preview.tif"),s1,width=7.5,height=8.75,device=ragg::agg_tiff,dpi=300,compression="lzw",bg="white")
save(lf,lcv,lasso_imp,rf,rf_imp,bf,boost_imp,xf,xgb_imp,shap_all,shap_top,x_top,rank_df,top6,baseline,overlap,
     file=file.path(panel_dir,"S1_corrected_analysis.RData"))
cat("DONE\n")
