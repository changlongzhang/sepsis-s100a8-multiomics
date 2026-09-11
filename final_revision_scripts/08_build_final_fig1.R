options(stringsAsFactors = FALSE)
.libPaths(c("C:/Users/ZCL/AppData/Local/R/win-library/4.5", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(patchwork); library(ragg); library(svglite); library(WGCNA)
})

root <- Sys.getenv("SEPSIS_PROJECT_ROOT", unset = normalizePath("..", winslash = "/", mustWork = TRUE))
wdir <- file.path(root, "work", "final_revision_20260908", "wgcna_allgenes_rerun")
out <- "C:/Users/ZCL/final_revision_outputs"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

ink <- "#20252A"; mid <- "#6F7880"; light <- "#D9DEE2"; pale <- "#F2F4F5"
red <- "#D66A5E"; blue <- "#4C78A8"; teal <- "#5B958C"; green <- "#59A14F"
ff <- "sans"
theme_pub <- function(base_size = 8) theme_classic(base_size = base_size, base_family = ff) +
  theme(text=element_text(colour=ink), axis.text=element_text(size=7.5,colour=ink),
        axis.title=element_text(size=8), axis.line=element_line(linewidth=.3,colour=ink),
        axis.ticks=element_line(linewidth=.3,colour=ink), plot.tag=element_text(size=11,face="bold"),
        plot.tag.position=c(0,1), plot.margin=margin(4,5,4,5), legend.text=element_text(size=7.5),
        legend.title=element_text(size=7.5))

deg <- read.csv(file.path(wdir, "corrected_DEG_results_used.csv"), check.names=FALSE)
mt <- read.csv(file.path(wdir, "module_trait_all_results.csv"), check.names=FALSE)
sizes <- read.csv(file.path(wdir, "module_sizes_all.csv"), check.names=FALSE)
res <- readRDS(file.path(wdir, "WGCNA_final_results.rds"))
go <- read.csv(file.path(wdir, "overlap_GO_all_results.csv"), check.names=FALSE)
kegg <- read.csv(file.path(wdir, "overlap_KEGG_all_results.csv"), check.names=FALSE)
black <- readLines(file.path(wdir, "primary_sepsis_module_black_genes.txt"), warn=FALSE)
overlap <- readLines(file.path(wdir, "primary_sepsis_module_corrected_DEG_overlap.txt"), warn=FALSE)
sig <- deg$Gene[deg$Regulation_final != "Not significant"]

pA <- ggplot(deg, aes(log2FoldChange, -log10(pmax(padj, 1e-300)))) +
  geom_point(aes(colour=Regulation_final), size=.65, alpha=.62) +
  geom_vline(xintercept=c(-1,1), colour=mid, linewidth=.3, linetype="dashed") +
  geom_hline(yintercept=-log10(.05), colour=mid, linewidth=.3, linetype="dashed") +
  geom_point(data=filter(deg,Gene=="S100A8"), shape=21, size=2.3, stroke=.6, fill=NA, colour="black") +
  geom_text_repel(data=filter(deg,Gene=="S100A8"), aes(label=Gene), size=2.7, fontface="bold",
                  min.segment.length=0, seed=123) +
  scale_colour_manual(values=c("Up"=red,"Down"=blue,"Not significant"="#C8CDD0"),
                      breaks=c("Up","Down","Not significant"),labels=c("Up","Down","NS")) +
  labs(x="Log2 fold change",y=expression(-log[10](FDR)),colour=NULL,tag="A") + theme_pub() +
  theme(legend.position=c(.82,.46),legend.background=element_blank())

mt_plot <- mt %>% mutate(module=factor(module,levels=rev(module)), label=sprintf("r = %.3f\nFDR = %.2g",correlation_sepsis,FDR))
pB <- ggplot(mt_plot,aes(x="Sepsis status",y=module,fill=correlation_sepsis)) +
  geom_tile(colour="white",linewidth=.7) + geom_text(aes(label=label),size=2.45,lineheight=.9) +
  scale_fill_gradient2(low=blue,mid="white",high=red,midpoint=0,limits=c(-.55,.55),oob=scales::squish) +
  labs(x=NULL,y="Module",fill="Pearson r",tag="B") + theme_minimal(base_size=8,base_family=ff) +
  theme(panel.grid=element_blank(),axis.text.x=element_text(size=7.5),axis.text.y=element_text(size=7.5),
        plot.tag=element_text(size=11,face="bold"),plot.tag.position=c(0,1),legend.position="right",
        plot.margin=margin(4,5,4,5))

meblack <- res$MEs[["MEblack"]]
sqc <- read.csv(file.path(wdir, "sample_qc_all_802.csv"), check.names=FALSE)
mi <- match(rownames(res$MEs), sqc$sample_id)
stopifnot(!anyNA(mi), identical(sqc$sample_id[mi], rownames(res$MEs)))
cdf <- data.frame(Group=factor(sqc$group[mi],levels=c("Normal","Sepsis")),MEblack=meblack)
wt <- wilcox.test(MEblack~Group,data=cdf,exact=FALSE)
pC <- ggplot(cdf,aes(Group,MEblack,fill=Group)) +
  geom_violin(trim=FALSE,scale="width",colour=NA,alpha=.72) +
  geom_boxplot(width=.20,outlier.shape=NA,fill="white",linewidth=.35) +
  geom_jitter(width=.10,size=.45,alpha=.28,colour=ink) +
  scale_fill_manual(values=c("Normal"=blue,"Sepsis"=red)) +
  annotate("text",x=1.5,y=max(cdf$MEblack)*1.04,label=sprintf("Wilcoxon P = %.2g",wt$p.value),size=2.6) +
  labs(x=NULL,y="Black-module eigengene",tag="C") + theme_pub() + theme(legend.position="none")

theta <- seq(0,2*pi,length.out=400)
vd <- bind_rows(data.frame(x=.40+.27*cos(theta),y=.50+.31*sin(theta),set="WGCNA black"),
                data.frame(x=.60+.27*cos(theta),y=.50+.31*sin(theta),set="DEGs"))
pD <- ggplot() + geom_polygon(data=vd,aes(x,y,group=set,fill=set),alpha=.78,colour="black",linewidth=.45) +
  annotate("text",x=.23,y=.50,label=length(setdiff(black,sig)),fontface="bold",size=3.5) +
  annotate("text",x=.50,y=.50,label=length(overlap),fontface="bold",size=3.5) +
  annotate("text",x=.77,y=.50,label=length(setdiff(sig,black)),fontface="bold",size=3.5) +
  annotate("text",x=.24,y=.85,label="WGCNA black",size=2.7) + annotate("text",x=.76,y=.85,label="DEGs",size=2.7) +
  scale_fill_manual(values=c("WGCNA black"=red,"DEGs"=blue)) + coord_fixed(xlim=c(.03,.97),ylim=c(.08,.93),clip="off") +
  labs(tag="D") + theme_void(base_family=ff) + theme(legend.position="none",plot.tag=element_text(size=11,face="bold"),plot.tag.position=c(0,1))

go_top <- go %>% filter(p.adjust<.05) %>% group_by(ONTOLOGY) %>% arrange(p.adjust,.by_group=TRUE) %>%
  slice_head(n=4) %>% ungroup() %>% mutate(lbl=stringr::str_wrap(Description,38)) %>% arrange(ONTOLOGY,Count) %>%
  mutate(lbl=factor(lbl,levels=unique(lbl)))
pE <- ggplot(go_top,aes(Count,lbl,fill=ONTOLOGY)) + geom_col(width=.68) +
  scale_fill_manual(values=c("BP"=red,"CC"=blue,"MF"=green)) +
  labs(x="Gene count",y=NULL,fill=NULL,tag="E") + theme_pub() +
  theme(axis.text.y=element_text(size=6.5),axis.ticks.y=element_blank(),legend.position="bottom",
        legend.direction="horizontal",legend.key.width=unit(3,"mm"),legend.key.height=unit(2.5,"mm"))

kegg_top <- kegg %>% arrange(p.adjust,pvalue) %>% slice_head(n=5) %>%
  mutate(lbl=stringr::str_wrap(Description,28)) %>% arrange(Count) %>% mutate(lbl=factor(lbl,levels=unique(lbl)))
pF <- ggplot(kegg_top,aes(Count,lbl)) + geom_col(width=.65,fill=teal) +
  labs(x="Gene count",y=NULL,tag="F",caption="Nominal terms shown;\nnone passed FDR < 0.05") + theme_pub() +
  theme(axis.text.y=element_text(size=6.5),axis.ticks.y=element_blank(),plot.caption=element_text(size=7,hjust=0,colour=mid),
        plot.margin=margin(4,10,4,5))

fig <- (pA | pB | pC) / (pD | pE | pF) + plot_layout(widths=c(1.05,1.25,1),heights=c(1,1))
ggsave(file.path(out,"Fig1.pdf"),fig,width=7.5,height=6.0,device=cairo_pdf,bg="white")
ggsave(file.path(out,"Fig1.svg"),fig,width=7.5,height=6.0,device=svglite,bg="white")
ggsave(file.path(out,"Fig1.tif"),fig,width=7.5,height=6.0,device=ragg::agg_tiff,dpi=600,compression="lzw",bg="white")
write.csv(data.frame(test="Wilcoxon rank-sum",comparison="Black-module eigengene: Sepsis vs Normal",p_value=wt$p.value),
          file.path(out,"Fig1C_statistics.csv"),row.names=FALSE)
cat(sprintf("Fig1 completed; Wilcoxon P=%.8g\n",wt$p.value))
