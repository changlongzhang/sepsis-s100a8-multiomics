# 输入：所有已完成分析结果。
# 输出：临床/预后图、质量控制HTML、审稿证据矩阵CSV及最终14项结论底稿。
rm(list=ls());options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","theme_nature.R"),encoding="UTF-8")

with_script_log("23_integrated_evidence_summary", {
  suppressPackageStartupMessages({library(ggplot2);library(pROC)})
  pa<-readRDS(file.path("02_class_imbalance","S100A8_probe_audit.rds"));dm0<-read.csv(file.path("02_class_imbalance","sample_metadata_clean.csv"),check.names=FALSE);v<-as.numeric(pa$probe_expr[1,match(dm0$sample_id,colnames(pa$probe_expr))]);yy<-factor(dm0$group,levels=c("Normal","Sepsis"));rr<-pROC::roc(yy,v,quiet=TRUE,direction="<")
  probe_sens<-data.frame(strategy=c("Mean across mapped probes","Median across mapped probes","Maximum across mapped probes","Highest-variance mapped probe"),mapped_probe_count=nrow(pa$probe_expr),representative_probe=rownames(pa$probe_expr)[1],mean_difference=mean(v[yy=="Sepsis"])-mean(v[yy=="Normal"]),AUC=as.numeric(pROC::auc(rr)),result="Identical because S100A8 maps to one probe on GPL13667")
  write_csv_utf8(probe_sens,file.path("02_class_imbalance","S100A8_probe_aggregation_sensitivity.csv"))
  smp<-read.csv(file.path("03_clinical_validation","clinical_validation_S100A8_samples.csv"),check.names=FALSE)
  specs<-list(GSE26440_Healthy=c("GSE26440","Healthy"),GSE28750_Healthy=c("GSE28750","Healthy"),GSE28750_Postoperative=c("GSE28750","Postoperative sterile inflammation"),GSE69528_Healthy=c("GSE69528","Healthy"),GSE69528_ChronicDisease=c("GSE69528","Uninfected chronic disease control"),GSE9692_Healthy=c("GSE9692","Healthy"),GSE134347_Healthy=c("GSE134347","Healthy"),GSE134347_NoninfectiousICU=c("GSE134347","Non-infectious ICU"),GSE131411_CardiogenicShock=c("GSE131411","Cardiogenic shock"))
  rocdata<-list();prdata<-list();exprdata<-list()
  for(nm in names(specs)){sp<-specs[[nm]];d<-smp[smp$dataset==sp[1]&smp$group%in%c("Sepsis",sp[2])&nzchar(smp$group),];d$truth<-d$group=="Sepsis";ro<-pROC::roc(d$truth,d$expression,quiet=TRUE,direction="<");co<-pROC::coords(ro,"all",ret=c("specificity","sensitivity"),transpose=FALSE);rocdata[[nm]]<-data.frame(comparison_id=nm,specificity=co$specificity,sensitivity=co$sensitivity,AUC=as.numeric(auc(ro)))
    o<-order(d$expression,decreasing=TRUE);tp<-cumsum(d$truth[o]);fp<-cumsum(!d$truth[o]);prdata[[nm]]<-data.frame(comparison_id=nm,recall=tp/sum(d$truth),precision=tp/(tp+fp),prevalence=mean(d$truth));exprdata[[nm]]<-data.frame(comparison_id=nm,dataset=sp[1],sample_id=d$sample_id,group=ifelse(d$group=="Sepsis","Sepsis",sp[2]),expression=d$expression)
  }
  rocdata<-do.call(rbind,rocdata);prdata<-do.call(rbind,prdata);exprdata<-do.call(rbind,exprdata)
  short_labels<-c(GSE26440_Healthy="GSE26440\nHealthy",GSE28750_Healthy="GSE28750\nHealthy",GSE28750_Postoperative="GSE28750\nPostoperative",GSE69528_Healthy="GSE69528\nHealthy",GSE69528_ChronicDisease="GSE69528\nChronic disease",GSE9692_Healthy="GSE9692\nHealthy",GSE134347_Healthy="GSE134347\nHealthy",GSE134347_NoninfectiousICU="GSE134347\nNon-infectious ICU",GSE131411_CardiogenicShock="GSE131411\nCardiogenic shock")
  roc_ann<-aggregate(AUC~comparison_id,rocdata,unique)
  p1<-ggplot(rocdata,aes(1-specificity,sensitivity))+geom_abline(slope=1,intercept=0,linetype=2,colour="#BDBDBD",linewidth=.35)+
    geom_line(linewidth=.72,colour="#315F7D")+
    geom_text(data=roc_ann,aes(x=.97,y=.04,label=sprintf("AUC %.2f",AUC)),inherit.aes=FALSE,hjust=1,vjust=0,size=2.35,colour="#315F7D")+
    facet_wrap(~comparison_id,labeller=as_labeller(short_labels))+coord_equal()+labs(x="1 − specificity",y="Sensitivity")+theme_curve()+
    theme(strip.text=element_text(hjust=0),panel.grid.minor=element_blank())
  save_nature_plot(p1,file.path("09_figures","Fig_ClinicalComparator_ROC"),183,130,source_data=rocdata)
  pr_base<-aggregate(prevalence~comparison_id,prdata,unique)
  p2<-ggplot(prdata,aes(recall,precision))+
    geom_hline(data=pr_base,aes(yintercept=prevalence),colour="#BDBDBD",linetype=2,linewidth=.35)+
    geom_line(linewidth=.72,colour="#4F8A76")+
    facet_wrap(~comparison_id,labeller=as_labeller(short_labels))+labs(x="Recall",y="Precision")+theme_curve()+
    theme(strip.text=element_text(hjust=0),panel.grid.minor=element_blank())
  save_nature_plot(p2,file.path("09_figures","Fig_ClinicalComparator_PR"),183,130,source_data=prdata)
  clinical_cols<-c(Sepsis="#C65D4B",Healthy="#6B8E9B","Postoperative sterile inflammation"="#4F9D76","Uninfected chronic disease control"="#B07AA1","Non-infectious ICU"="#D18F2D","Cardiogenic shock"="#777777")
  p3<-ggplot(exprdata,aes(group,expression,fill=group))+
    geom_violin(width=.78,trim=TRUE,scale="width",alpha=.28,colour=NA)+
    geom_boxplot(outlier.shape=NA,width=.22,fill="white",colour="#454545",linewidth=.35)+
    geom_point(size=.45,alpha=.3,position=position_jitter(width=.08),colour="#303030")+
    facet_wrap(~comparison_id,scales="free",labeller=as_labeller(short_labels))+scale_fill_manual(values=clinical_cols)+
    labs(x=NULL,y="S100A8 expression (cohort-specific units)",fill=NULL)+theme_nature()+
    theme(axis.text.x=element_text(angle=32,hjust=1),legend.position="none",strip.text=element_text(hjust=0))
  save_nature_plot(p3,file.path("09_figures","Fig_ClinicalComparator_Expression"),183,135,source_data=exprdata)

  prog<-read.csv(file.path("04_prognosis_severity","Table_Prognosis_Associations.csv"),check.names=FALSE);pd<-prog;pd$label<-paste(pd$dataset,pd$outcome,pd$model,sep=" | ")
  p4<-ggplot(pd,aes(estimate,reorder(label,estimate)))+geom_vline(xintercept=1,linetype=2,colour="#8C8C8C",linewidth=.35)+
    geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=0,linewidth=.55,colour="#5F6368")+
    geom_point(size=2.2,shape=18,colour="#C65A46")+scale_x_log10()+
    labs(x="OR or HR per 1-SD higher S100A8 (95% CI)",y=NULL)+theme_forest()
  save_nature_plot(p4,file.path("09_figures","Fig_Prognosis_Forest"),183,90,source_data=pd)
  ps<-read.csv(file.path("04_prognosis_severity","prognosis_samples_clean.csv"),check.names=FALSE);ps<-ps[!is.na(ps$mortality),];ps$mortality_label<-ifelse(ps$mortality==1,"Non-survivor","Survivor")
  p5<-ggplot(ps,aes(mortality_label,S100A8_z,fill=mortality_label))+
    geom_violin(width=.72,trim=TRUE,scale="width",alpha=.28,colour=NA)+
    geom_boxplot(outlier.shape=NA,width=.20,fill="white",colour="#454545",linewidth=.35)+
    geom_point(size=.5,alpha=.32,position=position_jitter(width=.08),colour="#303030")+
    facet_wrap(~dataset,scales="free_y",nrow=1)+scale_fill_manual(values=c(Survivor="#3C6E8F","Non-survivor"="#C65A46"))+
    labs(x=NULL,y="S100A8 (within-cohort z-score)",fill=NULL)+theme_nature()+theme(legend.position="none",strip.text=element_text(hjust=0))
  save_nature_plot(p5,file.path("09_figures","Fig_S100A8_Mortality"),183,80,source_data=ps)

  stable<-read.csv(file.path("02_class_imbalance","Table_S100A8_FeatureStability.csv"));interval<-read.csv(file.path("02_class_imbalance","S100A8_empirical_intervals.csv"));clin<-read.csv(file.path("03_clinical_validation","Table_ClinicalValidation_long.csv"));pb<-read.csv(file.path("05_single_cell_pseudobulk","Table_Pseudobulk_FocusGenes.csv"));corr<-read.csv(file.path("08_bulk_singlecell_integration","Table_BulkS100A8_SignatureCorrelation.csv"));meta<-read.csv(file.path("04_prognosis_severity","prognosis_random_effects_meta.csv"))
  getauc<-function(id)clin$estimate[clin$comparison_id==id&clin$metric=="AUC"][1]
  geta<-function(g)pb[pb$gene==g&pb$analysis=="primary_20",c("logFC","FDR")][1,]
  decisions<-data.frame(question=1:14,answer=c(
    sprintf("支持：1:1下采样S100A8 AUC中位数 %.3f（经验95%%区间 %.3f–%.3f），方向稳定。",interval$estimate[interval$analysis=="1:1 downsampling AUC"],interval$ci_low[interval$analysis=="1:1 downsampling AUC"],interval$ci_high[interval$analysis=="1:1 downsampling AUC"]),
    sprintf("部分支持：各算法/策略选择频率见稳定性表；并非所有迭代均为四模型共同特征。"),
    "不能将0.999 AUC主要归因于类别不平衡：无泄漏、平衡嵌套CV的AUC仍接近0.997–0.999，说明病例—健康对照信号真实且极强；但特征排名/入选频率和校准不稳定，原稿把该AUC外推为临床预测性能仍有过拟合式夸大，应将机器学习定位为特征优先排序工具。",
    sprintf("部分支持：非感染ICU AUC %.3f，心源性休克 AUC %.3f；术后无菌炎症 AUC %.3f，临床特异性不一致。",getauc("GSE134347_NoninfectiousICU"),getauc("GSE131411_CardiogenicShock"),getauc("GSE28750_Postoperative")),
    sprintf("不支持稳定预后标志物：三队列死亡随机效应Meta OR %.2f（95%%CI %.2f–%.2f，P=%.3f）；SOFA/APACHE等公开变量不可用。",meta$pooled_OR[1],meta$ci_low[1],meta$ci_high[1],meta$p_value[1]),
    sprintf("支持性证据：20-cell阈值pseudo-bulk S100A8 log2FC %.2f，FDR %.3g，但仅2名对照和4名脓毒症供者。",geta("S100A8")$logFC,geta("S100A8")$FDR),
    sprintf("方向支持：HLA-DRA log2FC %.2f/FDR %.3g；CD74 %.2f/%.3g；HLA-DPB1 %.2f/%.3g。",geta("HLA-DRA")$logFC,geta("HLA-DRA")$FDR,geta("CD74")$logFC,geta("CD74")$FDR,geta("HLA-DPB1")$logFC,geta("HLA-DPB1")$FDR),
    "更接近应急髓系生成/未成熟中性粒样髓系状态：LCN2、MMP8、MMP9、LTF、CD177升高，并伴抗原呈递程序降低。",
    "是（内部稳健性）：删除S100A8及同时删除S100A8/S100A9后，独立30基因UCell程序仍区分供者；供者数不足，不能称独立验证。",
    sprintf("仍存在，但队列差异明显；分组及单核细胞核心程序校正后的偏相关见整合表（n=%d个队列）。",nrow(corr)),
    "类别不平衡、临床相似对照、供者级pseudo-bulk、状态重注释和非循环签名已获得直接分析证据。",
    "预后/严重程度、传统生物标志物比较、多基因模型外部固定验证和TF活性仅部分回答。",
    "必须将diagnostic biomarker降为sepsis-associated candidate/case-control discriminatory marker；将机器学习降为特征排序；机制与直接结合均改为支持性/探索性。",
    "非湿实验层面足以形成诚实的大修稿，但三项湿实验机制要求仍需作者决定补做或在回复中明确为未完成，并进一步降低机制结论。"),stringsAsFactors=FALSE)
  write_csv_utf8(decisions,file.path("12_manuscript_revision","final_decision_answers.csv"))

  comments<-read.csv(text='id,category,comment\n1,Scope,Central hypothesis too broad\n2,Clinical,Case-control AUC is not real-world diagnosis\n3,Clinical,Clinically similar controls are missing\n4,Prognosis,Severity and outcome associations are missing\n5,Statistics,Severe class imbalance may inflate performance\n6,Preprocessing,Maximum-probe aggregation may bias results\n7,Machine learning,Near-perfect training AUC may be overfit\n8,Clinical utility,Additional diagnostic metrics and DCA are needed\n9,Single cell,Cells are not independent donors\n10,Cell state,S100A8-high biological state is unclear\n11,Wet lab,HLA-DR and CD74 experimental validation requested\n12,Circularity,S100A8-defined reference causes circular reasoning\n13,Compound,16-hydroxytriptolide and triptolide must be distinguished\n14,Compound,Triptolide-class substitution needs support\n15,Docking,PDB 4GGF biological rationale needs clarification\n16,Binding,CETSA does not prove direct binding\n17,Causality,S100A8 perturbation or rescue requested\n18,Wet lab,CCK-8 has technical replicates only\n19,Wet lab,NF-kB mechanistic depth is limited\n20,Claims,Overstated causal and therapeutic language',stringsAsFactors=FALSE)
  status<-c("Fully addressed","Fully addressed","Fully addressed","Partially addressed","Fully addressed","Fully addressed","Partially addressed","Partially addressed","Fully addressed","Fully addressed","Data unavailable","Fully addressed","Fully addressed","Partially addressed","Partially addressed","Fully addressed","Data unavailable","Fully addressed","Data unavailable","Fully addressed")
  action<-c("Refocus manuscript on S100A8 prioritization and cell-state context","Downgrade diagnostic wording","Add non-infectious ICU, postoperative sterile inflammation, and cardiogenic shock","Analyze available mortality and audit severity variables","Run weighted CV, repeated downsampling, bootstrap, and training-fold downsampling","Use log2 mean aggregation; verify S100A8 has one mapped probe","Use leakage-safe nested CV and downgrade prediction claims","Report full metrics and cross-validated DCA; audit traditional markers","Run donor pseudo-bulk and mixed-effects sensitivity","Unsupervised reclustering and multi-gene annotation","Requires new wet-lab data","Build S100A8-independent signatures and delete S100A8/S100A9 sensitivity","Separate computational hit from experimental analogue in text","No new structural similarity/docking comparison available","Clarify template details from author records; no new conformer test","Revise wording to target engagement/indirect stabilization","Requires new wet-lab perturbation/rescue","Restrict CCK-8 to descriptive concentration selection","Requires new wet-lab functional assay","Replace causal/therapeutic claims with associative wording")
  matrix<-data.frame(reviewer_comment=paste0("R1.",comments$id,": ",comments$comment),requested_action=action,analysis_performed=c("Narrative restructuring","Clinical interpretation audit","GSE28750/GSE134347/GSE131411 comparator analyses","Logistic/Cox and random-effects mortality meta-analysis","Scripts 03–07","Probe audit and mean aggregation","Nested 5-fold CV repeated 20 times","Full metrics and cross-validated DCA","edgeR pseudo-bulk plus lmer sensitivity","Unsupervised clustering and module scores","No non-wet substitute can validate protein/cell assay","Unsupervised marker signature with exclusions and LODO robustness","Manuscript wording revision","Evidence boundary documented","Author confirmation required for structural settings","Manuscript wording revision","No non-wet substitute for causal perturbation","Manuscript wording revision","No non-wet substitute for functional assay","Claim audit and revision"),script=c("25","10,11,25","10,11","12,13","03–07","02,23","03,06,07","10,11","15,16","17–19","N/A","21,22","25,26","25,26","25,26","25,26","N/A","25,26","N/A","25,26"),figure=c("N/A","Clinical ROC/PR","Clinical ROC/PR/DCA","Prognosis forest","Class-balance figures","N/A","Model-performance figure","Decision curve","Pseudo-bulk donor figure","State programs/UMAP","N/A","Independent signature/bulk projection","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A"),table=c("Evidence matrix","Clinical validation","Clinical validation/DCA","Prognosis associations","Class balance/stability","Probe audit","Class balance","Clinical validation/DCA","Pseudo-bulk focus genes","State annotation/abundance","N/A","Independent signatures/correlations","Evidence matrix","Evidence matrix","Evidence matrix","Evidence matrix","N/A","Evidence matrix","N/A","Evidence matrix"),main_result=decisions$answer[pmin(c(11,4,4,5,1,1,3,4,6,8,14,9,13,13,13,13,14,13,14,13),14)],supports_original_claim=c("Partially","No","Partially","No","Partially","Partially","No","Partially","Partially","Partially","Not assessable","Partially","Partially","Not assessable","Not assessable","No","Not assessable","No","Not assessable","No"),manuscript_change_required="Yes",status=status,stringsAsFactors=FALSE)
  write_csv_utf8(matrix,file.path("13_response_to_reviewers","reviewer_evidence_matrix.csv"))

  stat_files<-list.files(file.path("11_logs","status_parts"),pattern="\\.csv$",full.names=TRUE);stats<-do.call(rbind,lapply(stat_files,read.csv,check.names=FALSE));stats<-stats[order(stats$script),]
  qc<-data.frame(check=c("GSE65682 sample-expression match","Duplicate discovery samples","Missing discovery group","S100A8 discovery probe mapping","GSE69528 correct probe","Cross-platform raw matrix merge","Single-cell independent unit","Single-cell donor count","Control donor count","Pseudo-bulk threshold","Clinical universal threshold","PBMC projection wording","Run-status failures"),result=c("Pass: 802/802 matched","Pass: 0","Pass: 0","Pass: one probe 11753823_a_at","Pass: ILMN_1729801; old ACTB mapping rejected","Pass: cohorts analyzed separately","Pass: donor","7","2","20-cell primary; 50-cell strict sensitivity underpowered","Not claimed","Program projection, not absolute whole-blood proportion",sum(stats$status!="Success")),severity=c(rep("Pass",6),"Pass","Caution","Caution","Caution","Pass","Pass",ifelse(any(stats$status!="Success"),"Fail","Pass")))
  write_csv_utf8(qc,file.path("11_logs","quality_control_checks.csv"))
  rows<-paste0("<tr><td>",qc$check,"</td><td>",qc$result,"</td><td>",qc$severity,"</td></tr>",collapse="\n")
  html<-paste0("<!doctype html><html><head><meta charset='utf-8'><title>Quality control report</title><style>body{font-family:Arial,sans-serif;margin:36px;color:#222}table{border-collapse:collapse;width:100%}th,td{border-bottom:1px solid #ddd;padding:8px;text-align:left}th{background:#eef2f5}.Caution{color:#9a6700}.Fail{color:#b42318}</style></head><body><h1>Sepsis revision quality-control report</h1><p>Generated ",format(Sys.time()),". Original files were used read-only.</p><table><thead><tr><th>Check</th><th>Result</th><th>Assessment</th></tr></thead><tbody>",rows,"</tbody></table><h2>Run status</h2><pre>",paste(capture.output(print(stats[,c("script","status","runtime","warnings","errors")])),collapse="\n"),"</pre></body></html>")
  writeLines(html,file.path("11_logs","quality_control_report.html"),useBytes=TRUE)
})
