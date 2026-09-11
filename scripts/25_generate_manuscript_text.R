# 输入：最终汇总CSV。
# 输出：analysis_report.docx、final_decision_report.docx、manuscript_revision_text.docx。
rm(list=ls());options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8");source(file.path("functions","io_functions.R"),encoding="UTF-8")
with_script_log("25_generate_manuscript_text", {
  py<-"C:/Users/ZCL/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe";assert_file(py);s<-file.path("scripts","build_revision_documents.py");assert_file(s)
  code<-system2(py,c(s,"manuscript"),stdout=file.path("11_logs","25_documents_console.log"),stderr=file.path("11_logs","25_documents_error.log"));if(code!=0)stop("DOCX生成失败: ",code)
  lapply(c("12_manuscript_revision/analysis_report.docx","12_manuscript_revision/final_decision_report.docx","12_manuscript_revision/manuscript_revision_text.docx"),assert_file,"DOCX output")
})
