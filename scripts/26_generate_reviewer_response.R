# 输入：审稿证据矩阵。
# 输出：逐点response_to_reviewers_draft.docx，明确已完成、部分完成和数据不可用项。
rm(list=ls());options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8");source(file.path("functions","io_functions.R"),encoding="UTF-8")
with_script_log("26_generate_reviewer_response", {
  py<-"C:/Users/ZCL/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe";assert_file(py);s<-file.path("scripts","build_revision_documents.py");assert_file(s)
  code<-system2(py,c(s,"response"),stdout=file.path("11_logs","26_response_console.log"),stderr=file.path("11_logs","26_response_error.log"));if(code!=0)stop("返修信DOCX生成失败: ",code)
  assert_file("13_response_to_reviewers/response_to_reviewers_draft.docx","response DOCX")
})
