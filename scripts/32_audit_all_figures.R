options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(magick)
  library(pdftools)
})

root <- "D:/\u684c\u9762/sepsis"
audit_root <- file.path(Sys.getenv("TEMP"), "sepsis_all_figure_audit_r")
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

fit_tile <- function(img, label, width = 800, height = 650) {
  img <- image_background(img, "white", flatten = TRUE)
  img <- tryCatch(image_trim(img), error = function(e) img)
  info <- image_info(img)
  max_w <- width - 30
  max_h <- height - 70
  scale <- min(max_w / info$width, max_h / info$height)
  img <- image_resize(img, sprintf("%dx%d!", floor(info$width * scale),
                                   floor(info$height * scale)), filter = "Lanczos")
  img <- image_extent(img, sprintf("%dx%d", width, height - 45),
                      gravity = "center", color = "white")
  img <- image_extent(img, sprintf("%dx%d", width, height),
                      gravity = "south", color = "white")
  image_annotate(img, label, gravity = "northwest", location = "+12+10",
                 size = 22, font = "Arial", weight = 600, color = "#222222")
}

render_file <- function(path, id) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "pdf") {
    np <- pdf_info(path)$pages
    lapply(seq_len(np), function(pg) {
      im <- image_read(pdf_render_page(path, page = pg, dpi = 110))
      fit_tile(im, sprintf("%s | %s | p%d/%d", id, basename(path), pg, np))
    })
  } else {
    im <- image_read(path, density = 110)
    list(fit_tile(im[1], sprintf("%s | %s", id, basename(path))))
  }
}

make_set <- function(paths, prefix) {
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = TRUE))
  ids <- sprintf("%s%03d", prefix, seq_along(paths))
  manifest <- data.frame(id = ids, path = paths)
  write.csv(manifest, file.path(audit_root, paste0(prefix, "_manifest.csv")), row.names = FALSE)
  ims <- unlist(Map(render_file, paths, ids), recursive = FALSE)
  page_id <- ceiling(seq_along(ims) / 9)
  for (p in unique(page_id)) {
    z <- ims[page_id == p]
    while (length(z) < 9) z[[length(z) + 1L]] <- image_blank(800, 650, "white")
    rows <- lapply(split(z, rep(1:3, each = 3)), function(r) image_append(do.call(c, r)))
    sheet <- image_append(do.call(c, rows), stack = TRUE)
    image_write(sheet, file.path(audit_root, sprintf("%s_sheet_%02d.png", prefix, p)))
  }
}

stage_root <- file.path(Sys.getenv("TEMP"), "sepsis_figure_audit_stage")
staged <- function(set) sort(list.files(file.path(stage_root, set),
                                          pattern = "\\.pdf$", full.names = TRUE))
original <- staged("ORIG")
bio <- staged("BIO")
experiments <- staged("EXP")
docking_md <- staged("MD")
current <- staged("CUR")

make_set(original, "ORIG")
make_set(bio, "BIO")
make_set(experiments, "EXP")
make_set(docking_md, "MD")
make_set(current, "CUR")

cat(audit_root)
