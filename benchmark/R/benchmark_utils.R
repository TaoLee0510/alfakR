resolve_repo_dir <- function() {
  candidates <- unique(c(
    tryCatch(dirname(knitr::current_input(dir = TRUE)), error = function(e) NA_character_),
    getwd()
  ))

  for (cand in candidates) {
    if (is.na(cand) || !nzchar(cand)) {
      next
    }
    cand_norm <- normalizePath(cand, winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(cand_norm, "DESCRIPTION")) &&
        dir.exists(file.path(cand_norm, "benchmark"))) {
      return(cand_norm)
    }
    parent2 <- normalizePath(file.path(cand_norm, "..", ".."), winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(parent2, "DESCRIPTION")) &&
        dir.exists(file.path(parent2, "benchmark"))) {
      return(parent2)
    }
  }

  stop("Could not locate the alfakR repository root from the current knitting context.")
}

sort_pid_levels <- function(x) {
  x <- unique(as.character(x))
  ord_num <- suppressWarnings(as.integer(sub("^P", "", x)))
  x[order(ifelse(is.na(ord_num), Inf, ord_num), x)]
}

pm_to_label <- function(pm) {
  format(pm, scientific = FALSE, trim = TRUE)
}

write_tsv_base <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  invisible(path)
}

save_table_bundle <- function(x, stem) {
  saveRDS(x, paste0(stem, ".rds"))
  if (is.data.frame(x)) {
    write_tsv_base(x, paste0(stem, ".tsv"))
  }
  invisible(stem)
}

load_saved_table <- function(stem) {
  rds_path <- paste0(stem, ".rds")
  if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
  NULL
}

render_tbl <- function(x, caption = NULL, digits = 4) {
  if (knitr::is_html_output()) {
    return(
      knitr::kable(
        x,
        format = "html",
        digits = digits,
        caption = caption,
        table.attr = 'class="three-line-table"'
      )
    )
  }

  knitr::kable(
    x,
    format = "latex",
    digits = digits,
    caption = caption,
    booktabs = TRUE,
    longtable = TRUE
  )
}

html_escape_attr <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

emit_report_image <- function(path, alt = "") {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(invisible(FALSE))
  }

  if (knitr::is_html_output()) {
    cat(
      sprintf(
        "<div class=\"report-image\"><img src=\"%s\" alt=\"%s\" style=\"max-width:100%%; height:auto; display:block;\" /></div>\n\n",
        knitr::image_uri(path),
        html_escape_attr(alt)
      )
    )
  } else {
    cat(sprintf("![](%s)\n\n", path))
  }

  invisible(TRUE)
}

empty_note_tbl <- function(note) {
  tibble::tibble(note = note)
}

alfak_log <- function(...) {
  msg <- paste0(...)
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message("[", timestamp, "] ", msg)
  flush.console()
}

safe_read_rds <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}
