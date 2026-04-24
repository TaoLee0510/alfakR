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

safe_divide <- function(num, den) {
  num <- suppressWarnings(as.numeric(num))
  den <- suppressWarnings(as.numeric(den))
  if (!is.finite(num) || !is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  num / den
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
    return(knitr::kable(
      x,
      format = "html",
      digits = digits,
      caption = caption,
      table.attr = 'class="three-line-table"'
    ))
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

emit_report_table <- function(x, caption = NULL, digits = 4) {
  tbl <- render_tbl(x, caption = caption, digits = digits)
  if (knitr::is_html_output()) {
    cat(as.character(tbl), "\n\n", sep = "")
  } else {
    print(tbl)
  }
  invisible(TRUE)
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

emit_report_image_grid <- function(items, cols = 2L) {
  if (is.null(items) || !nrow(items)) {
    return(invisible(FALSE))
  }

  items <- tibble::as_tibble(items)
  if (!"title" %in% names(items)) {
    items$title <- rep("", nrow(items))
  }
  if (!"path" %in% names(items)) {
    stop("emit_report_image_grid() requires a `path` column.")
  }
  if (!"alt" %in% names(items)) {
    items$alt <- items$title
  }

  cols <- max(1L, as.integer(cols))

  if (knitr::is_html_output()) {
    cat(
      sprintf(
        "<div class=\"report-image-grid\" style=\"grid-template-columns: repeat(%d, minmax(0, 1fr));\">\n",
        cols
      )
    )
    for (i in seq_len(nrow(items))) {
      rr <- items[i, , drop = FALSE]
      title <- html_escape_attr(rr$title[[1]])
      path <- as.character(rr$path[[1]])
      alt <- html_escape_attr(rr$alt[[1]])

      cat("<div class=\"report-image-grid-item\">\n")
      if (nzchar(title)) {
        cat(sprintf("<div class=\"report-image-grid-title\">%s</div>\n", title))
      }
      if (!is.na(path) && nzchar(path) && file.exists(path)) {
        cat(
          sprintf(
            "<img src=\"%s\" alt=\"%s\" style=\"max-width:100%%; height:auto; display:block;\" />\n",
            knitr::image_uri(path),
            alt
          )
        )
      } else {
        cat("<div class=\"report-image-grid-missing\">Image not available.</div>\n")
      }
      cat("</div>\n")
    }
    cat("</div>\n\n")
  } else {
    for (i in seq_len(nrow(items))) {
      rr <- items[i, , drop = FALSE]
      if (nzchar(as.character(rr$title[[1]]))) {
        cat("### ", as.character(rr$title[[1]]), "\n\n", sep = "")
      }
      emit_report_image(as.character(rr$path[[1]]), alt = as.character(rr$alt[[1]]))
    }
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
