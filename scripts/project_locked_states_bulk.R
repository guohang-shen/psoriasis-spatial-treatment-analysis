## Project the locked single-cell state modules into independent bulk clinical
## cohorts. This is an orthogonal pharmacodynamic check, not cell-level TNF
## localization and not outcome-driven feature selection.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ .libPaths(c("C:/codex_r_lib_psoriasis", .libPaths())) })
suppressPackageStartupMessages({ library(GEOquery); library(Biobase) })
root <- "."; data_dir <- "data"; res_dir <- "results"
modules <- list(
  KC_IL17_IL36 = c("KRT16", "KRT17", "KRT6A", "KRT6B", "S100A7", "S100A8", "S100A9", "DEFB4A", "IL36G", "CCL20", "SERPINB3", "SERPINB4"),
  KC_differentiation = c("FLG", "LOR", "IVL", "KRT1", "KRT10", "KRT2", "KRT3"),
  T17_like = c("CCR6", "RORA", "KLRB1", "IL7R", "IL23R", "IL32", "CD3D", "TRBC1"),
  T_resident = c("CD69", "ITGAE", "CXCR6", "CD3D", "TRBC1", "IL7R"),
  IFN_response = c("ISG15", "IFIT1", "IFIT3", "IFITM1", "MX1", "OAS1", "OASL", "STAT1"),
  myeloid_inflammatory = c("S100A8", "S100A9", "IL1B", "CXCL8", "TNF", "FCGR3A", "NLRP3", "CTSD"),
  APC_DC = c("HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "FCER1A", "CD1C", "CLEC10A"),
  fibroblast_inflammatory = c("CCL19", "CCL13", "CXCL14", "CXCL12", "FAP", "C7", "IL6", "COL3A1"))
get_meta <- function(pd, key) {
  cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
  apply(pd[, cols, drop = FALSE], 1, function(row) { hit <- grep(paste0("^", key, ":\\s*"), row, value = TRUE, ignore.case = TRUE); if (!length(hit)) NA_character_ else sub("^[^:]+:\\s*", "", hit[[1]]) })
}
get_expr <- function(eset, platform, genes) {
  mat <- exprs(eset)
  gpl <- getGEO(platform, destdir = data_dir, AnnotGPL = TRUE)
  tab <- Table(gpl)
  sym_col <- grep("gene symbol|symbol", colnames(tab), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sym_col)) stop("No Gene Symbol column in ", platform)
  mapped <- data.frame(PROBEID = as.character(tab$ID), SYMBOL = as.character(tab[[sym_col]]), stringsAsFactors = FALSE)
  mapped <- mapped[mapped$PROBEID %in% rownames(mat) & !is.na(mapped$SYMBOL), , drop = FALSE]
  out <- matrix(NA_real_, length(genes), ncol(mat), dimnames = list(genes, colnames(mat)))
  for (g in genes) {
    probes <- mapped$PROBEID[vapply(mapped$SYMBOL, function(z) g %in% trimws(unlist(strsplit(z, "///", fixed = TRUE))), logical(1))]
    probes <- intersect(probes, rownames(mat)); if (length(probes)) out[g, ] <- apply(mat[probes, , drop = FALSE], 2, median, na.rm = TRUE)
  }
  out
}
score <- function(x, genes) { present <- genes[rowSums(!is.na(x[genes, , drop = FALSE])) > 0]; if (length(present) < max(3, ceiling(.5 * length(genes)))) return(rep(NA_real_, ncol(x))); colMeans(x[present, , drop = FALSE], na.rm = TRUE) }
analyse <- function(accession, platform) {
  e <- getGEO(filename = file.path(data_dir, paste0(accession, "_series_matrix.txt.gz")), getGPL = FALSE)
  pd <- pData(e)
  if (accession == "GSE117239") {
    m <- data.frame(gsm = pd$geo_accession, subject = get_meta(pd, "subjectid"), tissue = get_meta(pd, "tissue"), week = get_meta(pd, "timepoint"), drug = get_meta(pd, "treatment"), response = get_meta(pd, "pasi75 response at week 12"), stringsAsFactors = FALSE)
  } else {
    m <- data.frame(gsm = pd$geo_accession, subject = sub("^(Subject [0-9]+).*", "\\1", pd$title), week = get_meta(pd, "timepoint"), drug = get_meta(pd, "treatment"), pasi = suppressWarnings(as.numeric(get_meta(pd, "pasi"))), stringsAsFactors = FALSE)
  }
  x <- get_expr(e, platform, unique(unlist(modules)))
  cov <- data.frame(accession = accession, module = names(modules), n_genes = vapply(modules, function(z) sum(z %in% rownames(x)[rowSums(!is.na(x)) > 0]), integer(1)), total_genes = lengths(modules), stringsAsFactors = FALSE)
  write.csv(cov, file.path(res_dir, paste0(accession, "_locked_state_bulk_coverage.csv")), row.names = FALSE)
  for (nm in names(modules)) m[[nm]] <- score(x, modules[[nm]])
  out <- list(); k <- 0L
  for (p in unique(m$subject)) {
    v <- m[m$subject == p, , drop = FALSE]
    if (accession == "GSE117239") { a <- which(v$tissue == "Lesional" & v$week == "Baseline"); b <- which(v$tissue == "Lesional" & v$week == "Week 1"); if (!length(a) || !length(b)) next; a <- a[1]; b <- b[1]; resp <- v$response[a] } else { a <- which(v$week == "LS"); b <- which(v$week == "WK 1"); z <- which(v$week == "WK 16"); if (!length(a) || !length(b) || !length(z)) next; a <- a[1]; b <- b[1]; z <- z[1]; if (!is.finite(v$pasi[a]) || !is.finite(v$pasi[z])) next; resp <- ifelse(v$pasi[z] / pmax(v$pasi[a], 1e-8) <= .25, "R", "NR") }
    rr <- data.frame(accession = accession, subject = p, drug = v$drug[a], response = resp, stringsAsFactors = FALSE)
    for (nm in names(modules)) {
      rr[[paste0(nm, "_baseline")]] <- v[[nm]][a]
      rr[[paste0(nm, "_delta")]] <- v[[nm]][a] - v[[nm]][b]
      ## GSE85034 also has a prespecified week-16 biopsy. Keep the week-1
      ## delta above for comparability with GSE117239, and store week-16
      ## separately so late pharmacodynamic recovery is not mixed with it.
      rr[[paste0(nm, "_week16_delta")]] <- if (accession == "GSE85034") v[[nm]][a] - v[[nm]][z] else NA_real_
    }
    k <- k + 1L; out[[k]] <- rr
  }
  do.call(rbind, out)
}
a <- analyse("GSE117239", "GPL570"); b <- analyse("GSE85034", "GPL10558"); all_scores <- rbind(a, b)
write.csv(all_scores, file.path(res_dir, "bulk_locked_state_patient_scores.csv"), row.names = FALSE)
delta_cols <- paste0(rep(names(modules), each = 1), "_delta")
sum_rows <- list(); k <- 0L
for (acc in unique(all_scores$accession)) for (drug in unique(all_scores$drug[all_scores$accession == acc])) for (nm in names(modules)) {
  d <- all_scores[all_scores$accession == acc & all_scores$drug == drug, , drop = FALSE]; y <- d[[paste0(nm, "_delta")]]; ok <- is.finite(y); y <- y[ok]
  if (length(y) < 3) next
  k <- k + 1L; sum_rows[[k]] <- data.frame(accession = acc, drug = drug, module = nm, n = length(y), median_delta = median(y), positive_fraction = mean(y > 0), stringsAsFactors = FALSE)
}
summary <- do.call(rbind, sum_rows); write.csv(summary, file.path(res_dir, "bulk_locked_state_drug_summary.csv"), row.names = FALSE)
print(summary[order(summary$accession, summary$drug, summary$module), ], row.names = FALSE)

