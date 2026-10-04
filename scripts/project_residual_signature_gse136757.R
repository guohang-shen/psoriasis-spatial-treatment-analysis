## Fourth drug class: TYK2/JAK1 inhibitor PF-06700841 (GSE136757, GPL570).
## Project the LOCKED residual-after-clearance signature (100 up / 100 down)
## onto a placebo-controlled, longitudinal, lesional/non-lesional biopsy set.
## Biological unit: patient. Hypothesis: a higher residual programme
## (incomplete molecular resolution) tracks with incomplete clinical
## resolution (smaller PASI improvement). No gene is reselected here.
options(stringsAsFactors = FALSE)
.libPaths(c("C:/codex_r_lib_psoriasis", .libPaths()))
suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(data.table)
})
base_dir <- "C:/Users/15082/Desktop/Psoriasis"; res_dir <- file.path(base_dir, "results")

## ---- locked residual signature ----
sig <- read.csv(file.path(res_dir, "residual_after_clearance_signature_genes.csv"))
up_g <- unique(sig$symbol[sig$direction == "up_in_open"])
dn_g <- unique(sig$symbol[sig$direction == "down_in_open"])
all_g <- unique(c(up_g, dn_g))

## ---- expression + metadata ----
eset <- getGEO(filename = file.path(base_dir, "data", "GSE136757",
                                    "GSE136757_series_matrix.txt.gz"), getGPL = FALSE)
pd <- pData(eset); expr <- exprs(eset)
stopifnot(identical(colnames(expr), rownames(pd)))
get_field <- function(key) {
  cc <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
  apply(pd[, cc, drop = FALSE], 1, function(z) {
    hit <- grep(paste0("^", key, ":\\s*"), z, value = TRUE, ignore.case = TRUE)
    if (length(hit)) sub("^[^:]+:\\s*", "", hit[1]) else NA_character_
  })
}
meta <- data.frame(gsm = pd$geo_accession, patient = get_field("subjid"),
                   tissue = get_field("tissue"), day = get_field("time"),
                   arm = get_field("treatment"),
                   pasi = suppressWarnings(as.numeric(get_field("pasi_total"))),
                   stringsAsFactors = FALSE)
stopifnot(!anyNA(meta[, c("gsm", "patient", "tissue", "day", "arm")]))

## ---- GPL570 probe -> symbol (cached annot file, no network) ----
annot_file <- file.path(base_dir, "data", "GPL570.annot.gz")
con <- gzfile(annot_file, "rt"); al <- readLines(con, warn = FALSE); close(con)
i0 <- grep("^!platform_table_begin", al); i1 <- grep("^!platform_table_end", al)
tab <- fread(text = paste(al[i0:(i1 - 1)], collapse = "\n"),
             sep = "\t", quote = "", showProgress = FALSE)
sym_col <- grep("gene symbol", colnames(tab), ignore.case = TRUE, value = TRUE)[1]
probe_ids <- as.character(tab$ID); symbols <- as.character(tab[[sym_col]])
valid <- probe_ids %in% rownames(expr) & !is.na(symbols)
sp <- strsplit(symbols[valid], "///", fixed = TRUE)
smap <- data.frame(probe = rep(probe_ids[valid], lengths(sp)),
                   symbol = trimws(unlist(sp, use.names = FALSE)), stringsAsFactors = FALSE)

## Collapse probes to a single gene value (median of probes), then restrict to
## signature genes. Genes absent from GPL570 are reported, never imputed.
gene_expr <- matrix(NA_real_, length(all_g), ncol(expr),
                    dimnames = list(all_g, colnames(expr)))
n_probe <- setNames(integer(length(all_g)), all_g)
for (g in all_g) {
  p <- intersect(smap$probe[smap$symbol == g], rownames(expr))
  n_probe[g] <- length(p)
  if (length(p)) gene_expr[g, ] <- apply(expr[p, , drop = FALSE], 2, median, na.rm = TRUE)
}
present <- all_g[rowSums(is.finite(gene_expr)) == ncol(expr)]
cat("signature genes mapped:", length(present), "/", length(all_g),
    " (up", sum(present %in% up_g), "/", length(up_g),
    "; down", sum(present %in% dn_g), "/", length(dn_g), ")\n")
write.csv(data.frame(gene = all_g, arm = ifelse(all_g %in% up_g, "up_in_open", "down_in_open"),
                     n_probes = as.integer(n_probe[all_g]), mapped = all_g %in% present),
          file.path(res_dir, "GSE136757_residual_signature_gene_coverage.csv"), row.names = FALSE)

## Per-gene z-score within this cohort, then arm means. Score = up mean - down
## mean, i.e. larger = more residual programme retained.
gz <- t(scale(t(gene_expr[present, , drop = FALSE])))
up_p <- intersect(present, up_g); dn_p <- intersect(present, dn_g)
score <- colMeans(gz[up_p, , drop = FALSE], na.rm = TRUE) -
         colMeans(gz[dn_p, , drop = FALSE], na.rm = TRUE)
meta$residual_score <- as.numeric(score)
meta$residual_up <- as.numeric(colMeans(gz[up_p, , drop = FALSE], na.rm = TRUE))
meta$residual_down <- as.numeric(colMeans(gz[dn_p, , drop = FALSE], na.rm = TRUE))
write.csv(meta, file.path(res_dir, "GSE136757_residual_signature_sample_scores.csv"), row.names = FALSE)

## ---- patient-level pairing (baseline LS -> post LS) ----
rows <- list()
for (pt in unique(meta$patient)) {
  d <- meta[meta$patient == pt, , drop = FALSE]
  bl <- d[d$tissue == "lesional skin" & d$day == "Day0", , drop = FALSE]
  if (nrow(bl) != 1L) next
  nls <- d[d$tissue == "non-lesional skin" & d$day == "Day0", , drop = FALSE]
  for (dy in c("Day14", "Day28")) {
    po <- d[d$tissue == "lesional skin" & d$day == dy, , drop = FALSE]
    if (nrow(po) != 1L) next
    rows[[length(rows) + 1L]] <- data.frame(
      patient = pt, arm = bl$arm, day = dy,
      pasi_baseline = bl$pasi, pasi_post = po$pasi,
      pasi_improvement_pct = if (is.finite(bl$pasi) && is.finite(po$pasi) && bl$pasi > 0)
        (bl$pasi - po$pasi) / bl$pasi * 100 else NA_real_,
      residual_score_baseline = bl$residual_score,
      residual_score_post = po$residual_score,
      residual_drop = bl$residual_score - po$residual_score,
      residual_score_NL = if (nrow(nls) == 1L) nls$residual_score else NA_real_,
      stringsAsFactors = FALSE)
  }
}
pat <- do.call(rbind, rows)
pat$residual_above_NL <- with(pat, ifelse(is.finite(residual_score_NL),
                                          residual_score_post > residual_score_NL, NA))
write.csv(pat, file.path(res_dir, "GSE136757_residual_signature_patient.csv"), row.names = FALSE)

## ---- tests (patient-level) ----
sp_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L) return(c(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  c(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}
grp_list <- list(Active = c("PF.30", "PF.100"), Placebo = "Placebo")
tests <- list()
for (nm in names(grp_list)) {
  grp <- grp_list[[nm]]; d <- pat[pat$arm %in% grp, , drop = FALSE]
  for (dy in c("Day14", "Day28")) {
    z <- d[d$day == dy, , drop = FALSE]
    for (v in c("residual_score_post", "residual_drop", "residual_score_baseline")) {
      r <- sp_test(z[[v]], z$pasi_improvement_pct)
      tests[[length(tests) + 1L]] <- data.frame(
        cohort = nm, day = dy, predictor = v, n = r[["n"]],
        spearman_rho = r[["rho"]], p = r[["p"]], stringsAsFactors = FALSE)
    }
    ## Does residual programme at the post timepoint distinguish PASI75?
    ok <- is.finite(z$residual_score_post) & is.finite(z$pasi_improvement_pct)
    if (sum(ok) >= 4L && length(unique(z$pasi_improvement_pct[ok] >= 75)) == 2L) {
      w <- suppressWarnings(wilcox.test(z$residual_score_post[ok][z$pasi_improvement_pct[ok] >= 75],
                                        z$residual_score_post[ok][z$pasi_improvement_pct[ok] < 75]))
      tests[[length(tests) + 1L]] <- data.frame(
        cohort = nm, day = dy, predictor = "residual_score_post_PASI75_vs_below",
        n = sum(ok), spearman_rho = NA_real_, p = w$p.value, stringsAsFactors = FALSE)
    }
  }
}
tt <- do.call(rbind, tests)
tt$fdr_within_cohort <- ave(tt$p, tt$cohort, FUN = function(p) p.adjust(p, "BH"))
write.csv(tt, file.path(res_dir, "GSE136757_residual_signature_tests.csv"), row.names = FALSE)

cat("\n== patient-level pairing ==\n"); print(pat, row.names = FALSE)
cat("\n== tests ==\n"); print(tt, row.names = FALSE)

## ---- placebo-controlled module sanity: is the score drug-responsive? ----
act <- pat[pat$arm != "Placebo" & is.finite(pat$residual_drop), ]
plc <- pat[pat$arm == "Placebo" & is.finite(pat$residual_drop), ]
if (nrow(act) && nrow(plc)) {
  w <- suppressWarnings(wilcox.test(act$residual_drop, plc$residual_drop))
  cat(sprintf("\nresidual programme drop, active vs placebo: median %.3f vs %.3f, Wilcoxon p=%.4g\n",
              median(act$residual_drop), median(plc$residual_drop), w$p.value))
}
cat("\nDONE\n")

