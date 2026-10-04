## ============================================================
## Independent test of the "interferon axis normalises faster than
## keratinocyte stress" claim, in a SECOND cohort and a SECOND drug class.
##
## Fig 5B (week-28 residual) currently rests on a single cohort:
## GSE278330 (risankizumab / IL-23, n=8, keratinocyte pseudobulk).
## GSE11903 is etanercept (anti-TNF) bulk skin with patient-matched
## non-lesional skin at baseline and lesional biopsies at weeks 0/1/2/4/12.
## Only within-cohort contrasts are used: no cross-platform and no
## cross-drug effect sizes are mixed.
## ============================================================
options(stringsAsFactors = FALSE)
.libPaths(c("C:/codex_r_lib_psoriasis", .libPaths()))
suppressPackageStartupMessages({library(GEOquery); library(Biobase)})

root <- "C:/Users/15082/Desktop/Psoriasis"
out  <- file.path(root, "results")

mods <- list(
  KC_stress    = c("KRT16","KRT17","KRT6A","KRT6B","S100A7","S100A8","S100A9","DEFB4A"),
  IFN_I        = c("ISG15","IFIT1","IFIT3","IFITM1","MX1","OAS1","OASL","STAT1"),
  IFN_G        = c("CXCL9","CXCL10","CXCL11","GBP1","GBP5","IDO1","STAT1","HLA-DRA"),
  myeloid_infl = c("S100A8","S100A9","IL1B","CXCL8","TNF","FCGR3A","NLRP3","CTSD"),
  Th2_negctrl  = c("CCL17","CCL22","CCL26","IL13","IL4R","GATA3","POSTN","TSLP"))

e  <- getGEO(filename = file.path(root, "data/GSE11903/GSE11903_series_matrix.txt.gz"), getGPL = FALSE)
pd <- pData(e); x <- exprs(e)
field <- function(col) trimws(sub("^[^:]+:\\s*", "", pd[[col]]))
meta <- data.frame(
  gsm       = pd$geo_accession,
  patient   = field("characteristics_ch1"),
  week      = as.integer(sub(".*?(\\d+).*", "\\1", field("characteristics_ch1.1"))),
  responder = field("characteristics_ch1.2") == "responder",
  tissue    = field("characteristics_ch1.3"))
cat("samples:", nrow(meta), " patients:", length(unique(meta$patient)), "\n")
print(table(meta$week, meta$tissue))

raw <- readLines(gzfile(file.path(root, "data/GPL571.annot.gz")))
h   <- which(startsWith(raw, "ID\tGene title\tGene symbol"))[1]
ann <- read.delim(text = paste(raw[h:length(raw)], collapse = "\n"),
                  check.names = FALSE, quote = "", fill = TRUE)
sp  <- strsplit(as.character(ann$`Gene symbol`), "///", fixed = TRUE)
aa  <- data.frame(ID = rep(ann$ID, lengths(sp)),
                  symbol = trimws(unlist(sp, use.names = FALSE)))
aa  <- aa[aa$ID %in% rownames(x) & nzchar(aa$symbol) & aa$symbol != "---", ]

sym_expr <- function(g) {
  ids <- unique(aa$ID[aa$symbol == g])
  if (!length(ids)) return(NULL)
  if (length(ids) == 1L) as.numeric(x[ids, ]) else
    as.numeric(apply(x[ids, , drop = FALSE], 2, median))
}
build_m <- function(genes) {
  cols <- list()
  for (g in genes) { v <- sym_expr(g); if (!is.null(v)) cols[[g]] <- v }
  if (!length(cols)) return(NULL)
  m <- do.call(cbind, cols); colnames(m) <- names(cols); rownames(m) <- colnames(x)
  m
}

nl0 <- meta$gsm[meta$week == 0 & meta$tissue == "non-lesional"]
ls0 <- meta$gsm[meta$week == 0 & meta$tissue == "lesional"]
cat("baseline NL:", length(nl0), " baseline LS:", length(ls0), "\n")

rows <- list(); gene_use <- list()
for (M in names(mods)) {
  m <- build_m(mods[[M]])
  if (is.null(m)) next
  mu  <- colMeans(m[nl0, , drop = FALSE], na.rm = TRUE)
  sdv <- apply(m[nl0, , drop = FALSE], 2, sd, na.rm = TRUE)
  ok  <- is.finite(sdv) & sdv > 0.05
  if (sum(ok) < 3L) { cat("module", M, ": too few testable genes\n"); next }
  m <- m[, ok, drop = FALSE]; mu <- mu[ok]; sdv <- sdv[ok]
  z  <- sweep(sweep(m, 2, mu, "-"), 2, sdv, "/")
  sc <- setNames(as.numeric(rowMeans(z, na.rm = TRUE)), rownames(z))
  gene_use[[M]] <- colnames(m)
  for (p in unique(meta$patient)) {
    d <- meta[meta$patient == p, ]
    nlp <- d$gsm[d$week == 0 & d$tissue == "non-lesional"]
    lsp <- d$gsm[d$week == 0 & d$tissue == "lesional"]
    if (length(nlp) != 1L || length(lsp) != 1L) next
    for (wk in c(0, 1, 2, 4, 12)) {
      post <- if (wk == 0) lsp else d$gsm[d$week == wk & d$tissue == "lesional"]
      if (length(post) != 1L) next
      dis <- sc[lsp] - sc[nlp]; res <- sc[post] - sc[nlp]
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = "GSE11903_etanercept_bulk", module = M, patient = p,
        responder = d$responder[1], week = wk, n_genes = ncol(m),
        disease_offset = dis, post_offset = res,
        residual_fraction = if (is.finite(dis) && abs(dis) > 0.15) res / dis else NA_real_)
    }
  }
}
pv <- do.call(rbind, rows)
pv$cohort <- as.character(pv$cohort)
write.csv(pv, file.path(out, "GSE11903_residual_patient_deltas.csv"), row.names = FALSE)
cat("\ngenes kept per module:\n"); print(sapply(gene_use, length))

summ <- do.call(rbind, lapply(split(pv, list(pv$module, pv$week)), function(d) {
  d <- d[is.finite(d$post_offset), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  fr <- d$residual_fraction[is.finite(d$residual_fraction)]
  wt <- tryCatch(wilcox.test(d$post_offset, mu = 0, exact = FALSE)$p.value,
                 error = function(e) NA_real_)
  data.frame(cohort = d$cohort[1], module = d$module[1], week = d$week[1], n = nrow(d),
    median_disease_offset = median(d$disease_offset, na.rm = TRUE),
    median_post_offset    = median(d$post_offset, na.rm = TRUE),
    p_post_vs_NL = wt, n_with_fraction = length(fr),
    median_residual_fraction = if (length(fr)) median(fr) else NA_real_)
}))
summ <- summ[order(summ$module, summ$week), , drop = FALSE]
summ$fdr_post_vs_NL <- p.adjust(summ$p_post_vs_NL, "BH")
write.csv(summ, file.path(out, "GSE11903_residual_by_module.csv"), row.names = FALSE)
cat("\n=== GSE11903 residual by module (bulk, etanercept) ===\n")
print(summ, row.names = FALSE, digits = 3)

wide <- reshape(pv[pv$week %in% c(4, 12), c("patient","responder","module","week","post_offset","residual_fraction")],
                idvar = c("patient","responder","week"), timevar = "module", direction = "wide")
pair_tests <- list()
for (wk in c(4, 12)) {
  w <- wide[wide$week == wk, , drop = FALSE]
  for (cmp in list(c("KC_stress","IFN_I"), c("KC_stress","IFN_G"), c("KC_stress","Th2_negctrl"))) {
    a <- w[[paste0("post_offset.", cmp[1])]]; b <- w[[paste0("post_offset.", cmp[2])]]
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 5L) next
    pw <- wilcox.test(a[ok], b[ok], paired = TRUE, exact = FALSE)$p.value
    pair_tests[[length(pair_tests) + 1L]] <- data.frame(week = wk,
      contrast = paste(cmp[1], "vs", cmp[2]), n = sum(ok),
      median_diff_KCa_minus_b = median(a[ok] - b[ok]),
      n_KC_gt_other = sum(a[ok] > b[ok]), p_paired_wilcoxon = pw)
  }
}
pt <- do.call(rbind, pair_tests); pt$fdr <- p.adjust(pt$p_paired_wilcoxon, "BH")
write.csv(pt, file.path(out, "GSE11903_residual_paired_contrasts.csv"), row.names = FALSE)
cat("\n=== paired within-patient contrasts (post-treatment offset vs non-lesional) ===\n")
print(pt, row.names = FALSE, digits = 3)

cat("\n=== week 12 per-patient residual offset (baseline-NL SD units) ===\n")
w12 <- wide[wide$week == 12, , drop = FALSE]
print(w12[, c("patient","responder","post_offset.KC_stress","post_offset.IFN_I","post_offset.IFN_G","post_offset.myeloid_infl","post_offset.Th2_negctrl")],
      row.names = FALSE, digits = 3)
