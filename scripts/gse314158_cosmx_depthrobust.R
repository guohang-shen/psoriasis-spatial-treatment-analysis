## ============================================================
## GSE314158 CosMx: depth-robust module scoring, NON-CIRCULAR typing
##   * keratinocyte markers exclude KRT6A/6B/16/17 and S100A7/8/9
##   * z-scores relative to the paired NL keratinocyte reference
##   * a LARGE matched-control set (20 controls per module gene)
##   * background model fitted in NL keratinocytes:
##       z_module ~ a + b*control + c*log10(nCount)
##     corrected score = z_module - fitted background
## ============================================================
options(stringsAsFactors = FALSE)
.libPaths(c("C:/codex_r_lib_psoriasis", .libPaths()))
suppressPackageStartupMessages(library(data.table))
root <- "C:/Users/15082/Desktop/Psoriasis"; res <- file.path(root, "results")
dd <- file.path(root, "data", "GSE314158")
expr_file <- file.path(dd, "GSE314158_Cx009_S1_TMA2_CJ_6k_CSPI_Sk_250606_exprMat_file.csv.gz")
meta_file <- file.path(dd, "GSE314158_Cx009_S1_TMA2_CJ_6k_CSPI_Sk_250606_metadata_file.csv.gz")

type_modules <- list(
  keratinocyte = c("KRT14","KRT5","KRT1","KRT10","KRT15","KRTDAP","LGALS7"),
  T_NK         = c("CD3D","CD3E","TRBC2","IL7R","NKG7","GNLY","CD8A"),
  myeloid_DC   = c("LST1","TYROBP","FCER1A","CD1C","CTSS","LYZ","FCGR3A"),
  fibroblast   = c("COL1A1","COL1A2","COL3A1","DCN","LUM","COL6A1"),
  endothelial  = c("PECAM1","VWF","EMCN","KDR","ENG"),
  mast         = c("KIT","TPSAB1","TPSB2","MS4A2"),
  B_plasma     = c("CD79A","MS4A1","MZB1","JCHAIN"))
state_modules <- list(
  KC_IL17_IL36 = c("KRT16","KRT17","KRT6A","KRT6B","S100A7","S100A8","S100A9","DEFB4A","IL36G","CCL20","SERPINB3","SERPINB4"),
  KC_stress    = c("KRT16","KRT17","KRT6A","KRT6B","S100A7","S100A8","S100A9","DEFB4A"),
  KC_diff      = c("FLG","LOR","IVL","KRT1","KRT10","KRT2"),
  IFN_response = c("ISG15","IFIT1","IFIT3","IFITM1","MX1","OAS1","OASL","STAT1"),
  IFN_gamma    = c("CXCL9","CXCL10","CXCL11","GBP1","GBP5","IDO1","STAT1","HLA-DRA"),
  myeloid_infl = c("S100A8","S100A9","IL1B","CXCL8","TNF","FCGR3A","NLRP3","CTSD"),
  Th2          = c("CCL17","CCL22","CCL26","IL13","IL4R","GATA3","POSTN","TSLP"),
  Th17         = c("IL17A","IL17F","IL22","CCL20","IL36G","DEFB4A","IL23A"),
  TNF_NFkB     = c("TNF","NFKB1","NFKBIA","TNFAIP3","IL6","CXCL8","SOCS3","CCL2"),
  keratin_cyc  = c("MKI67","TOP2A","CCNB1","CDC20","UBE2C","CENPF"))
MULT <- 20L

header <- names(fread(expr_file, nrows = 0L, showProgress = FALSE))
genes_all <- setdiff(header, c("fov", "cell_ID"))
keys <- fread(expr_file, select = c("fov","cell_ID"), showProgress = FALSE)
meta <- fread(meta_file, select = c("fov","cell_ID","nCount_RNA","nFeature_RNA"), showProgress = FALSE)
setkey(meta, fov, cell_ID)
d0 <- meta[keys]
stopifnot(nrow(d0) == nrow(keys))
sm <- fread(file.path(res, "GSE314158_spatial_fov_sample_map.csv"))
d0$group   <- setNames(sm$group,   as.character(sm$fov))[as.character(d0$fov)]
d0$patient <- setNames(sm$patient, as.character(sm$fov))[as.character(d0$fov)]
ncount <- as.numeric(d0$nCount_RNA); ncount[!is.finite(ncount) | ncount <= 0] <- 1
ld <- log10(ncount)
ok <- is.finite(d0$nCount_RNA) & d0$nCount_RNA > 0 & !is.na(d0$group) & d0$group %in% c("LS","NL")
cat("matrix cells:", nrow(d0), "| LS/NL usable:", sum(ok),
    "| patients:", length(unique(d0$patient[ok])), "\n")

## ---- cell typing from stress-free markers ----
tm <- intersect(unique(unlist(type_modules)), genes_all)
X <- fread(expr_file, select = c("fov","cell_ID", tm), showProgress = FALSE)
Xn <- as.matrix(X[, ..tm]); rm(X); gc()
Xn <- log1p(sweep(Xn, 1L, ncount / 1e4, "/"))
ts <- sapply(type_modules, function(g) { ids <- intersect(g, tm)
  if (!length(ids)) return(rep(NA_real_, nrow(Xn))); rowMeans(Xn[, ids, drop = FALSE]) })
top <- apply(ts, 1L, max, na.rm = TRUE)
second <- apply(ts, 1L, function(x) sort(x, decreasing = TRUE)[2L])
ctype <- colnames(ts)[max.col(ts, ties.method = "first")]
ctype[!is.finite(top) | top < 0.10 | (top - second) < 0.02] <- "ambiguous"
d0$cell_type <- ctype
rm(X, Xn, ts); gc()
cat("\ncell types (all cells):\n"); print(table(d0$cell_type))

ref <- which(ok & d0$group == "NL" & d0$cell_type == "keratinocyte")
cat("reference (NL keratinocyte) cells:", length(ref), "\n")

chunks <- split(genes_all, ceiling(seq_along(genes_all) / 1500))
mu <- setNames(rep(NA_real_, length(genes_all)), genes_all); sdv <- mu; dtr <- mu
for (ch in seq_along(chunks)) {
  gg <- chunks[[ch]]
  Z <- fread(expr_file, select = c("fov","cell_ID", gg), showProgress = FALSE)
  Zn <- as.matrix(Z[, ..gg]); rm(Z); gc()
  Zn <- log1p(sweep(Zn, 1L, ncount / 1e4, "/"))
  Zr <- Zn[ref, , drop = FALSE]
  mu[gg] <- colMeans(Zr); sdv[gg] <- apply(Zr, 2L, sd); dtr[gg] <- colMeans(Zr > 0)
  rm(Zn, Zr); gc(); cat("  chunk", ch, "of", length(chunks), "done\n")
}
usable <- is.finite(mu) & is.finite(sdv) & sdv >= 0.05 & dtr >= 0.01
cat("usable genes:", sum(usable), "of", length(genes_all), "\n")
write.csv(data.frame(gene = genes_all, mean_ref = mu, sd_ref = sdv, det_ref = dtr, usable = usable),
          file.path(res, "GSE314158_cosmx_gene_stats.csv"), row.names = FALSE)

## ---- module genes + LARGE matched control sets ----
set.seed(2026)
qb <- cut(mu[usable], breaks = quantile(mu[usable], seq(0,1,length.out = 6)), labels = FALSE, include.lowest = TRUE)
db <- cut(dtr[usable], breaks = quantile(dtr[usable], seq(0,1,length.out = 6)), labels = FALSE, include.lowest = TRUE)
bin_id <- setNames(paste(qb, db, sep = "."), genes_all[usable])
univ <- split(names(bin_id), bin_id)

mod_use <- list(); ctrl_list <- list(); cov_rows <- list()
for (nm in names(state_modules)) {
  g <- intersect(state_modules[[nm]], genes_all)
  gu <- g[usable[match(g, genes_all)]]
  mod_use[[nm]] <- gu
  ctrl <- unique(unlist(lapply(gu, function(gn) {
    pool <- setdiff(univ[[bin_id[gn]]], g)
    if (!length(pool)) return(character(0)); sample(pool, min(MULT, length(pool)))
  })))
  ctrl_list[[nm]] <- ctrl
  cov_rows[[length(cov_rows) + 1L]] <- data.frame(module = nm, requested = length(state_modules[[nm]]),
    in_panel = length(g), usable = length(gu), ctrl = length(ctrl))
}
cov <- do.call(rbind, cov_rows); print(cov)
write.csv(cov, file.path(res, "GSE314158_cosmx_module_coverage.csv"), row.names = FALSE)

need <- unique(c(unlist(mod_use), unlist(ctrl_list)))
Z2 <- fread(expr_file, select = c("fov","cell_ID", need), showProgress = FALSE)
Zn <- as.matrix(Z2[, ..need]); rm(Z2); gc()
Zn <- log1p(sweep(Zn, 1L, ncount / 1e4, "/"))
Zs <- sweep(sweep(Zn, 2L, mu[need], "-"), 2L, sdv[need], "/")
score <- matrix(NA_real_, nrow(Zs), length(mod_use), dimnames = list(NULL, names(mod_use)))
naive <- score; ctrl <- score
for (nm in names(mod_use)) {
  ids <- mod_use[[nm]]; ci <- ctrl_list[[nm]]
  score[, nm] <- rowMeans(Zs[, ids, drop = FALSE])
  naive[, nm] <- rowMeans(Zn[, ids, drop = FALSE])
  if (length(ci) >= 2L) ctrl[, nm] <- rowMeans(Zs[, ci, drop = FALSE])
}
rm(Zn, Zs); gc()
sub1 <- score - ctrl

## background model fitted in non-lesional keratinocytes
corr <- score; bg <- list()
for (nm in colnames(score)) {
  fit <- lm(score[ref, nm] ~ ctrl[ref, nm] + ld[ref])
  b <- coefficients(fit); bg[[nm]] <- b
  corr[, nm] <- score[, nm] - (b[1] + b[2] * ctrl[, nm] + b[3] * ld)
}
cat("\nbackground model coefficients (NL keratinocytes):\n")
print(do.call(rbind, bg), digits = 3)

kz <- which(d0$cell_type == "keratinocyte" & ok)
cat("\n=== depth confound, keratinocytes (n =", length(kz), ") ===\n")
diag_tab <- do.call(rbind, lapply(colnames(score), function(nm) data.frame(module = nm,
  rho_naive = cor(naive[kz, nm], ld[kz], method = "spearman"),
  rho_z       = cor(score[kz, nm], ld[kz], method = "spearman"),
  rho_sub1    = cor(sub1[kz, nm], ld[kz], method = "spearman"),
  rho_bgmodel = cor(corr[kz, nm], ld[kz], method = "spearman"))))
print(diag_tab, row.names = FALSE, digits = 3)
write.csv(diag_tab, file.path(res, "GSE314158_cosmx_depth_confound.csv"), row.names = FALSE)

ki <- which(ok)
out <- cbind(d0[ki, .(fov, cell_ID, patient, group, cell_type, nCount_RNA)],
  setNames(as.data.frame(round(score[ki, , drop = FALSE], 4)), paste0(names(mod_use), "_z")),
  setNames(as.data.frame(round(sub1[ki, , drop = FALSE], 4)), paste0(names(mod_use), "_subctrl")),
  setNames(as.data.frame(round(corr[ki, , drop = FALSE], 4)), paste0(names(mod_use), "_corr")),
  setNames(as.data.frame(round(naive[ki, , drop = FALSE], 4)), paste0(names(mod_use), "_naive")))
fwrite(out, file.path(res, "GSE314158_cosmx_depthrobust_cells.csv"))
cat("cells written:", nrow(out), "x", ncol(out), "\n")

scores <- c(paste0(names(mod_use), "_naive"), paste0(names(mod_use), "_z"),
            paste0(names(mod_use), "_subctrl"), paste0(names(mod_use), "_corr"))
agg <- out[, c(lapply(.SD, median, na.rm = TRUE), list(n = .N)),
           by = .(patient, group, cell_type), .SDcols = scores]
fwrite(agg, file.path(res, "GSE314158_cosmx_depthrobust_patient_summary.csv"))
p <- merge(agg[group == "LS"], agg[group == "NL"], by = c("patient","cell_type"), suffixes = c("_LS","_NL"))
rows <- list(); k <- 0L
for (ct in unique(p$cell_type)) for (s in scores) {
  z <- p[cell_type == ct]; d <- z[[paste0(s,"_LS")]] - z[[paste0(s,"_NL")]]; d <- d[is.finite(d)]
  if (length(d) < 3L) next
  k <- k + 1L
  rows[[k]] <- data.frame(cell_type = ct, score = s, n_patients = length(d), median_delta = median(d),
    pos_frac = mean(d > 0), p_value = suppressWarnings(wilcox.test(d, mu = 0, exact = FALSE)$p.value))
}
st <- do.call(rbind, rows); st$fdr <- p.adjust(st$p_value, "BH")
write.csv(st, file.path(res, "GSE314158_cosmx_depthrobust_LS_vs_NL_stats.csv"), row.names = FALSE)
cat("\n=== CosMx keratinocytes: paired LS - NL ===\n")
print(st[st$cell_type == "keratinocyte", ], row.names = FALSE, digits = 3)
cat("\ncosmx depth-robust done\n")

