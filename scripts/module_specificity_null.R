## ============================================================
## Module specificity null: are the resolution effects specific
## to the 4 disease programs, or would any matched gene set do?
## For each cohort: 2000 random gene sets matched on expression
## and detection; empirical one-sided p for the observed module.
## ============================================================
options(stringsAsFactors = FALSE)
.libPaths(c("C:/codex_r_lib_psoriasis", .libPaths()))
res <- "C:/Users/15082/Desktop/Psoriasis/results"
NDRAW <- 2000L
mods <- list(
  KC_stress    = c("KRT16","KRT17","KRT6A","KRT6B","S100A7","S100A8","S100A9","DEFB4A"),
  IFN_I        = c("ISG15","IFIT1","IFIT3","IFITM1","MX1","OAS1","OASL","STAT1"),
  IFN_G        = c("CXCL9","CXCL10","CXCL11","GBP1","GBP5","IDO1","STAT1","HLA-DRA"),
  myeloid_infl = c("S100A8","S100A9","IL1B","CXCL8","TNF","FCGR3A","NLRP3","CTSD"),
  Th2_negctrl  = c("CCL17","CCL22","CCL26","IL13","IL4R","GATA3","POSTN","TSLP"))

null_test <- function(lc, base_idx, post_idx, cohort, tag, ens2sym) {
  sym <- if (!is.null(ens2sym)) ens2sym[rownames(lc)] else rownames(lc)
  ok <- !is.na(sym) & nzchar(sym) & !duplicated(sym)
  lc <- lc[ok, , drop = FALSE]; sym <- sym[ok]
  mu  <- setNames(as.numeric(rowMeans(lc)), sym)
  det <- setNames(as.numeric(rowMeans(lc > 0)), sym)
  qb <- as.integer(cut(rank(mu, ties.method="first"), breaks=10, labels=FALSE))
  db <- as.integer(cut(rank(det, ties.method="first"), breaks=10, labels=FALSE))
  bins <- split(seq_along(sym), paste(qb, db, sep="."))
  base <- lc[, base_idx, drop=FALSE]; post <- lc[, post_idx, drop=FALSE]
  raw_stat <- function(idx) {
    d <- rowMeans(post[idx, , drop=FALSE]) - rowMeans(base[idx, , drop=FALSE])
    median(d)
  }
  set.seed(2026)
  out <- list()
  for (nm in names(mods)) {
    gu <- intersect(mods[[nm]], sym)
    if (length(gu) < 3L) next
    ri <- match(gu, sym)
    obs <- raw_stat(ri)
    ## null: same number of genes, drawn from the same expression/detection bins
    nd <- numeric(NDRAW)
    for (b in seq_len(NDRAW)) {
      pick <- unlist(lapply(ri, function(i) {
        pool <- setdiff(bins[[paste(qb[i], db[i], sep=".")]], ri)
        if (!length(pool)) return(sample(setdiff(seq_along(sym), ri), 1))
        sample(pool, 1)}))
      nd[b] <- raw_stat(pick)
    }
    out[[nm]] <- data.frame(cohort=cohort, tag=tag, module=nm, n_genes=length(gu), n_patients=ncol(base),
      observed=obs, null_median=median(nd), null_sd=sd(nd),
      z_vs_null=(obs-median(nd))/sd(nd),
      p_emp_left=mean(nd <= obs), p_emp_abs=mean(abs(nd) >= abs(obs)))
  }
  do.call(rbind, out)
}

results <- list()

## ---- GSE183047 (IL-17A) and GSE278330 (risankizumab), keratinocyte ----
for (f in c("GSE183047_pseudobulk_by_celltype.rds","GSE278330_pseudobulk_by_celltype.rds")) {
  x <- readRDS(file.path(res, f))
  syms <- x$symbols; names(syms) <- rownames(x$counts[[1]])
  m <- x$counts[["keratinocyte"]]
  cn <- colnames(m)
  is_post <- grepl("post", cn)
  lc <- log2(t(t(m)/(colSums(m)/1e6)) + 1)
  lab <- f
  r <- null_test(lc, which(!is_post), which(is_post), lab, "keratinocyte", syms)
  results[[length(results)+1L]] <- r
}

## ---- GSE228421 (IL-23) keratinocyte state pseudobulk ----
x <- readRDS(file.path(res,"GSE228421_sample_pseudobulk.rds"))
y <- readRDS(file.path(res,"GSE228421_state_pseudobulk.rds"))
ga <- read.csv(file.path(res,"GSE228421_gene_annotation.csv"), stringsAsFactors=FALSE)
e2s <- setNames(ga$symbol, ga$ensembl)
kc <- grep("\\|keratinocyte$", colnames(y$counts), value=TRUE)
base_s <- sub("\\|keratinocyte$","",kc)
mi <- match(base_s, as.character(x$meta$sample))
vis <- as.character(x$meta$visit[mi]); tis <- as.character(x$meta$tissue[mi])
sel <- tis == "L" & vis %in% c("V1","V3")
m <- y$counts[, kc[sel], drop=FALSE]
lc <- log2(t(t(m)/(colSums(m)/1e6)) + 1)
r <- null_test(lc, which(vis[sel]=="V1"), which(vis[sel]=="V3"), "GSE228421_IL23_d14", "keratinocyte", e2s)
results[[length(results)+1L]] <- r

st <- do.call(rbind, results)
st$fdr_emp_left <- p.adjust(st$p_emp_left, "BH")
st$fdr_emp_abs  <- p.adjust(st$p_emp_abs, "BH")
write.csv(st, file.path(res,"module_specificity_null.csv"), row.names=FALSE)
cat("\n=== module specificity vs 2000 matched random gene sets ===\n")
print(st[order(st$cohort, st$p_emp_left), ], row.names=FALSE, digits=3)

