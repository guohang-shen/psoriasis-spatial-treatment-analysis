## ============================================================
## Cross-biologic convergence of the spatial resolution axis.
## Assemble per-patient module deltas from 5 treatment arms,
## build the convergence table, Stouffer meta across treatment arms,
## and the consensus figure.
## ============================================================
options(stringsAsFactors = FALSE)
.libPaths(c("C:/codex_r_lib_psoriasis", .libPaths()))
res <- "C:/Users/15082/Desktop/Psoriasis/results"
MODS <- c("KC_stress","IFN_I","IFN_G","myeloid_infl","Th2_negctrl")
lp <- list()

## ---- Arm 1: GSE228421, IL-23 blockade, keratinocyte pseudobulk ----
a <- read.csv(file.path(res,"GSE228421_spatial_module_scores.csv"), stringsAsFactors=FALSE)
a <- a[a$tag == "keratinocyte", ]
w <- reshape(a[, c("module","patient","vt","score")], idvar=c("module","patient"),
             timevar="vt", direction="wide")
lp[[1]] <- data.frame(cohort="GSE228421_IL23_d14", drug="IL-23", tag="keratinocyte",
  patient=w$patient, module=w$module, delta=w[["score.V3_L"]] - w[["score.V1_L"]])

## ---- Arms 2/3: GSE183047 (IL-17A) and GSE278330 (risankizumab) ----
b <- read.csv(file.path(res,"GSE183047_GSE278330_module_scores.csv"), stringsAsFactors=FALSE)
b$week <- ifelse(grepl("week", b$col), sub(".*week([0-9]+).*","\\1", b$col), "12")
b$tp2 <- ifelse(b$tp == "baseline", "baseline", paste0("wk", b$week))
b <- b[b$tag %in% c("keratinocyte","all_cells"), ]
mk <- function(coh, tps, tagname, drug) {
  s <- b[b$cohort == coh & b$tag == tagname & b$tp2 %in% tps, ]
  ww <- reshape(s[, c("module","patient","tp2","score")], idvar=c("module","patient"),
                timevar="tp2", direction="wide")
  data.frame(cohort=coh, drug=drug, tag=tagname, patient=ww$patient, module=ww$module,
             delta=ww[[paste0("score.", tps[2])]] - ww[["score.baseline"]])
}
lp[[2]] <- mk("GSE183047_IL17A", c("baseline","wk12"), "keratinocyte", "IL-17A")
lp[[3]] <- mk("GSE278330_risankizumab", c("baseline","wk28"), "keratinocyte", "IL-23")
lp[[4]] <- mk("GSE183047_IL17A", c("baseline","wk12"), "all_cells", "IL-17A")
lp[[5]] <- mk("GSE278330_risankizumab", c("baseline","wk28"), "all_cells", "IL-23")

## ---- Arms 6/7: GSE117239 bulk week1, Ustekinumab and Etanercept ----
c1 <- read.csv(file.path(res,"GSE117239_week1_module_delta_scores.csv"), stringsAsFactors=FALSE)
c1 <- c1[c1$drug == "Ustekinumab 90 mg", c("subject","drug","module","delta_score")]
c1 <- data.frame(cohort="GSE117239_UST_wk1", drug="IL-12/23", tag="bulk", patient=c1$subject,
                 module=c1$module, delta=c1$delta_score)
c2 <- read.csv(file.path(res,"GSE117239_week1_module_delta_scores.csv"), stringsAsFactors=FALSE)
c2 <- c2[c2$drug == "Etanercept", c("subject","drug","module","delta_score")]
c2 <- data.frame(cohort="GSE117239_ETN_wk1", drug="TNF", tag="bulk", patient=c2$subject,
                 module=c2$module, delta=c2$delta_score)
lp[[6]] <- c1; lp[[7]] <- c2
pv <- do.call(rbind, lp)
pv <- pv[is.finite(pv$delta), ]
write.csv(pv, file.path(res,"cross_biologic_patient_deltas.csv"), row.names=FALSE)

## ---- convergence table ----
tests <- list()
key <- paste(pv$cohort, pv$tag, sep = "||")
for (k in unique(key)) for (nm in MODS) {
  s <- pv[key == k & pv$module == nm, ]
  d <- s$delta[is.finite(s$delta)]
  if (length(d) < 4L) next
  tests[[length(tests)+1L]] <- data.frame(cohort=s$cohort[1], drug=s$drug[1], tag=s$tag[1],
    module=nm, n=length(d), median_delta=median(d), frac_negative=mean(d<0),
    p=suppressWarnings(wilcox.test(d, mu=0, exact=FALSE)$p.value))
}
st <- do.call(rbind, tests)
st$fdr <- p.adjust(st$p,"BH")
st <- st[order(st$module, st$cohort), ]
write.csv(st, file.path(res,"cross_biologic_convergence_tests.csv"), row.names=FALSE)

## ---- Stouffer meta across treatment arms (sign-consistent) ----
indep <- c("GSE228421_IL23_d14","GSE183047_IL17A","GSE278330_risankizumab","GSE117239_UST_wk1")
meta <- list()
for (nm in MODS) for (tg in c("keratinocyte","all_cells","bulk")) {
  s <- st[st$module == nm & st$tag == tg, ]
  if (tg != "bulk") s <- s[s$cohort %in% c("GSE228421_IL23_d14","GSE183047_IL17A","GSE278330_risankizumab"), ]
  else s <- s[s$cohort %in% c("GSE117239_UST_wk1","GSE117239_ETN_wk1"), ]
  if (nrow(s) < 2L) next
  ## signed z from one-sided Wilcoxon on the expected (negative) direction
  z <- qnorm(pmax(s$p/2, 1e-300), lower.tail = s$median_delta < 0)
  z <- ifelse(s$median_delta < 0, abs(z), -abs(z))
  wts <- sqrt(s$n); Z <- sum(wts*z)/sqrt(sum(wts^2))
  meta[[length(meta)+1L]] <- data.frame(tag=tg, module=nm, n_arms=nrow(s),
    n_patients=sum(s$n), n_same_sign=sum(sign(s$median_delta)==sign(median(s$median_delta))),
    median_of_medians=median(s$median_delta), Z=Z, p_one_sided=pnorm(Z, lower.tail=FALSE))
}
mt <- do.call(rbind, meta)
mt$fdr <- p.adjust(mt$p_one_sided,"BH")
write.csv(mt, file.path(res,"cross_biologic_stouffer.csv"), row.names=FALSE)
cat("\n=== convergence per arm ===\n"); print(st, row.names=FALSE, digits=3)
cat("\n=== Stouffer meta ===\n"); print(mt, row.names=FALSE, digits=3)

