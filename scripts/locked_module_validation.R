suppressPackageStartupMessages({library(stats)})

base_dir <- "."
res_dir <- file.path(base_dir, "results")

sig <- read.csv(file.path(res_dir, "GSE117239_perturbation_signature_paired_deltas.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
scores <- read.csv(file.path(res_dir, "GSE117239_patient_scores.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)

## The existing perturbation file stores baseline minus week-1 for the
## treatment comparison.  We keep that direction and name it explicitly:
## positive values mean a fall in the program after treatment.
sig <- sig[sig$comparison == "Week1_LS_minus_BL_LS", , drop = FALSE]
sig$value_baseline_minus_week1 <- sig$value

wide <- reshape(sig[, c("subject", "drug", "response", "signature",
                        "value_baseline_minus_week1")],
                idvar = c("subject", "drug", "response"),
                timevar = "signature", direction = "wide")
names(wide) <- sub("value_baseline_minus_week1\\.", "", names(wide))

ker_cols <- intersect(c("Ker_Th17", "Ker_Th1_17", "Ker_Th22"), names(wide))
fib_cols <- intersect(c("Fib_Th17", "Fib_Th1_17", "Fib_Th22"), names(wide))
wide$Ker_perturbation_mean <- rowMeans(wide[, ker_cols, drop = FALSE], na.rm = TRUE)
wide$Fib_perturbation_mean <- rowMeans(wide[, fib_cols, drop = FALSE], na.rm = TRUE)
wide$Structural_perturbation_mean <- rowMeans(wide[, c("Ker_perturbation_mean",
                                                        "Fib_perturbation_mean")],
                                               na.rm = TRUE)

scores <- scores[, c("subject", "drug", "response", "IL17_baseline", "IL17_delta",
                     "IFN_baseline", "IFN_delta"), drop = FALSE]
dat <- merge(wide, scores, by = c("subject", "drug", "response"), all = FALSE)

auc_rank <- function(x, y) {
  ok <- is.finite(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  if (length(unique(y)) < 2L) return(NA_real_)
  n1 <- sum(y == "R"); n0 <- sum(y == "NR")
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  (sum(rank(x)[y == "R"]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

metric_cols <- c("IL17_baseline", "IL17_delta", "IFN_baseline", "IFN_delta",
                 "Ker_perturbation_mean", "Fib_perturbation_mean",
                 "Structural_perturbation_mean")
out <- do.call(rbind, lapply(split(dat, dat$drug, drop = TRUE), function(d) {
  do.call(rbind, lapply(metric_cols, function(v) {
    data.frame(drug = d$drug[1], predictor = v, n = nrow(d),
               responders = sum(d$response == "R"),
               nonresponders = sum(d$response == "NR"),
               auc_R_vs_NR = auc_rank(d[[v]], d$response),
               median_R = median(d[[v]][d$response == "R"], na.rm = TRUE),
               median_NR = median(d[[v]][d$response == "NR"], na.rm = TRUE),
               wilcox_p = tryCatch(wilcox.test(d[[v]] ~ d$response,
                                               exact = FALSE)$p.value,
                                   error = function(e) NA_real_),
               stringsAsFactors = FALSE)
  }))
}))
out$FDR_within_drug <- ave(out$wilcox_p, out$drug,
                           FUN = function(x) p.adjust(x, method = "BH"))
write.csv(dat, file.path(res_dir, "GSE117239_locked_module_patient_scores.csv"), row.names = FALSE)
write.csv(out, file.path(res_dir, "GSE117239_locked_module_validation.csv"), row.names = FALSE)
print(out, row.names = FALSE)

cat("\\nInterpretation rule: this is a locked, perturbation-derived module test; no predictor was selected using the response labels.\\n")

