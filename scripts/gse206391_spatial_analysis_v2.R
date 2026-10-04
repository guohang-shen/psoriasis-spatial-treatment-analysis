## GSE206391 spatial layer -- final statistics + figure  (v2 inputs, robust version)
options(stringsAsFactors = FALSE)
.libPaths(c("C:/codex_r_lib_psoriasis", .libPaths()))
root <- "C:/Users/15082/Desktop/Psoriasis"; out <- file.path(root, "results")
r <- read.csv(file.path(out, "GSE206391_spot_axis_score_v2.csv"))
r <- r[!is.na(r$disease), ]
ord <- c("psoriasis","atopic dermatitis","lichen planus")
cols <- c("psoriasis" = "#C0392B", "atopic dermatitis" = "#2471A3", "lichen planus" = "#1E8449")
vars <- c("axis_score","axis_corrected","KC_stress","KC_stress_corr","KC_basal_corr","KC_supra_corr",
          "fibroblast_corr","endothelial_corr","Tcell_corr","myeloid_corr","IFN_corr")
vars <- vars[vars %in% colnames(r)]

pm <- NULL
for (v in vars) {
  a <- aggregate(r[[v]], by = list(disease = r$disease, patient = r$patient, biopsy_type = r$biopsy_type),
                 FUN = function(z) median(z, na.rm = TRUE))
  colnames(a)[4] <- v
  pm <- if (is.null(pm)) a else merge(pm, a, by = c("disease","patient","biopsy_type"), all = TRUE)
}
write.csv(pm, file.path(out, "GSE206391_patient_medians_v2.csv"), row.names = FALSE)
cat("patient-level rows:", nrow(pm), "\n")

paired_delta <- function(v) {
  L <- pm[pm$biopsy_type == "LESIONAL", c("disease","patient",v)]
  N <- pm[pm$biopsy_type == "NON LESIONAL", c("disease","patient",v)]
  m <- merge(L, N, by = c("disease","patient")); colnames(m)[3:4] <- c("L","N")
  m$delta <- m$L - m$N; m
}
perm_p <- function(delta, grp, n = 20000) {
  obs <- median(delta[grp == "type1/17"]) - median(delta[grp == "AD"])
  set.seed(7)
  s <- replicate(n, { g <- sample(grp); median(delta[g == "type1/17"]) - median(delta[g == "AD"]) })
  (sum(s >= obs) + 1) / (n + 1)
}
tab <- do.call(rbind, lapply(vars, function(v) {
  m <- paired_delta(v); m <- m[!is.na(m$delta), ]
  g <- ifelse(m$disease == "atopic dermatitis", "AD", "type1/17")
  data.frame(score = v,
    pso_delta = round(median(m$delta[m$disease == "psoriasis"]), 4),
    lp_delta  = round(median(m$delta[m$disease == "lichen planus"]), 4),
    ad_delta  = round(median(m$delta[m$disease == "atopic dermatitis"]), 4),
    n_type117 = sum(g == "type1/17"), n_AD = sum(g == "AD"),
    p_wilcox = signif(suppressWarnings(wilcox.test(m$delta[g == "type1/17"], m$delta[g == "AD"], alternative = "greater")$p.value), 3),
    p_perm = round(perm_p(m$delta, g), 4))
}))
cat("\n=== paired lesional minus non-lesional, corrected scores ===\n"); print(tab, row.names = FALSE)
write.csv(tab, file.path(out, "GSE206391_cross_disease_delta_tests.csv"), row.names = FALSE)
md <- paired_delta("axis_corrected"); saveRDS(md, file.path(out, "_gse206391_delta_corrected.rds"))
write.csv(md, file.path(out, "GSE206391_patient_delta_corrected.csv"), row.names = FALSE)

cat("\ndepth confound: rho(axis, log10 det) =", round(cor(r$axis_score, log10(r$genes_detected), method="spearman"), 3),
    "| R2(axis ~ matched control) =", round(cor(r$axis_score, r$axis_ctrl)^2, 3),
    "| rho(corrected, log10 det) =", round(cor(r$axis_corrected, log10(r$genes_detected), method="spearman"), 3), "\n")

## ---------- figure ----------
png(file.path(out, "GSE206391_spatial_layer.png"), width = 1900, height = 1350, res = 130)
par(mfrow = c(2,2), mar = c(4.6,4.6,3,1))
s <- r[sample(nrow(r), 20000), ]
plot(log10(s$genes_detected), s$axis_score, pch = 16, cex = 0.35, col = paste0(cols[s$disease],"55"),
     xlab = "log10 genes detected per spot", ylab = "naive axis score",
     main = sprintf("A  Naive score is depth-driven (rho=%.2f)", cor(r$axis_score, log10(r$genes_detected), method="spearman")))
abline(lm(axis_score ~ log10(genes_detected), data = r), lwd = 2)
legend("topright", ord, pch = 16, col = cols[ord], bty = "n", cex = 0.85)
plot(s$axis_ctrl, s$axis_score, pch = 16, cex = 0.35, col = "#33333322",
     xlab = "expression-matched control score", ylab = "naive axis score",
     main = sprintf("B  Matched control explains R2=%.2f", cor(r$axis_score, r$axis_ctrl)^2))
abline(0, 1, lty = 2, col = "red", lwd = 2)
g <- ifelse(md$disease == "atopic dermatitis", "AD", "type-1/17\n(psoriasis+LP)")
pv <- suppressWarnings(wilcox.test(md$delta[g == "type-1/17\n(psoriasis+LP)"], md$delta[g == "AD"], alternative = "greater")$p.value)
boxplot(delta ~ g, data = data.frame(delta = md$delta, g = g), col = c("#8E44AD66","#2471A366"),
        outline = FALSE, ylab = "paired delta (corrected axis)",
        main = sprintf("C  Corrected axis, lesional - non-lesional\nWilcoxon p=%.3f,  permutation p=%.4f", pv, tab$p_perm[tab$score=="axis_corrected"]))
stripchart(delta ~ g, data = data.frame(delta = md$delta, g = g), vertical = TRUE, method = "jitter", pch = 21, bg = "white", add = TRUE)
abline(h = 0, lty = 3)
par(mfrow = c(2,2), mar = c(4.6,4.6,3,1))
for (d in ord) {
  z <- r[r$disease == d & r$biopsy_type == "LESIONAL", ]
  lib <- names(sort(tapply(z$axis_corrected, z$library_id, median), decreasing = TRUE))[1]
  z <- z[z$library_id == lib, ]
  lim <- quantile(r$axis_corrected, c(0.02, 0.98), na.rm = TRUE)
  cc <- colorRampPalette(c("#2166AC","grey95","#B2182B"))(64)
  ii <- pmin(64, pmax(1, round((z$axis_corrected - lim[1]) / diff(lim) * 63) + 1))
  plot(z$array_col, max(z$array_row) - z$array_row, pch = 15, cex = 0.7, asp = 1, col = cc[ii],
       xlab = "", ylab = "", main = sprintf("%s lesional (%s)", d, lib))
}
dev.off()
cat("figure written\n")


