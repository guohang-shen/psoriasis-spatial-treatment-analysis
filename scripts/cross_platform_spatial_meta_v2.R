## ============================================================
## Cross-platform spatial meta-analysis (Visium GSE206391 v3 +
## CosMx GSE314158 LOSO-validated), direction-preserving Stouffer
## ============================================================
options(stringsAsFactors = FALSE)
.libPaths(c("C:/codex_r_lib_psoriasis", .libPaths()))
root <- "C:/Users/15082/Desktop/Psoriasis"; res <- file.path(root, "results")
vis <- read.csv(file.path(res, "GSE206391_cross_disease_delta_tests_v3.csv"), check.names = FALSE)
cos <- read.csv(file.path(res, "GSE314158_cosmx_loso_deltas.csv"))
map <- data.frame(
  visium = c("KC_stress_corr","IFN_type1_corr","IFN_g_corr","myeloid_corr","Th2_corr","Th17_corr","TNF_NFkB_corr","keratin_cyc_corr"),
  cosmx  = c("KC_stress","IFN_response","IFN_gamma","myeloid_infl","Th2","Th17","TNF_NFkB","keratin_cyc"),
  short  = c("KC stress","IFN type-I","IFN-gamma","Myeloid infl.","Th2 (neg. control)","Th17","TNF-NFkB","KC cycle"))
map$full <- c("Keratinocyte stress / alarmin (KRT6/16/17, S100A7/8/9, DEFB4A)",
              "Type-I IFN response (ISG15/OAS1/MX1/IFIT*)",
              "IFN-gamma (CXCL9/10/11, GBP1/5, IDO1)",
              "Myeloid inflammation (S100A8/9, IL1B, CXCL8)",
              "Th2 / AD axis (CCL17/22/26, IL13/IL4R, POSTN) - negative control",
              "Th17 (IL17A/F, IL22, CCL20)",
              "TNF-NFkB (TNF, NFKB1/IA, TNFAIP3, IL6)",
              "Keratinocyte cell cycle (MKI67/TOP2A/CCNB1)")
n_vis <- 12L; n_cos <- 7L
rows <- list()
for (i in seq_len(nrow(map))) {
  vr <- vis[vis$score == map$visium[i], ]
  cd <- cos[cos$module == map$cosmx[i], ]
  d_vis <- vr$pso_delta - vr$ad_delta
  p_vis <- vr$p_t117_vs_AD
  d_cos <- median(cd$delta_loso)
  alt <- if (d_cos > 0) "greater" else "less"
  p_cos <- wilcox.test(cd$delta_loso, mu = 0, alternative = alt, exact = FALSE)$p.value
  z_vis <- qnorm(1 - p_vis); z_cos <- qnorm(1 - p_cos)
  wv <- sqrt(n_vis); wc <- sqrt(n_cos)
  z_cmb <- (wv * z_vis + wc * z_cos) / sqrt(wv^2 + wc^2)
  rows[[i]] <- data.frame(module = map$full[i], short = map$short[i],
    visium_module = map$visium[i], cosmx_module = map$cosmx[i],
    n_vis = n_vis, delta_vis = d_vis, p_vis = p_vis, z_vis = z_vis,
    n_cos = n_cos, delta_cos = d_cos, p_cos = p_cos, z_cos = z_cos,
    z_combined = z_cmb, p_combined = 1 - pnorm(z_cmb),
    sign_agree = sign(d_vis) == sign(d_cos))
}
out <- do.call(rbind, rows)
out$fdr_combined <- p.adjust(out$p_combined, "BH")
out$direction <- ifelse(out$z_combined > 0, "up in type-1/17 lesional skin", "not up / down")
write.csv(out, file.path(res, "GSE_crossplatform_spatial_meta.csv"), row.names = FALSE)
cat("=== cross-platform spatial meta (Visium n=12 | CosMx n=7) ===\n")
print(out[, c("short","delta_vis","p_vis","delta_cos","p_cos","z_combined","p_combined","fdr_combined","sign_agree")],
      row.names = FALSE, digits = 3)

## ---------------- figure ----------------
png(file.path(res, "Fig_spatial_meta.png"), width = 2600, height = 1500, res = 230)
op <- par(mar = c(4.4, 12.5, 3.4, 7.5), family = "sans")
d <- out[order(out$z_combined), ]
d$y <- seq_len(nrow(d))
sig <- d$fdr_combined < 0.05
plot(NA, xlim = c(-4.4, 4.4), ylim = c(0.35, nrow(d) + 0.75), yaxt = "n", xaxt = "n",
     xlab = "signed meta-Z   (positive = enriched in type-1/17 lesional skin)", ylab = "", bty = "n")
rect(-4.4, d$y[sig] - 0.34, 4.4, d$y[sig] + 0.34, col = "#FBEFF0", border = NA)
abline(v = c(-1.96, 0, 1.96), col = c("#BBBBBB","#777777","#BBBBBB"), lty = c(2,1,2), lwd = c(1,1.3,1))
segments(d$z_vis, d$y, d$z_cos, d$y, col = "#E0E0E0", lwd = 3)
points(d$z_vis, d$y, pch = 21, bg = "#4C72B0", col = "white", cex = 1.5)
points(d$z_cos, d$y, pch = 22, bg = "#DD8452", col = "white", cex = 1.5)
points(d$z_combined, d$y, pch = 23, bg = ifelse(sig, "#9E2A2B", "#C9A227"), col = "white", cex = 2.1)
axis(1, at = c(-4, -2, 0, 2, 4))
axis(2, at = d$y, labels = d$short, las = 1, tick = FALSE, cex.axis = 1.05,
     font = ifelse(sig, 2, 1))
txt <- ifelse(d$fdr_combined < 0.001, sprintf("q=%.1e", d$fdr_combined), sprintf("q=%.3f", d$fdr_combined))
mtext(txt, side = 4, at = d$y, las = 1, line = 0.3, cex = 0.78,
      col = ifelse(sig, "#9E2A2B", "#444444"), font = ifelse(sig, 2, 1))
legend("bottomright", legend = c("Visium (n=12)", "CosMx (n=7)", "Stouffer combined"),
       pch = c(21, 22, 23), pt.bg = c("#4C72B0", "#DD8452", "#9E2A2B"), col = "white",
       bty = "n", cex = 1.0, inset = c(0.01, 0.02))
title("Cross-platform spatial replication of the type-1/17 lesional programme", cex.main = 1.12)
mtext("shaded rows: combined BH q < 0.05", side = 3, line = 0.1, cex = 0.78, col = "#666666")
par(op); dev.off()
cat("figure written\n")

