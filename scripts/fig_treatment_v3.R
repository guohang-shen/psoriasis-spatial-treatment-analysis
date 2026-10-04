## ============================================================
## Fig 4/5 (v3): cross-biologic resolution of the spatial axis
##   Fig_treatment_convergence.png : A heatmap + colour bar,
##     B Stouffer forest, C matched-gene-set null
##   Fig_kinetics_residual.png     : D kinetics small multiples,
##     E week-28 residual, F negative-control summary
## ============================================================
options(stringsAsFactors = FALSE)
.libPaths(c("C:/codex_r_lib_psoriasis", .libPaths()))
res <- "C:/Users/15082/Desktop/Psoriasis/results"
st  <- read.csv(file.path(res,"cross_biologic_convergence_tests.csv"))
mt  <- read.csv(file.path(res,"cross_biologic_stouffer.csv"))
nl  <- read.csv(file.path(res,"module_specificity_null.csv"))
kk  <- read.csv(file.path(res,"resolution_kinetics_scores.csv"))
rt  <- read.csv(file.path(res,"GSE278330_residual_tests.csv"))

MODS <- c("KC_stress","IFN_I","IFN_G","myeloid_infl","Th2_negctrl")
LAB  <- c(KC_stress="KC stress", IFN_I="IFN-I", IFN_G="IFN-gamma",
          myeloid_infl="Myeloid infl.", Th2_negctrl="Th2 (control)")
DISEASE_RESP <- c("KC_stress","IFN_I","IFN_G","myeloid_infl")
arms  <- c("GSE183047_IL17A","GSE278330_risankizumab","GSE228421_IL23_d14",
           "GSE117239_UST_wk1","GSE117239_ETN_wk1")
arLab <- c("IL-17A\nwk12\n(KC)","IL-23\nwk28\n(KC)","IL-23\nd14\n(KC)",
           "IL-12/23\nwk1\n(bulk)","TNF\nwk1\n(bulk)")

M <- matrix(NA_real_, length(MODS), length(arms), dimnames=list(MODS, arms))
P <- M
for (i in seq_along(MODS)) for (j in seq_along(arms)) {
  s <- st[st$module==MODS[i] & st$cohort==arms[j] & st$tag %in% c("keratinocyte","bulk"), ]
  if (nrow(s)) { M[i,j] <- median(s$median_delta); P[i,j] <- min(s$p) }
}
R <- M
for (j in seq_len(ncol(M))) {
  sc <- max(abs(M[DISEASE_RESP, j]), na.rm=TRUE)
  if (!is.finite(sc) || sc <= 0) sc <- 1
  R[, j] <- M[, j]/sc
}
pal <- colorRampPalette(c("#2E6DA4","#A9C4E0","#F7F4F1","#E8A598","#B03A2E"))(101)

png(file.path(res,"Fig_treatment_convergence.png"), width=3400, height=2150, res=250, family="sans")
layout(matrix(c(1,2,1,3), nrow=2, byrow=TRUE), widths=c(1.18,1))

## ---------------- Panel A ----------------
par(mar=c(6.5,10.5,5.0,1.5))
plot(NA, xlim=c(0.5,length(arms)+0.5), ylim=c(0.35,6.30), axes=FALSE, xlab="", ylab="",
     main="A  Resolution response across treatment arms", cex.main=1.08)
for (i in seq_len(nrow(R))) for (j in seq_len(ncol(R))) {
  v <- R[i,j]; if (is.na(v)) next
  col <- pal[round((v+1)/2*100)+1]
  rect(j-0.5, i-0.5, j+0.5, i+0.5, col=col, border="white", lwd=2)
  star <- if (is.na(P[i,j])) "" else if (P[i,j]<0.01) "**" else if (P[i,j]<0.05) "*" else ""
  text(j, i, sprintf("%+.2f%s", v, star), cex=1.05, col=if (abs(v)>0.55) "white" else "#222222")
}
axis(1, at=1:length(arms), labels=arLab, tick=FALSE, cex.axis=0.95, line=-0.6)
axis(2, at=1:length(MODS), labels=LAB[MODS], tick=FALSE, las=1, cex.axis=1.15)
abline(h=length(MODS)-0.5, col="#333333", lwd=2, lty=2)
text(0.55, length(MODS)+0.42, "negative control", adj=0, cex=0.9, col="#666666")
mtext("response relative to the strongest disease response in the same arm", side=1, line=4.4, cex=0.85, col="#555555")
## colour bar (top strip)
cb <- seq(0.5, length(arms)+0.5, length.out=102)
for (k in seq_len(101)) rect(cb[k], 5.72, cb[k+1], 5.97, col=pal[k], border=NA)
rect(0.5, 5.72, length(arms)+0.5, 5.97, border="#888888")
text(c(0.5, 3.0, 5.5), 6.10, c("-1","0","+1"), cex=0.9)
text(6.0, 5.85, "normalised\nresponse", adj=0, cex=0.85, col="#555555")
text(3.0, 6.22, "collapsed  <--  relative effect  -->  retained", cex=0.8, col="#555555")

## ---------------- Panel B ----------------
par(mar=c(5.5,10.5,5.0,4.5))
m2 <- mt[mt$tag != "all_cells", ]
m2$key <- paste0(LAB[m2$module], ifelse(m2$tag=="bulk"," [bulk]"," [KC]"))
m2 <- m2[order(match(m2$module, rev(MODS))), ]
y <- seq_len(nrow(m2))
plot(NA, xlim=c(-0.5, 5.6), ylim=c(0.3, nrow(m2)+0.9), axes=FALSE,
     xlab="Stouffer Z  (positive = consistent resolution)", ylab="",
     main="B  Weighted meta-analysis across treatment arms", cex.main=1.25)
abline(v=0, lty=3)
cols <- ifelse(m2$fdr<0.05, "#1F6FB4", "#BBBBBB")
segments(0, y, m2$Z, y, col=cols, lwd=3.4); points(m2$Z, y, pch=19, col=cols, cex=1.7)
axis(1); axis(2, at=y, labels=m2$key, las=1, tick=FALSE, cex.axis=1.0)
text(5.55, y, sprintf("q=%.2g", m2$fdr), adj=1, cex=0.8, col="#555555")
text(m2$Z + 0.28, y, sprintf("%d/%d", m2$n_same_sign, m2$n_arms), adj=0, cex=0.78, col="#777777")
mtext("arms agreeing in sign", side=3, at=2.6, line=0.2, cex=0.82, col="#666666")

## ---------------- Panel C ----------------
par(mar=c(11.5,10.5,4.0,4.5))
sel <- nl[nl$module %in% MODS, ]
sel$cohort2 <- ifelse(grepl("183047", sel$cohort), "IL-17A",
               ifelse(grepl("278330", sel$cohort), "IL-23 wk28", "IL-23 d14"))
sel$lbl <- paste0(LAB[sel$module], "\n", sel$cohort2)
sel <- sel[order(match(sel$module, MODS), sel$cohort2), ]
bp <- barplot(-sel$z_vs_null, names.arg=sel$lbl, las=2, cex.names=0.8,
        col=ifelse(sel$p_emp_left<0.01,"#1F6FB4","#C9C9C9"), border=NA,
        ylab="effect vs matched random gene sets (z)",
        main="C  Specificity against matched random gene sets", cex.main=1.25)
abline(h=0)
text(bp, -sel$z_vs_null + 0.35, sprintf("p=%.3g", sel$p_emp_left), cex=0.62, col="#444444", srt=90, adj=0)
legend("topright", legend=c("empirical p < 0.01","not significant"), fill=c("#1F6FB4","#C9C9C9"),
       border=NA, bty="n", cex=0.9)
dev.off()
cat("Fig_treatment_convergence.png written\n")

## ============================================================
## Fig 5: kinetics small multiples + week-28 residual
## ============================================================
kk$rel <- NA_real_
for (cu in unique(kk$cohort)) for (md in MODS) {
  idx <- which(kk$cohort==cu & kk$module==md)
  if (!length(idx)) next
  b0 <- kk$score[idx][kk$day[idx]==0]
  sd0 <- sd(b0); if (!is.finite(sd0) || sd0 <= 0) sd0 <- 1
  kk$rel[idx] <- (kk$score[idx]-mean(b0))/sd0
}
coh <- c("GSE228421 (IL-23)","GSE183047 (IL-17A)","GSE278330 (IL-23)")
cols2 <- c("GSE228421 (IL-23)"="#1F6FB4","GSE183047 (IL-17A)"="#C0442E","GSE278330 (IL-23)"="#2E8B57")
pch2  <- c("GSE228421 (IL-23)"=19,"GSE183047 (IL-17A)"=17,"GSE278330 (IL-23)"=15)

ag <- do.call(rbind, lapply(coh, function(cu) do.call(rbind, lapply(MODS, function(md) {
  s <- kk[kk$cohort==cu & kk$module==md, ]; if (!nrow(s)) return(NULL)
  a <- aggregate(rel ~ day, data=s, FUN=function(z) c(mean=mean(z), se=sd(z)/sqrt(length(z))))
  data.frame(cohort=cu, module=md, day=a$day, m=a$rel[,"mean"], se=a$rel[,"se"])
}))))
ymin <- min(ag$m - 1.96*ag$se, na.rm=TRUE); ymax <- max(ag$m + 1.96*ag$se, na.rm=TRUE)
pad <- 0.35; ylim <- c(ymin-pad, ymax+pad)

png(file.path(res,"Fig_kinetics_residual.png"), width=3600, height=2000, res=250, family="sans")
## Fixed figure coordinates keep the two major panels balanced.  A occupies
## the left 62% and B occupies the right 32%, with no detached blank block.
plot.new()

draw_kinetics <- function(md, fig, show_y=FALSE, show_x=FALSE) {
  par(fig=fig, mar=c(if (show_x) 3.4 else 1.0,
                    if (show_y) 4.1 else 1.0,
                    2.5, 0.5), new=TRUE)
  plot(NA, xlim=c(log10(1), log10(210)), ylim=ylim, axes=FALSE,
       xlab="", ylab="", main=LAB[md],
       cex.main=1.0, col.main="#111111")
  abline(h=0, lty=3, col="#999999")
  for (cu in coh) {
    s <- ag[ag$cohort==cu & ag$module==md, ]; if (!nrow(s)) next
    lt <- if (md=="Th2_negctrl") 2 else 1
    segments(log10(s$day+1), s$m-1.96*s$se, log10(s$day+1), s$m+1.96*s$se,
             col=cols2[cu], lwd=1.0)
    lines(log10(s$day+1), s$m, col=cols2[cu], lwd=1.9, lty=lt)
    points(log10(s$day+1), s$m, col=cols2[cu], pch=pch2[cu], cex=1.1, lty=lt)
  }
  ticks <- log10(c(1,4,15,85,197))
  if (show_x) axis(1, at=ticks, labels=c(0,3,14,84,196), cex.axis=0.70)
  else axis(1, at=ticks, labels=FALSE, tck=-0.025)
  if (show_y) {
    axis(2, cex.axis=0.70)
    mtext("change from baseline\n(baseline SD units)", side=2, line=2.6, cex=0.66)
  }
  box(col="#AAAAAA")
}

## A: three panels above, two panels below, with the legend in the sixth slot.
par(fig=c(0.04,0.64,0.93,0.995), mar=c(0,0,0,0), new=TRUE)
plot.new(); text(0.01,0.42,
  "A  Program-resolution kinetics across longitudinal treatment cohorts",
  adj=0, cex=1.15, font=2)
draw_kinetics("KC_stress",    c(0.04,0.23,0.53,0.90), show_y=TRUE,  show_x=FALSE)
draw_kinetics("IFN_I",        c(0.245,0.435,0.53,0.90), show_y=FALSE, show_x=FALSE)
draw_kinetics("IFN_G",        c(0.46,0.65,0.53,0.90), show_y=FALSE, show_x=FALSE)
draw_kinetics("myeloid_infl", c(0.04,0.23,0.08,0.45), show_y=TRUE,  show_x=TRUE)
draw_kinetics("Th2_negctrl",  c(0.245,0.435,0.08,0.45), show_y=FALSE, show_x=TRUE)
par(fig=c(0.46,0.65,0.08,0.45), mar=c(0,0,0,0), new=TRUE); plot.new()
legend("center", legend=c("IL-23 day 14 (KC)","IL-17A week 12 (KC)",
                          "IL-23 week 28 (KC)","vertical bar = 95% CI"),
       col=c(cols2,"#BBBBBB"), lwd=c(1.9,1.9,1.9,5), pch=c(19,17,15,NA),
       bty="n", cex=0.76, seg.len=1.9)

## B: week-28 residual, shown as a full-height right-hand panel.
## Short labels preserve plotting width in the narrow right column.
lab_short <- c(KC_stress="KC stress", IFN_I="IFN-I", IFN_G="IFN-gamma",
               myeloid_infl="Myeloid", Th2_negctrl="Th2 control")
par(fig=c(0.68,0.98,0.08,0.91), mar=c(5.0,6.4,3.8,1.0), new=TRUE)
r2 <- rt[rt$tag=="keratinocyte", ]; r2 <- r2[match(MODS, r2$module), ]
bp <- barplot(rbind(r2$med_disease_delta, r2$med_residual_delta), beside=TRUE, horiz=TRUE,
  names.arg=lab_short[r2$module], las=1, cex.names=0.82, col=c("#C0442E","#4A7FB5"), border=NA,
  xlim=c(min(-0.4, r2$med_residual_delta)-0.05, max(r2$med_disease_delta)*1.85),
  xlab="difference vs non-lesional skin (log2 units)",
  main="B  Week-28 residual\nrisankizumab (n=8)", cex.main=0.98)
abline(v=0)
xr <- pmax(r2$med_disease_delta, r2$med_residual_delta) + 0.07
ratio <- if ("residual_fraction" %in% names(r2)) r2$residual_fraction else r2$med_residual_delta / r2$med_disease_delta  # residual offset / baseline offset
lab <- ifelse(is.na(ratio), "",
       ifelse(abs(ratio) < 0.15, sprintf("residual ~0%% (%.0f%%)", 100*ratio),
       ifelse(sign(r2$med_disease_delta) != sign(r2$med_residual_delta),
              sprintf("%.1f\u00d7 baseline; sign-flipped", ratio),
              sprintf("residual %.0f%%", 100*ratio))))
text(xr, bp[2,], lab, adj=0, cex=0.62, col="#31558A")
## Put the B legend in a dedicated title strip so it cannot collide with
## the effect-size annotations.
par(fig=c(0.68,0.98,0.915,0.995), mar=c(0,0,0,0), new=TRUE); plot.new()
legend("center", horiz=TRUE,
       legend=c("baseline lesion - NL (disease)","week 28 lesion - NL (residual)"),
       fill=c("#C0442E","#4A7FB5"), border=NA, bty="n", cex=0.62,
       seg.len=1.0)
mtext("residual = week 28 lesion - non-lesional; percentage = residual offset / baseline offset; sign-flipped = week-28 lesion is on the other side of non-lesional skin",
      side=1, line=3.2, cex=0.68, col="#666666", outer=TRUE)
mtext("A  Program-resolution kinetics across longitudinal treatment cohorts", side=3, outer=TRUE,
      at=0.30, line=0.45, cex=1.15, font=2)
dev.off()
cat("Fig_kinetics_residual.png written\n")

