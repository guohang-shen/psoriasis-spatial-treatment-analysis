## Supplementary Figure S6: etanercept residual kinetics in GSE11903.
## The plot is descriptive: it shows median residual fractions by week and
## does not imply a patient-level confidence interval for each point.
options(stringsAsFactors = FALSE)
res <- "C:/Users/15082/Desktop/Psoriasis/results"
out <- "C:/Users/15082/Desktop/Psoriasis/results/SuppFigure_S6_Etanercept_residual_kinetics.png"
d <- read.csv(file.path(res, "GSE11903_residual_by_module.csv"), stringsAsFactors = FALSE)
mods <- c("KC_stress", "IFN_I", "IFN_G", "myeloid_infl", "Th2_negctrl")
labs <- c(KC_stress="KC stress", IFN_I="IFN-I", IFN_G="IFN-gamma",
          myeloid_infl="Myeloid inflammation", Th2_negctrl="Th2 control")
cols <- c(KC_stress="#C0442E", IFN_I="#2E6DA4", IFN_G="#7A3E9D",
          myeloid_infl="#2E8B57", Th2_negctrl="#777777")
png(out, width=2400, height=1500, res=230, family="sans")
par(mar=c(5.2,5.8,4.2,1.5), mgp=c(2.7,0.8,0), bty="l")
plot(NA, xlim=c(0,12), ylim=c(-0.15,1.35), axes=FALSE,
     xlab="Etanercept treatment week", ylab="Residual fraction of baseline offset",
     main="Supplementary Figure S6 | Etanercept residual kinetics (GSE11903)",
     cex.main=1.15)
abline(h=c(0,1), lty=c(3,2), col=c("#999999","#BBBBBB"))
for (md in mods) {
  s <- d[d$module == md, ]
  s <- s[order(s$week), ]
  lines(s$week, s$median_residual_fraction, col=cols[md], lwd=2.5,
        lty=ifelse(md=="Th2_negctrl",2,1))
  points(s$week, s$median_residual_fraction, col=cols[md], pch=19, cex=1.25)
}
axis(1, at=c(0,1,2,4,12), labels=c(0,1,2,4,12))
axis(2, las=1)
legend("topright", legend=unname(labs[mods]), col=cols[mods],
       lwd=2.5, lty=c(1,1,1,1,2), pch=19, bty="n", cex=0.9)
mtext("fraction = post-treatment lesional-to-non-lesional offset / baseline offset; week-4 paired ordering n=14",
      side=1, line=3.8, cex=0.78, col="#555555")
dev.off()
cat("Supplementary Figure S6 written\n")

