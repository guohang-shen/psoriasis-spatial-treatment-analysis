## ============================================================
## Fig 6: placebo-controlled contrast, GSE136757
##   PF-06700841 (TYK2/JAK1) vs placebo, patient-level paired
##   module drops.  Source tables already computed by
##   pilot_gse136757_tyk2_placebo.R  -- this script only plots.
## ============================================================
options(stringsAsFactors = FALSE)
.libPaths(c("C:/codex_r_lib_psoriasis", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2); library(gridExtra); library(grid)
})
res <- "C:/Users/15082/Desktop/Psoriasis/results"

st <- read.csv(file.path(res, "GSE136757_placebo_module_stats.csv"), stringsAsFactors = FALSE)
ps <- read.csv(file.path(res, "GSE136757_placebo_paired_module_scores.csv"), stringsAsFactors = FALSE)

LAB <- c(keratinocyte_activation        = "Keratinocyte activation",
         KC_IL17_IL36                   = "KC IL-17 / IL-36",
         KC_IL17_response_no_keratins   = "KC IL-17 response",
         IFN_response                   = "Interferon response",
         myeloid_inflammatory           = "Myeloid inflammatory",
         keratinocyte_identity_clean    = "KC identity (neg. ctrl)",
         KC_differentiation             = "KC differentiation (neg. ctrl)")
ORD <- c("keratinocyte_activation","KC_IL17_IL36","KC_IL17_response_no_keratins",
         "IFN_response","myeloid_inflammatory",
         "keratinocyte_identity_clean","KC_differentiation")

st$module   <- factor(st$module, levels = ORD)
st$label    <- factor(LAB[as.character(st$module)], levels = LAB[ORD])
st$role     <- ifelse(as.character(st$module) %in% names(LAB)[6:7], "negative control", "disease program")
st$fdr_txt  <- ifelse(is.na(st$fdr), "na",
                      ifelse(st$fdr < 1e-3, sprintf("FDR %.1e", st$fdr),
                             sprintf("FDR %.3f", st$fdr)))

dumbbell <- function(dd, ttl) {
  ggplot(dd, aes(y = label)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    geom_segment(aes(x = median_drop_placebo, xend = median_drop_active,
                     yend = label), colour = "grey40", linewidth = 0.7) +
    geom_point(aes(x = median_drop_placebo), colour = "grey35", size = 2.6) +
    geom_point(aes(x = median_drop_active),  colour = "#2E6DA4", size = 2.6) +
    geom_text(aes(x = pmax(median_drop_placebo, median_drop_active),
                  label = fdr_txt), hjust = -0.12, size = 2.6, colour = "grey25") +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.42))) +
    facet_grid(role ~ ., scales = "free_y", space = "free_y") +
    labs(title = ttl, x = "patient-level module drop (baseline lesional - post)",
         y = NULL,
         subtitle = sprintf("active n=%d, placebo n=%d", dd$n_active[1], dd$n_placebo[1])) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey93", colour = NA),
          plot.title = element_text(face = "bold", size = 9.5),
          plot.subtitle = element_text(size = 8, colour = "grey30"))
}
p14 <- dumbbell(st[st$day == "Day14", ], "A  Placebo-controlled module resolution, Day 14")
p28 <- dumbbell(st[st$day == "Day28", ], "B  Same contrast, Day 28")

## ---- patient-level strip for the two strongest contrasts -----------------
ps14 <- ps[ps$day == "Day14", ]
ps14$grp <- ifelse(ps14$arm == "Placebo", "Placebo", "PF-06700841")
long <- rbind(
  data.frame(grp = ps14$grp, module = "Keratinocyte activation", drop = ps14$keratinocyte_activation_drop),
  data.frame(grp = ps14$grp, module = "Interferon response",      drop = ps14$IFN_response_drop))
long <- long[is.finite(long$drop), ]
long$grp <- factor(long$grp, levels = c("Placebo", "PF-06700841"))
pC <- ggplot(long, aes(x = grp, y = drop, colour = grp)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  ggbeeswarm::geom_quasirandom(width = 0.18, size = 1.9, alpha = 0.9) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "grey20", linewidth = 0.4) +
  facet_wrap(~ module, nrow = 1) +
  scale_colour_manual(values = c(Placebo = "grey40", "PF-06700841" = "#2E6DA4"), guide = "none") +
  labs(title = "C  Individual patients, Day 14",
       x = NULL, y = "module drop (baseline - post)") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 9.5))

g <- arrangeGrob(p14, p28, pC, layout_matrix = rbind(c(1, 2), c(3, 3)), heights = c(1.15, 0.85))
ggsave(file.path(res, "Fig6_placebo_controlled.png"), g,
       width = 11, height = 8.4, dpi = 320, bg = "white")

cat("\n=== GSE136757 placebo-controlled module drops ===\n")
print(st[order(st$day, st$module), c("day","module","n_active","n_placebo",
                                     "median_drop_active","median_drop_placebo",
                                     "wilcox_p_two_sided","fdr")], row.names = FALSE, digits = 3)
cat("\nwrote:", file.path(res, "Fig6_placebo_controlled.png"), "\n")


