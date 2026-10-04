## Patient-paired cross-disease negative control for locked psoriasis modules.
## The n here is patients, not spots or duplicate sections. These comparisons are
## exploratory because the module was originally selected in psoriasis data.

base_dir <- "C:/Users/15082/Desktop/Psoriasis"
res_dir <- file.path(base_dir, "results")
dat <- read.csv(file.path(res_dir, "GSE206391_cross_disease_paired_state_deltas.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)

exact_permutation_mean <- function(x, y) {
  z <- c(x, y)
  n <- length(x)
  observed <- mean(x) - mean(y)
  index <- combn(seq_along(z), n)
  null <- apply(index, 2, function(i) mean(z[i]) - mean(z[-i]))
  ## Inclusive exact p: observed allocation is part of the null distribution.
  c(mean_difference = observed,
    p_two_sided = mean(abs(null) >= abs(observed) - 1e-12),
    n_allocations = length(null))
}

states <- c("KC_IL17_IL36", "KC_IL17_response_no_keratins",
            "keratinocyte_activation", "keratinocyte_identity_clean")
rows <- list()
for (st in states) {
  pso <- dat[[st]][dat$disease == "psoriasis"]
  for (comparator in c("atopic dermatitis", "lichen planus")) {
    other <- dat[[st]][dat$disease == comparator]
    permutation <- exact_permutation_mean(pso, other)
    rows[[length(rows) + 1L]] <- data.frame(
      state = st, contrast = paste("psoriasis", "vs", comparator),
      n_psoriasis = length(pso), n_comparator = length(other),
      median_psoriasis = median(pso), median_comparator = median(other),
      mean_difference = unname(permutation["mean_difference"]),
      exact_p_two_sided = unname(permutation["p_two_sided"]),
      n_permutations = unname(permutation["n_allocations"]),
      stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(res_dir, "GSE206391_cross_disease_exact_permutation_v5.csv"),
          row.names = FALSE)
print(out, row.names = FALSE)

