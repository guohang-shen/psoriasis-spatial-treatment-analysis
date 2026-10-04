soft_file <- file.path("data", "GSE206391", "GSE206391_family.soft.gz")
out_file <- file.path("results", "GSE206391_sample_metadata.csv")

lines <- readLines(gzfile(soft_file), warn = FALSE)
starts <- grep("^\\^SAMPLE =", lines)
ends <- c(starts[-1] - 1L, length(lines))

get_value <- function(block, key) {
  hit <- grep(paste0("^!Sample_", key, " = "), block, value = TRUE)
  if (!length(hit)) return(NA_character_)
  sub(paste0("^!Sample_", key, " = "), "", hit[1])
}

parse_characteristics <- function(block, key) {
  hit <- grep("^!Sample_characteristics_ch1 =", block, value = TRUE)
  if (!length(hit)) return(NA_character_)
  vals <- sub("^!Sample_characteristics_ch1 = ", "", hit)
  one <- vals[grepl(paste0("^", key, ":"), vals)]
  if (!length(one)) return(NA_character_)
  sub(paste0("^", key, ":\\s*"), "", one[1])
}

rows <- lapply(seq_along(starts), function(i) {
  block <- lines[starts[i]:ends[i]]
  data.frame(
    gsm = sub("^\\^SAMPLE = ", "", block[1]),
    title = get_value(block, "title"),
    disease = parse_characteristics(block, "disease"),
    tissue_type = parse_characteristics(block, "tissue type"),
    patient = parse_characteristics(block, "patient"),
    stringsAsFactors = FALSE
  )
})

meta <- do.call(rbind, rows)
meta$sample_key <- sub(".*\\[([^]]+)\\].*", "\\1", meta$title)
meta$patient_primary <- sub(",.*", "", meta$patient)
meta$is_psoriasis <- grepl("psoriasis", meta$disease, ignore.case = TRUE)
meta$is_matched_psoriasis <- meta$is_psoriasis & grepl("lesional|non lesional", meta$tissue_type, ignore.case = TRUE)

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
write.csv(meta, out_file, row.names = FALSE, na = "")

cat("samples:", nrow(meta), "\n")
print(with(meta, table(disease, tissue_type, useNA = "ifany")))
cat("psoriasis matched samples:", sum(meta$is_matched_psoriasis), "\n")
cat("psoriasis patients:", length(unique(meta$patient_primary[meta$is_matched_psoriasis])), "\n")

