library(tidyverse)
# for GCS authentication
library(googleCloudStorageR, include.only = c("gcs_list_objects", "gcs_auth"))
library(gargle, include.only = "token_fetch")
library(fs)
library(glue)

x1 <- here::here("nogit/cromwell/e5b51e30-3952-47f6-9b39-e17fb81cf929_metadata.json")
x2 <- here::here("nogit/cromwell/af805f60-ce68-495b-ab12-5886ca292e0e_metadata.json")

cromwell_read_metadata <- function(x) {
  j <- jsonlite::read_json(x)
  calls <- j$calls$GatherSampleEvidenceBatch.GatherSampleEvidence
  o <- purrr::map(calls, function(k) {
    m <- k[["subWorkflowMetadata"]][["calls"]]
    tibble::tibble(
      sample_id = k[["inputs"]][["sample_id"]],
      cram = m[["GatherSampleEvidence.CramToBam"]][[1]][["inputs"]][["cram_file"]],
      bam = m[["GatherSampleEvidence.CramToBam"]][[1]][["outputs"]][["bam_file"]],
      manta = m[["GatherSampleEvidence.Manta"]][[1]][["outputs"]][["vcf"]],
      delly = m[["GatherSampleEvidence.Delly"]][[1]][["outputs"]][["vcf"]],
      wham = m[["GatherSampleEvidence.Whamg"]][[1]][["outputs"]][["vcf"]],
      pesr_split = m[["GatherSampleEvidence.PESRCollection"]][[1]][["outputs"]][["split_out"]],
      pesr_disc = m[["GatherSampleEvidence.PESRCollection"]][[1]][["outputs"]][["disc_out"]],
      # sample_metrics = k[["outputs"]][["sample_metrics_files"]]
    )
  }) |>
    dplyr::bind_rows()
  o
}

m1 <- cromwell_read_metadata(x1)
m2 <- cromwell_read_metadata(x2)

m2 |>
  pivot_longer(cram:pesr_disc) |>
  filter(!grepl("gs://", value)) |>
  arrange(sample_id) |>
  select(sample_id, name)
  # pivot_wider(names_from = name, values_from=value)

# TODO: get size of cram/bam
m1 |>
  select(sample_id, cram, bam)

guess_file_type <- function(x) {
  case_when(
    grepl("\\.bam$", x, ignore.case = TRUE) ~ "BAM",
    grepl("\\.bai$", x, ignore.case = TRUE) ~ "BAMindex",
    grepl("\\.cram$", x, ignore.case = TRUE) ~ "CRAM",
    grepl("\\.crai$", x, ignore.case = TRUE) ~ "CRAMindex",
    grepl("\\.fastq.gz$", x, ignore.case = TRUE) ~ "FASTQ",
    grepl("\\.fastq$", x, ignore.case = TRUE) ~ "FASTQ",
    grepl("\\.fq$", x, ignore.case = TRUE) ~ "FASTQ",
    grepl("\\.fq\\.gz$", x, ignore.case = TRUE) ~ "FASTQ",
    grepl("manifest\\.txt$", x, ignore.case = TRUE) ~ "Manifest",
    grepl("\\.md5$", x, ignore.case = TRUE) ~ "MD5",
    grepl("md5\\.txt$", x, ignore.case = TRUE) ~ "MD5txt",
    grepl("\\.vcf$", x, ignore.case = TRUE) ~ "VCF_unz",
    grepl("\\.g\\.vcf\\.gz$", x, ignore.case = TRUE) ~ "GVCF",
    grepl("\\.vcf\\.gz$", x, ignore.case = TRUE) ~ "VCF",
    grepl("\\.tbi$", x, ignore.case = TRUE) ~ "VCFindex",
    grepl("\\.csv$", x, ignore.case = TRUE) ~ "CSV",
    grepl("\\.json$", x, ignore.case = TRUE) ~ "JSON",
    TRUE ~ "OTHER")
}

gcs_list_objects2 <- function(b, prefix) {
  googleCloudStorageR::gcs_list_objects(
    bucket = b,
    prefix = prefix,
    detail = "summary"
  ) %>%
    as_tibble() %>%
    mutate(
      name = glue("gs://{b}/{name}"),
      size = sub(" bytes", "", size), # else returns NA
      size = as_fs_bytes(size),
      ftype = guess_file_type(name)
    )
}


scope <- c("https://www.googleapis.com/auth/cloud-platform")
token <- token_fetch(scopes = scope)
gcs_auth(token = token)
cram_archive <- gcs_list_objects2("cpg-tob-wgs-archive", "cram")

crams <- cram_archive |>
  filter(ftype == "CRAM") |>
  mutate(basename = basename(name),
         sample = sub(".cram", "", basename),
         batch = sub("batch(.)", "\\1", basename(dirname(name)))) |>
  select(sample, fullname = name, batch, size, basename)

crams |>
  summarise(
    min = as_fs_bytes(min(size)),
    max = as_fs_bytes(max(size)),
    q1 = as_fs_bytes(quantile(size, 0.25)),
    median = as_fs_bytes(median(size)),
    q3 = as_fs_bytes(quantile(size, 0.75)),
    total = as_fs_bytes(sum(size)),
    .groups = "drop") %>%
  pivot_longer(cols = everything()) %>%
  knitr::kable(caption = glue("Size metrics for {nrow(crams)} CRAM files in archive."))


size_cutoff <- 40000000000 # ~30G

crams |>
  mutate(
    size_clean = str_trim(size),
    outlier = if_else(
      (!sample %in% failed_samples & batch == "14"),
      # size > size_cutoff,
      glue("{sample} ({size_clean}, batch{batch})"),
      NA_character_
    )
  ) %>%
  ggplot(aes(x = "", y = size, label = outlier, fill = batch, group = "")) +
  geom_violin(fill = "transparent", colour = "grey80", alpha = 0.4) +
  geom_sina(shape = 21, seed = 42, alpha = 0.7) +
  geom_text_repel(na.rm = TRUE, hjust = -0.6, seed = 42) +
  scale_y_continuous(labels = comma, breaks = breaks_pretty(10)) +
  theme(panel.grid.minor = element_blank()) +
  # facet_wrap(~batch) +
  ggtitle(glue("Size (bytes) for {nrow(crams)} CRAM files in 'archive'.")) +
  xlab("")


bams <- read_table(here::here("nogit/bam_sizes.txt"),
                   col_names = c("size", "gb", "name"),
                   col_types = "ccc") |>
  mutate(size = paste0(size, gb),
         size = fs::as_fs_bytes(size),
         basename = basename(name),
         sample = sub(".bam", "", basename)) |>
  select(sample, fullname = name, size, basename)

bams |>
  mutate(
    size_clean = str_trim(size),
    outlier = if_else(sample %in% failed_samples,
                      glue("{sample} ({size_clean})"),
                      NA_character_
    )) |>
  ggplot(aes(x = "", y = size, label = outlier, group = "")) +
  geom_violin(fill = "transparent", colour = "grey80", alpha = 0.4) +
  geom_sina(shape = 21, seed = 42, alpha = 0.7) +
  geom_text_repel(na.rm = TRUE, seed = 42) +
  scale_y_continuous(labels = comma, breaks = breaks_pretty(10)) +
  theme(panel.grid.minor = element_blank()) +
  # facet_wrap(~batch) +
  ggtitle(glue("Size (bytes) for {nrow(bams)} BAM files.")) +
  xlab("")

