library(tidyverse)
# for GCS authentication
# library(googleCloudStorageR, include.only = c("gcs_list_objects", "gcs_auth"))
# library(gargle, include.only = "token_fetch")
library(fs)
library(glue)
library(here)

x1 <- here("nogit/cromwell/01_e5b51e30-3952-47f6-9b39-e17fb81cf929_metadata.json")
x2 <- here("nogit/cromwell/02_af805f60-ce68-495b-ab12-5886ca292e0e_metadata.json")
x3 <- here("nogit/cromwell/03_e6aeef3c-878f-4df3-b257-4fd6eebaf7d6_metadata.json")
x4 <- here("nogit/cromwell/04_04da63d0-03a5-431d-8ec2-449cada8f647_metadata.json")

cromwell_read_metadata <- function(x) {
  j <- jsonlite::read_json(x)
  calls <- j$calls$GatherSampleEvidenceBatch.GatherSampleEvidence
  clean_sample_metrics <- function(d) {
    if (nrow(d) == 0) {
      return(list(metrics = NA))
    }
    out <- d |>
      mutate(
        b = basename(value),
        group = case_when(
          grepl("delly_.*\\.vcf\\.tsv", b) ~ "metrics_delly",
          grepl("manta_.*\\.vcf\\.tsv", b) ~ "metrics_manta",
          grepl("wham_.*\\.vcf\\.tsv", b) ~ "metrics_wham",
          grepl(".*\\.sr-file\\.tsv", b) ~ "metrics_sr",
          grepl(".*\\.pe-file\\.tsv", b) ~ "metrics_pe",
          grepl(".*\\.raw-counts\\.tsv", b) ~ "metrics_raw"
        )) |>
      select(group, value) |>
      pivot_wider(names_from = group, values_from = value)
    list(metrics = out)
  }
  o <- purrr::map(calls, function(k) {
    m <- k[["subWorkflowMetadata"]][["calls"]]
    tibble::tibble(
      sample_id = k[["inputs"]][["sample_id"]],
      cram = m[["GatherSampleEvidence.CramToBam"]][[1]][["inputs"]][["cram_file"]],
      bam = m[["GatherSampleEvidence.CramToBam"]][[1]][["outputs"]][["bam_file"]],
      bami = m[["GatherSampleEvidence.CramToBam"]][[1]][["outputs"]][["bam_index"]],
      manta = m[["GatherSampleEvidence.Manta"]][[1]][["outputs"]][["vcf"]],
      mantai = m[["GatherSampleEvidence.Manta"]][[1]][["outputs"]][["index"]],
      delly = m[["GatherSampleEvidence.Delly"]][[1]][["outputs"]][["vcf"]],
      dellyi = m[["GatherSampleEvidence.Delly"]][[1]][["outputs"]][["index"]],
      wham = m[["GatherSampleEvidence.Whamg"]][[1]][["outputs"]][["vcf"]],
      whami = m[["GatherSampleEvidence.Whamg"]][[1]][["outputs"]][["index"]],
      counts = m[["GatherSampleEvidence.CollectCounts"]][[1]][["outputs"]][["counts"]],
      pesr_split = m[["GatherSampleEvidence.PESRCollection"]][[1]][["outputs"]][["split_out"]],
      pesr_spliti = m[["GatherSampleEvidence.PESRCollection"]][[1]][["outputs"]][["split_out_index"]],
      pesr_disc = m[["GatherSampleEvidence.PESRCollection"]][[1]][["outputs"]][["disc_out"]],
      pesr_disci = m[["GatherSampleEvidence.PESRCollection"]][[1]][["outputs"]][["disc_out_index"]],
      sample_metrics = m[["GatherSampleEvidence.GatherSampleEvidenceMetrics"]][[1]][["outputs"]][["sample_metrics_files"]] |> unlist() |> as_tibble_col() |> clean_sample_metrics()
    )
  }) |>
    dplyr::bind_rows() |>
    tidyr::unnest(sample_metrics)
    # select(-sample_metrics)
  o
}

m1 <- cromwell_read_metadata(x1) |> select(-sample_metrics)
m2 <- cromwell_read_metadata(x2) |> select(-sample_metrics)
m3 <- cromwell_read_metadata(x3)
m4 <- cromwell_read_metadata(x4)

m1 |>
  pivot_longer(cram:metrics_raw) |>
  filter(!grepl("gs://", value)) |>
  filter(!grepl("metrics|whami|dellyi", name)) |>
  arrange(sample_id) |>
  select(sample_id, name)
  # pivot_wider(names_from = name, values_from=value)

# m1 |>
#   select(-cram) |>
#   pivot_longer(bam:metrics_raw) |>
#   filter(!is.na(value)) |>
#   mutate(name2 = case_when(
#     grepl("^metrics", name) ~ "svmetrics",
#     grepl("bam", name) ~ "bam",
#     grepl("^manta", name) ~ "manta",
#     grepl("^delly", name) ~ "delly",
#     grepl("^wham", name) ~ "wham",
#     grepl("^pesr", name) ~ "pesr",
#     grepl("counts", name) ~ "covcounts",
#     TRUE ~ "other"
#   )) |>
#   select(sample_id, name = name2, from = value) |>
#   mutate(to = glue("gs://cpg-tob-wgs-test/pdiakumis/gatksv/{sample_id}/{name}/")) |>
#   mutate(cmd = glue("gsutil mv {from} {to}")) |>
#   select(cmd) |>
#   write_tsv("nogit/transfer_outputs.sh", col_names = FALSE)

bind_rows(m3, m4) |>
  select(-cram) |>
  pivot_longer(bam:metrics_raw) |>
  filter(!is.na(value)) |>
  mutate(name2 = case_when(
    grepl("^metrics", name) ~ "svmetrics",
    grepl("bam", name) ~ "bam",
    grepl("^manta", name) ~ "manta",
    grepl("^delly", name) ~ "delly",
    grepl("^wham", name) ~ "wham",
    grepl("^pesr", name) ~ "pesr",
    grepl("counts", name) ~ "covcounts",
    TRUE ~ "other"
  )) |>
  select(sample_id, name = name2, from = value) |>
  filter(name != "bam") |>
  mutate(to = glue("gs://cpg-tob-wgs-test/pdiakumis/gatksv/{sample_id}/{name}/")) |>
  mutate(cmd = glue("gsutil mv {from} {to}")) |>
  select(cmd) |>
  write_tsv("nogit/transfer_outputs2.sh", col_names = FALSE)

m2 |> filter(sample_id == "CPG11049") |>
  select(-cram) |>
  pivot_longer(bam:metrics_raw) |>
  mutate(name2 = case_when(
    grepl("^metrics", name) ~ "svmetrics",
    grepl("bam", name) ~ "bam",
    grepl("^manta", name) ~ "manta",
    grepl("^delly", name) ~ "delly",
    grepl("^wham", name) ~ "wham",
    grepl("^pesr", name) ~ "pesr",
    grepl("counts", name) ~ "covcounts",
    TRUE ~ "other"
  )) |>
  select(sample_id, name = name2, from = value) |>
  filter(name %in% c("delly", "svmetrics")) |>
  mutate(to = glue("gs://cpg-tob-wgs-test/pdiakumis/gatksv/{sample_id}/{name}/")) |>
  mutate(cmd = glue("gsutil mv {from} {to}")) |>
  select(cmd) |>
  write_tsv("nogit/transfer_outputs3.sh", col_names = FALSE)

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


read_tsv("nogit/cromwell/bucket_contents.txt", col_types = "c", col_names = "path") |>
  separate(path, into = c("sample_id", "group", "fname"), sep = "/") |>
  count(sample_id) |> View()



x <- here("nogit/cromwell/melt_a360b090-fcee-4118-b604-8e9b4c3b7ca7_metadata.json")
cromwell_read_melt_metadata <- function(x) {
  j <- jsonlite::read_json(x)
  calls <- j$calls$GatherSampleEvidenceBatch.GatherSampleEvidence
  o <- purrr::map(calls, function(k) {
    m <- k[["subWorkflowMetadata"]][["calls"]]
    melt <- m[["GatherSampleEvidence.MELT"]][[1]][["subWorkflowMetadata"]][["calls"]][["MELT.RunMELT"]][[1]]
    metrics <- m[["GatherSampleEvidence.GatherSampleEvidenceMetrics"]][[1]][["outputs"]][["sample_metrics_files"]]
    tibble::tibble(
      sample_id = k[["inputs"]][["sample_id"]],
      melt_read_length = melt[["inputs"]][["read_length"]],
      melt_insert_size = melt[["inputs"]][["insert_size"]],
      melt_coverage = melt[["inputs"]][["coverage"]],
      melt_chimeras = melt[["inputs"]][["pct_chimeras"]],
      melt_total_reads = melt[["inputs"]][["total_reads"]],
      melt_pf_reads_improper_pairs = melt[["inputs"]][["pf_reads_improper_pairs"]],
      melt_metrics = metrics |> unlist() |> grep(pattern = "melt_.*\\.vcf\\.tsv", value = TRUE),
      melt_vcf = melt[["outputs"]][["vcf"]],
      melt_vcfi = melt[["outputs"]][["index"]]
    )
  }) |>
     dplyr::bind_rows() #|>
    # tidyr::unnest(sample_metrics)
    # select(-sample_metrics)
  o
}

y <- cromwell_read_melt_metadata(x)

y |>
  select(sample_id, melt_metrics, melt_vcf, melt_vcfi) |>
  pivot_longer(c(melt_metrics, melt_vcf, melt_vcfi), values_to = "from") |>
  mutate(name = if_else(grepl("melt_metrics", name), "svmetrics", "melt")) |>
  mutate(to = glue("gs://cpg-tob-wgs-test/pdiakumis/gatksv/{sample_id}/{name}/")) |>
  mutate(cmd = glue("gsutil mv {from} {to}")) |>
  select(cmd) |>
  write_tsv("nogit/transfer_outputs_melt.sh", col_names = FALSE)

read_tsv("nogit/cromwell/bucket_contents2.txt", col_types = "c", col_names = "path") |>
  separate(path, into = c("gs", "foo", "bucket", "me", "gatksv", "sample_id", "group", "fname"), sep = "/") |>
  count(sample_id) |>
  filter(n != 22) |>
  pull(sample_id) |>
  sort()

