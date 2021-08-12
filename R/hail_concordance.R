#' Convert Hail Concordance Raw Results to Matrix
#'
#' Converts raw string with Hail concordance results to a 5x5 matrix.
#'
#' @param x String containing a JSON array of 5 elements, each with 5 subelements.
#' @param rowname Name suffix for matrix rows.
#' @param colname Name suffix for matrix columns.
#'
#' @return A 5x5 matrix with Hail concordance results for a single sample.
#'
#' @examples
#' f1 <- system.file("extdata/hail/concordance/single.tsv", package = "cpgr")
#' s <- readr::read_tsv(f1, col_types = "ccd") |> dplyr::pull(concordance)
#' (m <- hail_conc_raw2matrix(s))
#'
#' @testexamples
#' expect_equal(m["het_left", "het_right"], 1394)
#' expect_equal(m["missing_variant_left", "homalt_right"], 21271)
#'
#' @export
hail_conc_raw2matrix <- function(x, rowname = "left", colname = "right") {
  m <- jsonlite::fromJSON(x)
  assertthat::assert_that(ncol(m) == 5, nrow(m) == 5, is.numeric(m))
  conc_names <- c("missing_variant", "missing_gt", "homref", "het", "homalt")
  colnames(m) <- glue::glue("{conc_names}_{colname}")
  rownames(m) <- glue::glue("{conc_names}_{rowname}")
  m
}

#' Summarise Hail concordance results
#'
#' Summarises Hail concordance results.
#'
#' @param m Named matrix (5x5) with Hail concordance results for a single sample.
#'
#' @return
#' @export
#'
#' @examples
hail_conc_summary <- function(m, colname = "right", rowname = "left") {
  conc_tot1 <- sum(m[3, 3], m[4, 4], m[5, 5])
  disc_tot1 <- sum(m[3:5, 3:5]) - conc_tot1
  disc_tot2 <- m["het_snp", "missing_variant_wgs"] + m["homalt_snp", "missing_variant_wgs"]
  disc_tot <- disc_tot1 + disc_tot2

  # missing, homrefs, hets and homalts
  snp_tot <- sum(m[c(2, 3, 4, 5), ])
  wgs_tot <- sum(m[, c(2, 3, 4, 5)])
  conc_over_snp <- conc_tot / snp_tot
  conc_over_wgs <- conc_tot / wgs_tot
  disc_over_snp <- disc_tot / snp_tot
  disc_over_wgs <- disc_tot / wgs_tot
  tibble(disc_tot = disc_tot,
         conc_tot = conc_tot,
         wgs_tot = wgs_tot,
         snp_tot = snp_tot,
         `conc / wgs` = conc_over_wgs,
         `conc / snp` = conc_over_snp,
         `disc / snp` = disc_over_snp,
         `disc / wgs` = disc_over_wgs)
}


#
#hail_concordance_summarise <- function(x) {
#
#  x |>
#    readr::read_tsv(col_types = "ccd") |>
#    dplyr::mutate(conc = purrr::map(concordance, ~jsonlite::fromJSON(.) |> mat_rename())) |>
#    tidyr::unnest(conc)
#}
##
###res <- params$RES_SAMPLES |>
###  read_tsv(col_types = "ccd") |>
###  mutate(conc2 = map(concordance, ~fromJSON(.) |> get_mat_stats())) |>
###  unnest(conc2) |>
###  select(sample = s, `conc / snp`, `conc / wgs`,
###         conc_tot, disc_tot,
###         snp_tot, wgs_tot,
###         `disc / snp`, `disc / wgs`)
###
