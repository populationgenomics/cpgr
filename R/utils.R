#' Make Directory
#'
#' Creates a directory.
#'
#' @param d Directory to create.
#'
#' @return If the directory exists, do nothing. If it doesn't, create it and
#' invisibly return a logical indicating if the operation succeeded.
#'
#' @export
mkdir <- function(d) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

#' Analysis-runner GCS bucket path
#'
#' Returns a full GCS path for the given bucket category and path.
#' This is useful for specifying input files, as in contrast to the
#' \code{\link{output_path}} function, `bucket_path` does _not_ take the "output"
#' parameter from the analysis-runner invocation into account.
#' Requires the `DATASET` and `ACCESS_LEVEL` environment variables to be set.
#'
#' @param path A path to append to the bucket.
#' @param bucket_category A category like "upload", "tmp", "web". If omitted,
#' defaults to the "main" and "test" buckets based on the access level. See
#' https://github.com/populationgenomics/team-docs/tree/main/storage_policies
#' for a full list of categories and their use cases.
#'
#' @return Full GCS path.
#'
#' @seealso \code{\link{output_path}}
#'
#' @examples
#' # Assuming that the analysis-runner has been invoked with
#' # `--dataset tob-wgs --access-level test --output snp/v1`:
#' Sys.setenv("DATASET" = "proj1", "ACCESS_LEVEL" = "test", "OUTPUT" = "dirA/v1")
#' bucket_path("dir/v1/file.txt") # gs://cpg-proj1-test/dir/v1/file.txt
#' bucket_path("dir/v1/report.html", "web") # gs://cpg-proj1-test-web/dir/v1/report.html
#'
#' @export
bucket_path <- function(path, bucket_category = NULL) {

  dataset <- Sys.getenv("DATASET")
  access_level <- Sys.getenv("ACCESS_LEVEL")
  assertthat::assert_that(nchar(dataset) > 0, nchar(access_level) > 0)

  namespace <- dplyr::if_else(access_level == "test", "test", "main")
  if (is.null(bucket_category)) {
    bucket_category <- namespace
  } else if (!bucket_category %in% c("archive", "upload")) {
    bucket_category <- glue::glue("{namespace}-{bucket_category}")
  }

  return(file.path("gs:/", glue::glue("cpg-{dataset}-{bucket_category}"), path))
}

#' Analysis-runner GCS output path
#'
#' In contrast to the \code{\link{bucket_path}} function, `output_path` takes
#' the "output" parameter from the analysis-runner invocation into account.
#'
#' @param path_suffix A suffix to append to the bucket + output directory.
#' @param bucket_category A category like "upload", "tmp", "web". If omitted,
#' defaults to the "main" and "test" buckets based on the access level. See
#' https://github.com/populationgenomics/team-docs/tree/main/storage_policies
#' for a full list of categories and their use cases.
#'
#' @return Full GCS path.
#'
#' @seealso \code{\link{bucket_path}}
#' @examples
#' # Assuming that the analysis-runner has been invoked with
#' # `--dataset tob-wgs --access-level test --output snp/v1`:
#' Sys.setenv("DATASET" = "proj1", "ACCESS_LEVEL" = "test", "OUTPUT" = "dirA/v1")
#' output_path("report.html", "web") # "gs://cpg-proj1-test-web/dirA/v1/report.html"
#' output_path("output.txt") # "gs://cpg-proj1-test/dirA/v1/output.txt"
#'
#' @export
output_path <- function(path_suffix, bucket_category = NULL) {
  output <- Sys.getenv("OUTPUT")
  assertthat::assert_that(nchar(output) > 0)
  return(bucket_path(file.path(output, path_suffix), bucket_category))
}
