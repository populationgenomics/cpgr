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
