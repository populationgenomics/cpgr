# Generated by roxytest: Do not edit by hand!

# File R/utils.R: @testexamples

test_that("Function bucket_path() @ L47", {
  
  # Assuming that the analysis-runner has been invoked with
  # `--dataset tob-wgs --access-level test --output snp/v1`:
  Sys.setenv("DATASET" = "proj1", "ACCESS_LEVEL" = "test", "OUTPUT" = "dirA/v1")
  (b1 <- bucket_path("dir/v1/file.txt")) # gs://cpg-proj1-test/dir/v1/file.txt
  (b2 <- bucket_path("dir/v1/report.html", "web")) # gs://cpg-proj1-test-web/dir/v1/report.html
  
  expect_equal(b1, "gs://cpg-proj1-test/dir/v1/file.txt")
  expect_equal(b2, "gs://cpg-proj1-test-web/dir/v1/report.html")
})


test_that("Function output_path() @ L91", {
  
  # Assuming that the analysis-runner has been invoked with
  # `--dataset tob-wgs --access-level test --output snp/v1`:
  Sys.setenv("DATASET" = "proj1", "ACCESS_LEVEL" = "test", "OUTPUT" = "dirA/v1")
  (o1 <- output_path("report.html", "web")) # gs://cpg-proj1-test-web/dirA/v1/report.html
  (o2 <- output_path("output.txt")) # gs://cpg-proj1-test/dirA/v1/output.txt
  
  expect_equal(o1, "gs://cpg-proj1-test-web/dirA/v1/report.html")
  expect_equal(o2, "gs://cpg-proj1-test/dirA/v1/output.txt")
})

