# =============================================================================
# read_10x(): one reader for CellRanger folders, h5 files and explicit triplets
# =============================================================================

# Write a minimal 10x triplet into `dir` and return the three paths.
write_triplet <- function(dir, gz = FALSE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  mtx <- c(
    "%%MatrixMarket matrix coordinate integer general",
    "%",
    "3 2 4",
    "1 1 1",
    "2 1 2",
    "3 2 4",
    "1 2 1"
  )
  feat <- c("ENSG1\tG1", "ENSG2\tG2", "ENSG3\tG2")
  bc <- c("CELL1", "CELL2")

  write_one <- function(lines, name) {
    path <- file.path(dir, name)
    con <- if (gz) gzfile(path, "wt") else file(path, "wt")
    writeLines(lines, con)
    close(con)
    path
  }
  list(
    mtx = write_one(mtx, if (gz) "matrix.mtx.gz" else "matrix.mtx"),
    barcodes = write_one(bc, if (gz) "barcodes.tsv.gz" else "barcodes.tsv"),
    features = write_one(feat, if (gz) "features.tsv.gz" else "features.tsv")
  )
}


test_that("read_10x reads a CellRanger folder", {
  dir <- tempfile("tenx_")
  write_triplet(dir)
  counts <- read_10x(dir)

  expect_equal(dim(counts), c(3L, 2L))
  expect_equal(colnames(counts), c("CELL1", "CELL2"))
  # Duplicate symbols are made unique, because Seurat requires unique rownames.
  expect_equal(rownames(counts), c("G1", "G2", "G2.1"))
  expect_equal(as.numeric(counts[, "CELL1"]), c(1, 2, 0))
  expect_equal(as.numeric(counts[, "CELL2"]), c(1, 0, 4))
  expect_s4_class(counts, "dgCMatrix")
})


test_that("read_10x reads gzip-compressed files", {
  dir <- tempfile("tenx_gz_")
  write_triplet(dir, gz = TRUE)
  counts <- read_10x(dir)

  expect_equal(dim(counts), c(3L, 2L))
  expect_equal(sum(counts), 8)
})


test_that("read_10x accepts the three files explicitly", {
  dir <- tempfile("tenx_paths_")
  paths <- write_triplet(dir)
  counts <- read_10x(
    mtx = paths$mtx,
    barcodes = paths$barcodes,
    features = paths$features
  )

  expect_equal(dim(counts), c(3L, 2L))
  expect_equal(colnames(counts), c("CELL1", "CELL2"))
})


test_that("read_10x aborts on an incomplete folder", {
  dir <- tempfile("tenx_missing_")
  dir.create(dir)
  writeLines("x", file.path(dir, "matrix.mtx"))

  expect_error(read_10x(dir), "No complete 10x triplet")
})


test_that("read_10x aborts on files that do not belong together", {
  dir <- tempfile("tenx_mismatch_")
  paths <- write_triplet(dir)
  writeLines(c("CELL1"), paths$barcodes)

  expect_error(
    read_10x(
      mtx = paths$mtx,
      barcodes = paths$barcodes,
      features = paths$features
    ),
    "do not belong together"
  )
})


test_that("read_10x rejects a path that is neither a folder nor .h5", {
  path <- tempfile(fileext = ".csv")
  writeLines("a,b", path)

  expect_error(read_10x(path), "neither a directory nor a .h5")
})
