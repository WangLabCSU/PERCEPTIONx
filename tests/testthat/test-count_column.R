# Regression test: a mapping that carries a per-clone abundance column
# (count / cells / n_cells / abundance) must expand to TRUE cell counts in
# prepare_data(), matching the Shiny app's upload behaviour. Before this was
# wired into prepare_data() the scripted path silently fell back to equal
# weights (1/n per clone), which changes patient-level scores.

test_that("prepare_data expands a mapping carrying a per-clone count column", {
  set.seed(1)
  genes  <- paste0("G", seq_len(20))
  clones <- c("P1_c1", "P1_c2", "P2_c1")
  expr <- matrix(runif(20 * 3), nrow = 20, dimnames = list(genes, clones))

  ct <- function(cc, pat, cl) {
    as.numeric(cc[cc$patients == pat, cl])
  }

  map_equal <- data.frame(cell_id = clones, patient_id = c("P1", "P1", "P2"),
                          stringsAsFactors = FALSE)
  map_count <- data.frame(cell_id = clones, patient_id = c("P1", "P1", "P2"),
                          count = c(90, 10, 50), stringsAsFactors = FALSE)

  p_eq <- suppressMessages(prepare_data(expression_matrix = expr,
                                        patient_mapping = map_equal,
                                        skip_clustering = TRUE))
  p_ct <- suppressMessages(prepare_data(expression_matrix = expr,
                                        patient_mapping = map_count,
                                        skip_clustering = TRUE))

  # without a count column every clone receives equal weight (1 cell)
  expect_equal(unname(ct(p_eq$clone_counts, "P1", c("c1", "c2"))), c(1, 1))

  # with a count column the true abundances are used
  expect_equal(unname(ct(p_ct$clone_counts, "P1", c("c1", "c2"))), c(90, 10))
  expect_equal(unname(ct(p_ct$clone_counts, "P2", "c1")), 50)

  # the mapping really was expanded to one row per cell
  expect_equal(nrow(p_ct$cell_clone_map), 150)

  # the abundance column must be consumed, not carried into the mapping
  expect_false("count" %in% colnames(p_ct$cell_clone_map))

  # alternate spellings are accepted
  map_alt <- data.frame(cell_id = clones, patient_id = c("P1", "P1", "P2"),
                        n_cells = c(90, 10, 50), stringsAsFactors = FALSE)
  p_alt <- suppressMessages(prepare_data(expression_matrix = expr,
                                         patient_mapping = map_alt,
                                         skip_clustering = TRUE))
  expect_equal(unname(ct(p_alt$clone_counts, "P1", c("c1", "c2"))), c(90, 10))
})
