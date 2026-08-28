# verify_range01.R — verify the constant-vector fix: range01 returns 0
# ("not expressed") instead of the misleading 0.5 for sparse/constant genes.
# Usage: Rscript dev/verify_range01.R
suppressMessages(devtools::load_all(".", quiet = TRUE))

# 1. truly constant vector
r <- PERCEPTIONx::range01(rep(5, 10))
stopifnot(all(r == 0))
cat("constant vector -> all 0\n")

# 2. sparse gene (expressed in <=5% cells; 5th==95th percentile == 0)
sparse <- c(rep(0, 95), 5, 6, 7, 8, 9)
r <- PERCEPTIONx::range01(sparse)
stopifnot(length(r) == length(sparse), all(r == 0))
cat("sparse gene (95 zeros + 5 expressed) -> all 0, not 0.5\n")

# 3. normal gene still scales to [0,1] with min 0 / max 1
set.seed(1)
x <- rnorm(100) + 10
r <- PERCEPTIONx::range01(x)
stopifnot(all(r >= 0 & r <= 1), min(r) == 0, max(r) == 1)
cat("normal gene -> [0,1] scale, min 0 max 1\n")

# 4. demo genes: none should produce a CONSTANT 0.5 any more
set.seed(42)
gene_names <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PIK3CA", "PTEN", "RB1",
                "APC", "BRAF", "CDH1", "CDKN2A", "ERBB2", "FGFR1", "ALK",
                "MET", "RET", "ROS1", "NRAS", "HRAS", "MAP2K1", "MAPK1",
                "JAK2", "STAT3", "MTOR", "AKT1", "AKT2", "CTNNB1", "SMAD4",
                "VHL", "NF1", "NF2", "STK11", "FBXW7", "ARID1A", "KDM5C",
                "KMT2D", "SETD2", "BAP1", "PBRM1", "NOTCH1", "NOTCH2",
                "NOTCH3", "JAK1", "JAK3", "SOX9", "IDH1", "IDH2", "FLT3")
n_cells <- 400
cell_names <- paste0("CELL_", sprintf("%04d", 1:n_cells))
patient_names <- paste0("PAT_", sprintf("%03d", 1:20))
patient_assignment <- rep(patient_names, each = ceiling(n_cells / 20))[1:n_cells]
clinical_response <- data.frame(patient = patient_names,
                                response = c(rep("Responder", 10), rep("Non-responder", 10)),
                                stringsAsFactors = FALSE)
is_responder_cell <- clinical_response$response[
  match(patient_assignment, clinical_response$patient)] == "Responder"
expr_matrix <- matrix(0.1, nrow = length(gene_names), ncol = n_cells)
rownames(expr_matrix) <- gene_names; colnames(expr_matrix) <- cell_names
for (g in gene_names[1:5]) {
  expr_matrix[g, is_responder_cell]  <- pmax(rnorm(sum(is_responder_cell), mean = 8, sd = 3), 0.1)
  expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 3, sd = 2), 0.1)
}
for (g in gene_names[6:10]) {
  expr_matrix[g, is_responder_cell]  <- pmax(rnorm(sum(is_responder_cell), mean = 3, sd = 2), 0.1)
  expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean = 8, sd = 3), 0.1)
}
for (g in gene_names[11:length(gene_names)]) expr_matrix[g, ] <- runif(n_cells, 0.5, 8)
storage.mode(expr_matrix) <- "numeric"

flat05 <- character(0)
for (g in gene_names) {
  v <- PERCEPTIONx::range01(expr_matrix[g, ])
  if (all(v == 0.5)) flat05 <- c(flat05, g)
}
stopifnot(length(flat05) == 0)
cat("demo: no gene flattens to 0.5 anymore\n")
cat("\nRANGE01 VERIFY PASSED\n")
