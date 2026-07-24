# Quick test of demo data structure
sink("c:/Users/Lenovo/Desktop/PERCEPTION/tmp_test_output.txt")
set.seed(42)
gene_names <- c("TP53","BRCA1","EGFR","MYC","KRAS","PIK3CA","PTEN","RB1","APC","BRAF",
                "CDH1","CDKN2A","ERBB2","FGFR1","ALK","MET","RET","ROS1","NRAS","HRAS")
n_cells <- 200
n_patients <- 5
cell_names <- paste0("CELL_", sprintf("%03d", 1:n_cells))
patient_names <- paste0("PAT_", sprintf("%03d", 1:n_patients))
clinical_response <- data.frame(patient=patient_names,
                                response=c("Responder","Responder","Non-responder","Responder","Non-responder"),
                                stringsAsFactors=FALSE)
patient_assignment <- rep(patient_names, each=ceiling(n_cells/n_patients))[1:n_cells]
is_responder_cell <- clinical_response$response[match(patient_assignment, clinical_response$patient)] == "Responder"

abemaciclib_markers <- gene_names[1:5]
erlotinib_markers <- gene_names[6:10]
noise_genes <- gene_names[11:length(gene_names)]

expr_matrix <- matrix(0.1, nrow=length(gene_names), ncol=n_cells)
rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- cell_names

for (g in abemaciclib_markers) {
  expr_matrix[g, is_responder_cell] <- pmax(rnorm(sum(is_responder_cell), mean=12, sd=2), 0.1)
  expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean=2, sd=1), 0.1)
}
for (g in erlotinib_markers) {
  expr_matrix[g, is_responder_cell] <- pmax(rnorm(sum(is_responder_cell), mean=2, sd=1), 0.1)
  expr_matrix[g, !is_responder_cell] <- pmax(rnorm(sum(!is_responder_cell), mean=12, sd=2), 0.1)
}
for (g in noise_genes) {
  expr_matrix[g,] <- runif(n_cells, 0.5, 8)
}

cat("=== Demo data structure verification ===\n")
cat("abemaciclib markers (should be HIGH in responders, LOW in non-responders):\n")
for (g in abemaciclib_markers) {
  cat(sprintf("  %s: R=%.2f, NR=%.2f\n", g,
              mean(expr_matrix[g, is_responder_cell]),
              mean(expr_matrix[g, !is_responder_cell])))
}
cat("erlotinib markers (should be LOW in responders, HIGH in non-responders):\n")
for (g in erlotinib_markers) {
  cat(sprintf("  %s: R=%.2f, NR=%.2f\n", g,
              mean(expr_matrix[g, is_responder_cell]),
              mean(expr_matrix[g, !is_responder_cell])))
}
cat("noise gene example (should be similar between groups):\n")
cat(sprintf("  %s: R=%.2f, NR=%.2f\n", noise_genes[1],
            mean(expr_matrix[noise_genes[1], is_responder_cell]),
            mean(expr_matrix[noise_genes[1], !is_responder_cell])))

# Test training data generation
make_structured_training <- function(marker_genes, direction, n_train = 80) {
  x_train <- matrix(0, nrow=length(gene_names), ncol=n_train)
  rownames(x_train) <- gene_names
  half <- n_train %/% 2
  responder_like <- seq_len(half)
  nonresponder_like <- seq.int(half + 1, n_train)
  for (g in marker_genes) {
    x_train[g, responder_like] <- pmax(rnorm(half, mean=12, sd=2), 0.1)
    x_train[g, nonresponder_like] <- pmax(rnorm(half, mean=2, sd=1), 0.1)
  }
  for (g in setdiff(gene_names, marker_genes)) {
    x_train[g, ] <- runif(n_train, 0.5, 8)
  }
  y_train <- direction * colMeans(x_train[marker_genes, , drop=FALSE])
  y_train <- y_train + rnorm(n_train, sd=0.3)
  train_df <- as.data.frame(cbind(y=y_train, t(x_train)))
  train_df
}

# Verify training data structure
train_abema <- make_structured_training(abemaciclib_markers, +1, 80)
cat("\n=== Training data verification (abemaciclib) ===\n")
cat(sprintf("y for responder-like samples (high marker expr): mean=%.2f\n",
            mean(train_abema$y[1:40])))
cat(sprintf("y for non-responder-like samples (low marker expr): mean=%.2f\n",
            mean(train_abema$y[41:80])))
cat(sprintf("Column names preserved: %s\n",
            paste(head(colnames(train_abema), 3), collapse=", ")))

# Train a quick model
if (requireNamespace("caret", quietly=TRUE) && requireNamespace("glmnet", quietly=TRUE)) {
  library(caret)
  set.seed(101)
  model <- train(y ~ ., data=train_abema, method="glmnet",
                 trControl=trainControl(method="cv", number=3),
                 tuneLength=3)
  cat("\n=== Model trained successfully ===\n")
  cat("Model coefnames (first 5):", head(model$coefnames, 5), "\n")
  cat("Total features:", length(model$coefnames), "\n")

  # Test prediction on the demo expression
  # (Use cell-level expression directly for testing — not rank-normalized,
  # but should still show the signal pattern)
  pred <- predict(model, newdata=as.data.frame(t(expr_matrix)))
  cat("\n=== Prediction test ===\n")
  cat(sprintf("Predictions for responder cells: mean=%.2f, sd=%.2f\n",
              mean(pred[is_responder_cell]), sd(pred[is_responder_cell])))
  cat(sprintf("Predictions for non-responder cells: mean=%.2f, sd=%.2f\n",
              mean(pred[!is_responder_cell]), sd(pred[!is_responder_cell])))
  cat(sprintf("Prediction range: [%.2f, %.2f]\n", min(pred), max(pred)))
  if (mean(pred[is_responder_cell]) > mean(pred[!is_responder_cell])) {
    cat("SUCCESS: Predictions are HIGHER for responder cells (positive correlation with markers)\n")
  } else {
    cat("WARNING: Predictions do not show expected pattern\n")
  }
} else {
  cat("\n(caret/glmnet not available — skipping model training test)\n")
}
sink()
