# Verify key fixes from the full-repo audit. Run:
#   E:/R-4.6.1/bin/Rscript.exe dev/verify_fixes.R
# Exits non-zero if any check FAILS.

# Load the package source directly (function bodies are lazy, so sourcing in
# alphabetical order is fine).
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  suppressMessages(source(f))
}

fails <- 0L
chk <- function(name, ok, detail = "") {
  cat(sprintf("[%s] %s%s\n", if (ok) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0(" -> ", detail) else ""))
  if (!ok) fails <<- fails + 1L
}

# 1. range01 constant vector -> 0.5 (no NaN)
r <- range01(rep(5, 10))
chk("range01 constant vector returns 0.5", isTRUE(all(r == 0.5)))

# 2. range01 normal vector stays in [0,1]
set.seed(1); x <- rnorm(100)
r2 <- range01(x)
chk("range01 normal vector in [0,1]", all(r2 >= 0 & r2 <= 1, na.rm = TRUE))

# 3. hypergeometric both tails finite; lower.tail with full overlap ~ 1
g <- letters[1:26]
p_up <- hypergeometric_test_for_twolists(letters[1:5], letters[1:5], g, lower.tail = FALSE)
p_lo <- hypergeometric_test_for_twolists(letters[1:5], letters[1:5], g, lower.tail = TRUE)
chk("hypergeometric tails finite", is.finite(p_up) && is.finite(p_lo))
chk("hypergeometric lower.tail(full overlap) ~ 1", p_lo > 0.9)

# 4. cor.test_trimmed_v0 formula interface (S3) still works after rename
#    (cor.test formula convention: one-sided "~ x + y")
set.seed(2); xx <- rnorm(60); yy <- xx + rnorm(60)
dd <- data.frame(x = xx, y = yy)
r3 <- tryCatch(cor.test_trimmed_v0(~ x + y, data = dd), error = function(e) e)
chk("cor.test_trimmed_v0.formula dispatch", inherits(r3, "htest") && is.finite(r3$p.value))

# 5. base cor.test(~x+y) is NOT hijacked (full htest with conf.int)
r4 <- tryCatch(cor.test(~ x + y, data = dd), error = function(e) e)
chk("base cor.test keeps conf.int", inherits(r4, "htest") && !is.null(r4$conf.int))

# 6. clone_mean_expression aligns patient_ids when mapping covers fewer cells
expr <- matrix(rnorm(30), nrow = 10, ncol = 3,
               dimnames = list(letters[1:10], paste0("C", 1:3)))
cmap <- data.frame(cell_id = c("C1", "C2"), clone_id = c("cl1", "cl2"))
pids <- c(C1 = "P1", C2 = "P1", C3 = "P2")   # length == raw expr columns
res <- tryCatch(clone_mean_expression(expr, cmap, patient_ids = pids),
                error = function(e) e)
chk("clone_mean_expression aligns patient_ids", !inherits(res, "error"))

# 7. stripall2match keeps underscores (documented behavior)
chk("stripall2match keeps underscore", identical(stripall2match("Drug_1"), "drug_1"))

# 8. each_patient_pseudo_bulk empty-clone guard returns all-NA
cc <- data.frame(patients = c("P1", "P2"), cl1 = c(5, 0), cl2 = c(3, 0),
                 check.names = FALSE)
clz <- matrix(1, nrow = 4, ncol = 2,
              dimnames = list(NULL, c("P1@@cl1", "P1@@cl2")))
out <- tryCatch(each_patient_pseudo_bulk(x = 2, comb_viability_df = NULL,
                                         Clone_Counts_per_patients = cc,
                                         clone_Level_z_expression_df = clz),
                error = function(e) e)
chk("each_patient_pseudo_bulk empty-clone guard",
    !inherits(out, "error") && all(is.na(out)))

cat(if (fails == 0L) "\nALL CHECKS PASSED\n"
    else sprintf("\n%d CHECK(S) FAILED\n", fails))
quit(save = "no", status = if (fails == 0L) 0L else 1L)
