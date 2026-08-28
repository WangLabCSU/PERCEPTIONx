# Verify the INSTALLED package (E:/rlibtest): S3 hijack removed, exports OK.
.libPaths(c('E:/rlibtest', .libPaths()))
suppressMessages(library(PERCEPTIONx))
fails <- 0L
chk <- function(name, ok, detail = "") {
  cat(sprintf("[%s] %s%s\n", if (ok) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0(" -> ", detail) else ""))
  if (!ok) fails <<- fails + 1L
}
set.seed(1); d <- data.frame(x = 1:10, y = 1:10 + rnorm(10))

# base cor.test must NOT be hijacked (returns full htest with conf.int)
r <- tryCatch(cor.test(~ x + y, data = d), error = function(e) e)
chk("installed: base cor.test not hijacked (conf.int)",
    inherits(r, "htest") && !is.null(r$conf.int))

# the trimmed generic's formula method still works
r2 <- tryCatch(cor.test_trimmed_v0(~ x + y, data = d), error = function(e) e)
chk("installed: cor.test_trimmed_v0.formula works",
    inherits(r2, "htest") && is.finite(r2$p.value))

# cor.test.formula must no longer be registered on the namespace
chk("installed: cor.test.formula NOT in namespace",
    !exists("cor.test.formula", envir = asNamespace("PERCEPTIONx"), inherits = FALSE))

# the app's background extract_meta worker depends on this internal function
chk("installed: internal extract_depmap_meta exists",
    exists("extract_depmap_meta", envir = asNamespace("PERCEPTIONx")))

# key exports
chk("installed: exports ok",
    all(c("load_depmap", "train_models", "prepare_data",
          "predict_drugs", "predict_patients") %in%
        getNamespaceExports("PERCEPTIONx")))

# callr declared in Imports
dsc <- packageDescription("PERCEPTIONx")
chk("installed: callr in Imports", grepl("callr", dsc$Imports, fixed = TRUE))

cat(if (fails == 0L) "\nINSTALL VERIFY ALL PASSED\n"
    else sprintf("\n%d CHECK(S) FAILED\n", fails))
quit(save = "no", status = if (fails == 0L) 0L else 1L)
