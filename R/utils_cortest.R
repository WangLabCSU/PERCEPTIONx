# C_pKendall and C_pRho are internal C symbols from the stats package
# used via .Call() - declare as global variables to suppress R CMD check note
utils::globalVariables(c("C_pKendall", "C_pRho"))

#' Fast correlation test (generic)
#'
#' A generic function for fast correlation testing. Dispatches methods based on
#' the class of the first argument. This is a trimmed-down version of the base R
#' `cor.test` function, optimized for speed in large-scale computations.
#'
#' @param x First argument (typically a numeric vector).
#' @param ... Additional arguments passed to methods.
#'
#' @return A list containing the p-value and correlation estimate.
#'         The exact structure depends on the method used.
#'
#' @seealso \code{\link{cor.test}} for the full version with confidence intervals.
#'
#' @export
cor.test_trimmed_v0 <- function(x, ...) UseMethod("cor.test_trimmed_v0")


#' Fast Pearson correlation test (default method)
#'
#' Computes the Pearson correlation coefficient between two numeric vectors and
#' returns only the p-value and estimate. Confidence intervals and full
#' statistics are omitted to improve computational speed for large-scale
#' analyses (e.g., feature ranking across thousands of genes).
#'
#' Supports three correlation methods: Pearson, Kendall, and Spearman, with
#' one-sided or two-sided alternative hypotheses.
#'
#' @param x Numeric vector.
#' @param y Numeric vector of the same length as x.
#' @param alternative Alternative hypothesis: "two.sided", "less", or "greater".
#'        Default is "two.sided".
#' @param method Correlation method: "pearson", "kendall", or "spearman".
#'        Default is "pearson".
#' @param exact Logical; whether to compute exact p-values for Kendall and
#'        Spearman methods. Ignored for Pearson.
#' @param continuity Logical; whether to use continuity correction for Kendall
#'        and Spearman methods.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list of class "htest" containing:
#' \item{p.value}{The p-value of the test.}
#' \item{estimate}{The estimated correlation coefficient.}
#'
#' @seealso \code{\link{cor.test}} for the full version with confidence intervals.
#'
#' @keywords internal
#' @noRd
#' @export
cor.test_trimmed_v0.default <-
  function(x, y, alternative = c("two.sided", "less", "greater"),
           method = c("pearson"), exact = NULL,
           # conf.level = 0.95,
           continuity = FALSE, ...)
  {
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

    if(length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    if(!is.numeric(x)) stop("'x' must be a numeric vector")
    if(!is.numeric(y)) stop("'y' must be a numeric vector")
    OK <- complete.cases(x, y)
    x <- x[OK]
    y <- y[OK]
    n <- length(x)

    NVAL <- 0
    # conf.int <- FALSE

    if(method == "pearson") {
      if(n < 3L)
        stop("not enough finite observations")
      method <- "Pearson's product-moment correlation"
      names(NVAL) <- "correlation"
      r <- cor(x, y)
      df <- n - 2L
      ESTIMATE <- c(cor = r)
      PARAMETER <- c(df = df)
      STATISTIC <- c(t = sqrt(df) * r / sqrt(1 - r^2))
      # Do not compute confidence Int
      # if(n > 3) { ## confidence int.
      #   if(!missing(conf.level) &&
      #      (length(conf.level) != 1 || !is.finite(conf.level) ||
      #       conf.level < 0 || conf.level > 1))
      #     stop("'conf.level' must be a single number between 0 and 1")
      #   conf.int <- TRUE
      #   z <- atanh(r)
      #   sigma <- 1 / sqrt(n - 3)
      #   cint <-
      #     switch(alternative,
      #            less = c(-Inf, z + sigma * qnorm(conf.level)),
      #            greater = c(z - sigma * qnorm(conf.level), Inf),
      #            two.sided = z +
      #              c(-1, 1) * sigma * qnorm((1 + conf.level) / 2))
      #   cint <- tanh(cint)
      #   attr(cint, "conf.level") <- conf.level
      # }
      PVAL <- switch(alternative,
                     "less" = pt(STATISTIC, df),
                     "greater" = pt(STATISTIC, df, lower.tail=FALSE),
                     "two.sided" = 2 * min(pt(STATISTIC, df),
                                           pt(STATISTIC, df, lower.tail=FALSE)))
    }
    else {
      if(n < 2)
        stop("not enough finite observations")
      PARAMETER <- NULL
      TIES <- (min(length(unique(x)), length(unique(y))) < n)
      if(method == "kendall") {
        method <- "Kendall's rank correlation tau"
        names(NVAL) <- "tau"
        r <- cor(x,y, method = "kendall")
        ESTIMATE <- c(tau = r)

        if(!is.finite(ESTIMATE)) {  # all x or all y the same
          ESTIMATE[] <- NA
          STATISTIC <- c(T = NA)
          PVAL <- NA
        }
        else {
          if(is.null(exact))
            exact <- (n < 50)
          if(exact && !TIES) {
            q <- round((r + 1) * n * (n - 1) / 4)
            STATISTIC <- c(T = q)
            pkendall <- function(q, n) .Call(C_pKendall, q, n)
            PVAL <-
              switch(alternative,
                     "two.sided" = {
                       if(q > n * (n - 1) / 4)
                         p <- 1 - pkendall(q - 1, n)
                       else
                         p <- pkendall(q, n)
                       min(2 * p, 1)
                     },
                     "greater" = 1 - pkendall(q - 1, n),
                     "less" = pkendall(q, n))
          } else {
            xties <- table(x[duplicated(x)]) + 1
            yties <- table(y[duplicated(y)]) + 1
            T0 <- n * (n - 1)/2
            T1 <- sum(xties * (xties - 1))/2
            T2 <- sum(yties * (yties - 1))/2
            S <- r * sqrt((T0 - T1) * (T0 - T2))
            v0 <- n * (n - 1) * (2 * n + 5)
            vt <- sum(xties * (xties - 1) * (2 * xties + 5))
            vu <- sum(yties * (yties - 1) * (2 * yties + 5))
            v1 <- sum(xties * (xties - 1)) * sum(yties * (yties - 1))
            v2 <- sum(xties * (xties - 1) * (xties - 2)) *
              sum(yties * (yties - 1) * (yties - 2))

            var_S <- (v0 - vt - vu) / 18 +
              v1 / (2 * n * (n - 1)) +
              v2 / (9 * n * (n - 1) * (n - 2))

            if(exact && TIES)
              warning("Cannot compute exact p-value with ties")
            if (continuity) S <- sign(S) * (abs(S) - 1)
            STATISTIC <- c(z = S / sqrt(var_S))
            PVAL <- switch(alternative,
                           "less" = pnorm(STATISTIC),
                           "greater" = pnorm(STATISTIC, lower.tail=FALSE),
                           "two.sided" = 2 * min(pnorm(STATISTIC),
                                                 pnorm(STATISTIC, lower.tail=FALSE)))
          }
        }
      } else {
        method <- "Spearman's rank correlation rho"
        if (is.null(exact))
          exact <- TRUE
        names(NVAL) <- "rho"
        r <- cor(rank(x), rank(y))
        ESTIMATE <- c(rho = r)
        if(!is.finite(ESTIMATE)) {  # all x or all y the same
          ESTIMATE[] <- NA
          STATISTIC <- c(S = NA)
          PVAL <- NA
        }
        else {
          ## Use the test statistic S = sum(rank(x) - rank(y))^2
          ## and AS 89 for obtaining better p-values than via the
          ## simple normal approximation.
          ## In the case of no ties, S = (1-rho) * (n^3-n)/6.
          pspearman <- function(q, n, lower.tail = TRUE) {
            if(n <= 1290 && exact) # n*(n^2 - 1) does not overflow
              .Call(C_pRho, round(q) + 2*lower.tail, n, lower.tail)
            else { # for large n: asymptotic t_{n-2}
              den <- (n*(n^2-1))/6 # careful for overflow
              ## Kendall et all (1939) p. 260
              if (continuity) den <- den + 1
              r <- 1 - q/den
              pt(r / sqrt((1 - r^2)/(n-2)), df = n-2,
                 lower.tail = !lower.tail)
            }
          }
          q <- (n^3 - n) * (1 - r) / 6
          STATISTIC <- c(S = q)
          if(TIES && exact){
            exact <- FALSE
            warning("Cannot compute exact p-value with ties")
          }
          PVAL <-
            switch(alternative,
                   "two.sided" = {
                     p <- if(q > (n^3 - n) / 6)
                       pspearman(q, n, lower.tail = FALSE)
                     else
                       pspearman(q, n, lower.tail = TRUE)
                     min(2 * p, 1)
                   },
                   "greater" = pspearman(q, n, lower.tail = TRUE),
                   "less" = pspearman(q, n, lower.tail = FALSE))
        }
      }
    }

    RVAL <- list(
      # statistic = STATISTIC,
      # parameter = PARAMETER,
      p.value = as.numeric(PVAL),
      estimate = ESTIMATE
      # ,
      # null.value = NVAL,
      # alternative = alternative,
      # method = method,
      # data.name = DNAME
    )
    # if(conf.int)
    #   RVAL <- c(RVAL, list(conf.int = cint))
    class(RVAL) <- "htest"
    RVAL
  }



#' Formula interface for fast correlation test
#'
#' Provides a formula interface to \code{\link{cor.test_trimmed_v0.default}},
#' allowing users to specify variables using R formula syntax.
#'
#' @param formula A formula of the form ~ x + y or ~ x + y + z (only first two
#'        variables are used). The response is not used; only the right-hand side
#'        is evaluated.
#' @param data A data frame or environment containing the variables in the formula.
#' @param subset An optional vector specifying a subset of observations to use.
#' @param na.action A function specifying how to handle missing values.
#' @param ... Additional arguments passed to \code{\link{cor.test_trimmed_v0.default}}.
#'
#' @return A list of class "htest" containing the p-value and correlation estimate.
#'
#' @seealso \code{\link{cor.test_trimmed_v0.default}}
#'
#' @keywords internal
#' @noRd
#' @exportS3Method stats::cor.test
cor.test.formula <-
  function(formula, data, subset, na.action, ...)
  {
    if(missing(formula)
       || !inherits(formula, "formula")
       || length(formula) != 2L)
      stop("'formula' missing or invalid")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, environment(formula))
    if(length(mf) != 2L)
      stop("invalid formula")
    DNAME <- paste(names(mf), collapse = " and ")
    names(mf) <- c("x", "y")
    y <- do.call("cor.test_trimmed_v0", c(mf, list(...)))
    y$data.name <- DNAME
    y
  }
