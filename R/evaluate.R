#' Evaluate the fidelity of a compressed posterior
#'
#' Compares draws regenerated from a compressed posterior against a
#' reference set of draws (typically the original MCMC output) using
#' two distribution-free, backend-independent two-sample diagnostics:
#'
#' \itemize{
#'   \item \strong{Energy distance} (Szekely & Rizzo, 2013): a
#'         scale-invariant distance between two samples that is zero
#'         iff their distributions match. The raw distance is anchored
#'         against a \emph{self-baseline} obtained by bootstrapping
#'         `reference_draws` to two independent samples matched in
#'         size to the reconstructed sample. The 90th percentile of
#'         this baseline distribution defines the Monte Carlo noise
#'         envelope; the reproduction score is
#'         `100 * min(1, noise_envelope / max(distance, noise_envelope))`,
#'         so 100\% means the reconstruction is within the bootstrap
#'         envelope and the score decays as 1/ratio for larger
#'         distances.
#'   \item \strong{Classifier two-sample test (C2ST)} (Lopez-Paz &
#'         Oquab, 2017): trains a classifier (random forest via
#'         `ranger::ranger()` under cross-validation
#'         to discriminate reference draws from reconstructed draws and
#'         reports the out-of-fold ROC AUC. AUC = 0.5 means the two
#'         samples are indistinguishable; the score is mapped to
#'         `100 * (1 - 2 * |AUC - 0.5|)`.
#' }
#'
#' Both metrics are computed only from samples and are independent of
#' the compression backend (`mclust`, `mvdens_gmm`, `mvdens_kde`), so
#' the score is not "circular": it does not depend on the same density
#' model that produced the compression.
#'
#' For a strict out-of-sample evaluation, hold out a fraction of the
#' draws \emph{before} fitting:
#'
#' \preformatted{
#'   idx  <- sample.int(nrow(draws), size = 0.8 * nrow(draws))
#'   comp <- compress_posterior(draws[idx, ], method = "mclust")
#'   evaluate_compression(comp, reference_draws = draws[-idx, ])
#' }
#'
#' @param comp A `posterior_compressed` object (or a path to an `.rds`
#'   file containing one), or the return value of [compress_brmsfit()] /
#'   [compress_sccomp()] (`compressed_fit`: `$compressed` is used).
#' @param reference_draws A draws matrix, data.frame, or
#'   `posterior::draws_*` object whose columns include all parameters
#'   in `comp`. Treated as the ground-truth distribution.
#' @param metric Character vector of metrics to compute. One or more of
#'   `"energy"` and `"c2st"`. Default: both.
#' @param n_eval Integer number of draws to regenerate from `comp`.
#'   Defaults to `min(max_n, nrow(reference_draws))`.
#' @param max_n Cap on the number of points used in any pairwise
#'   distance / classifier computation. Both samples are subsampled to
#'   at most this many rows. Default `2000`.
#' @param n_self_reps Number of self-baseline replicates for the energy
#'   metric. Default `20`.
#' @param classifier Classifier for C2ST. `"ranger"` (default) uses a
#'   random forest (`ranger::ranger()`). `"knn"` uses a k-NN probability
#'   estimate instead (no random forest).
#' @param cv_folds Cross-validation folds for C2ST. Default `5`.
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical; print progress messages.
#'
#' @return An S3 object of class `compression_fidelity` containing:
#' \describe{
#'   \item{`reproduction_pct`}{Headline reproduction score in `[0, 100]`,
#'     averaged across requested metrics.}
#'   \item{`metrics`}{Named list with detailed per-metric results.}
#'   \item{`n_reference`, `n_eval`, `n_params`}{Sample sizes used.}
#' }
#'
#' @references
#' Szekely, G. J. & Rizzo, M. L. (2013). Energy statistics: A class of
#' statistics based on distances. *Journal of Statistical Planning and
#' Inference*, 143(8), 1249-1272.
#'
#' Lopez-Paz, D. & Oquab, M. (2017). Revisiting Classifier Two-Sample
#' Tests. *ICLR*.
#'
#' @examples
#' set.seed(1)
#' draws <- matrix(rnorm(2000 * 3), ncol = 3,
#'                 dimnames = list(NULL, c("alpha", "beta", "sigma")))
#' comp <- compress_posterior(draws, method = "mclust", n_components = 2)
#' fidelity <- evaluate_compression(comp, reference_draws = draws, seed = 1L)
#' fidelity
#'
#' @importFrom stats predict
#'
#' @export
evaluate_compression <- function(
    comp,
    reference_draws,
    metric       = c("energy", "c2st"),
    n_eval       = NULL,
    max_n        = 2000L,
    n_self_reps  = 20L,
    classifier   = c("ranger", "knn"),
    cv_folds     = 5L,
    seed         = NULL,
    verbose      = FALSE) {
  comp <- .resolve_compressed(comp)
  if (!inherits(comp, "posterior_compressed")) {
    stop(
      "`comp` must be a posterior_compressed object, a path to one, ",
      "or a compressed_fit from compress_brmsfit() / compress_sccomp().",
      call. = FALSE
    )
  }

  metric     <- match.arg(metric, choices = c("energy", "c2st"),
                          several.ok = TRUE)
  classifier <- match.arg(classifier)

  ref_mat <- .as_draws_matrix(reference_draws, variables = comp$param_names)
  if (ncol(ref_mat) == 0L) {
    stop("`reference_draws` has no columns matching the compressed posterior.",
         call. = FALSE)
  }
  if (nrow(ref_mat) < 4L) {
    stop("`reference_draws` must contain at least 4 rows for evaluation.",
         call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  if (is.null(n_eval)) n_eval <- min(max_n, nrow(ref_mat))
  if (verbose) {
    .message("Generating ", n_eval, " draws from compressed posterior ...")
  }
  recon_mat <- sample_posterior(comp, n_draws = n_eval)
  recon_mat <- recon_mat[, comp$param_names, drop = FALSE]

  ref_used   <- .subsample_rows(ref_mat,   max_n)
  recon_used <- .subsample_rows(recon_mat, max_n)

  scaling      <- .scaling_from(ref_used)
  ref_scaled   <- .apply_scaling(ref_used,   scaling)
  recon_scaled <- .apply_scaling(recon_used, scaling)

  metrics_out <- list()

  if ("energy" %in% metric) {
    if (verbose) .message("Computing energy distance ...")
    metrics_out$energy <- .energy_metric(
      ref_scaled, recon_scaled, n_self_reps = n_self_reps
    )
  }
  if ("c2st" %in% metric) {
    if (verbose) .message("Computing classifier two-sample test ...")
    metrics_out$c2st <- .c2st_metric(
      ref_scaled, recon_scaled,
      classifier = classifier,
      cv_folds   = cv_folds
    )
  }

  pcts <- vapply(metrics_out, function(m) m$reproduction_pct, numeric(1))
  reproduction_pct <- if (length(pcts) > 0L) mean(pcts) else NA_real_

  out <- list(
    reproduction_pct = reproduction_pct,
    metrics          = metrics_out,
    n_reference      = nrow(ref_used),
    n_eval           = nrow(recon_used),
    n_params         = ncol(ref_used),
    method           = comp$method,
    param_names      = comp$param_names
  )
  class(out) <- c("compression_fidelity", "list")
  out
}


#' @rdname evaluate_compression
#' @param x A `compression_fidelity` object.
#' @param ... Currently unused.
#' @export
print.compression_fidelity <- function(x, ...) {
  cat("<compression_fidelity>\n")
  cat(sprintf(
    "  method        : %s\n  parameters    : %d\n  reference n   : %d\n  eval n        : %d\n",
    x$method, x$n_params, x$n_reference, x$n_eval
  ))
  cat("  ----------------------------------------\n")
  if (!is.null(x$metrics$energy)) {
    e <- x$metrics$energy
    cat(sprintf(
      "  energy        : %5.1f%% reproduction\n",
      e$reproduction_pct
    ))
    cat(sprintf(
      "    distance    : %.4f   noise envelope (q90): %.4f   ratio: %.2fx\n",
      e$energy, e$noise_envelope, e$ratio
    ))
  }
  if (!is.null(x$metrics$c2st)) {
    c2 <- x$metrics$c2st
    cat(sprintf(
      "  C2ST          : %5.1f%% reproduction\n",
      c2$reproduction_pct
    ))
    cat(sprintf(
      "    AUC         : %.3f   classifier: %s   cv_folds: %d\n",
      c2$auc, c2$classifier, c2$cv_folds
    ))
  }
  cat("  ----------------------------------------\n")
  cat(sprintf("  reproduction  : %.1f%%\n", x$reproduction_pct))
  invisible(x)
}


# -- internals ---------------------------------------------------------------

#' @keywords internal
#' @noRd
.subsample_rows <- function(x, max_n) {
  if (nrow(x) <= max_n) return(x)
  idx <- sample.int(nrow(x), size = max_n, replace = FALSE)
  x[idx, , drop = FALSE]
}

#' @keywords internal
#' @noRd
.scaling_from <- function(x) {
  mu <- colMeans(x)
  sd <- apply(x, 2, stats::sd)
  sd[!is.finite(sd) | sd == 0] <- 1
  list(center = mu, scale = sd)
}

#' @keywords internal
#' @noRd
.apply_scaling <- function(x, scaling) {
  sweep(sweep(x, 2, scaling$center, "-"), 2, scaling$scale, "/")
}


# -- energy distance ---------------------------------------------------------

#' @keywords internal
#' @noRd
.energy_metric <- function(ref, recon, n_self_reps = 20L) {
  ed_test <- .energy_distance(ref, recon)

  m <- nrow(recon)
  n <- nrow(ref)
  ed_self_vals <- replicate(n_self_reps, {
    idx_a <- sample.int(n, size = m, replace = TRUE)
    idx_b <- sample.int(n, size = m, replace = TRUE)
    .energy_distance(
      ref[idx_a, , drop = FALSE],
      ref[idx_b, , drop = FALSE]
    )
  })
  ed_self_mean <- mean(ed_self_vals)
  ed_self_q90  <- as.numeric(stats::quantile(ed_self_vals, 0.9, names = FALSE))

  envelope <- max(ed_self_q90, .Machine$double.eps)
  ratio    <- ed_test / envelope
  pct      <- 100 * min(1, 1 / max(1, ratio))

  list(
    reproduction_pct = pct,
    energy           = ed_test,
    self_baseline    = ed_self_mean,
    noise_envelope   = ed_self_q90,
    ratio            = ratio,
    n_self_reps      = n_self_reps
  )
}

#' @keywords internal
#' @noRd
.energy_distance <- function(x, y) {
  d_xy <- .mean_pdist_cross(x, y)
  d_xx <- .mean_pdist_within(x)
  d_yy <- .mean_pdist_within(y)
  ed2 <- 2 * d_xy - d_xx - d_yy
  sqrt(max(0, ed2))
}

#' @keywords internal
#' @noRd
.mean_pdist_within <- function(x) {
  n <- nrow(x)
  if (n < 2L) return(0)
  d <- as.numeric(stats::dist(x))
  (2 / (n^2)) * sum(d)
}

#' @keywords internal
#' @noRd
.mean_pdist_cross <- function(x, y) {
  xx <- rowSums(x^2)
  yy <- rowSums(y^2)
  d2 <- outer(xx, yy, "+") - 2 * tcrossprod(x, y)
  d2[d2 < 0] <- 0
  mean(sqrt(d2))
}

# -- C2ST --------------------------------------------------------------------

#' @keywords internal
#' @noRd
.c2st_metric <- function(ref, recon,
                         classifier = "ranger",
                         cv_folds   = 5L) {
  Z <- rbind(ref, recon)
  labels <- c(rep(0L, nrow(ref)), rep(1L, nrow(recon)))

  n <- nrow(Z)
  cv_folds <- max(2L, min(as.integer(cv_folds), n - 1L))
  folds <- sample(rep(seq_len(cv_folds), length.out = n))

  preds <- numeric(n)
  for (k in seq_len(cv_folds)) {
    test_idx  <- which(folds == k)
    train_idx <- which(folds != k)
    if (length(test_idx) == 0L || length(train_idx) == 0L) next
    preds[test_idx] <- .classifier_predict(
      classifier,
      Z[train_idx, , drop = FALSE],
      labels[train_idx],
      Z[test_idx, , drop = FALSE]
    )
  }

  auc <- .auc(labels, preds)
  pct <- 100 * max(0, 1 - 2 * abs(auc - 0.5))

  list(
    reproduction_pct = pct,
    auc              = auc,
    classifier       = classifier,
    cv_folds         = cv_folds
  )
}

#' @keywords internal
#' @noRd
.classifier_predict <- function(classifier, x_train, y_train, x_test) {
  if (classifier == "ranger") {
    # Formula interface rejects many brms/Stan parameter names (e.g. brackets).
    xx_train <- as.matrix(x_train)
    xx_test <- as.matrix(x_test)
    p <- ncol(xx_train)
    safe <- sprintf("v%04d", seq_len(p))
    colnames(xx_train) <- safe
    colnames(xx_test) <- safe
    y_fac <- factor(y_train, levels = c(0L, 1L))
    mod <- ranger::ranger(
      x           = xx_train,
      y           = y_fac,
      probability = TRUE,
      num.trees   = 200L,
      verbose     = FALSE
    )
    df_test <- as.data.frame(xx_test)
    pr <- predict(mod, df_test)$predictions
    cn <- colnames(pr)
    pr[, if ("1" %in% cn) "1" else cn[length(cn)]]
  } else {
    k <- max(3L, min(50L, floor(sqrt(nrow(x_train)))))
    .knn_proba(x_train, y_train, x_test, k = k)
  }
}

#' @keywords internal
#' @noRd
.knn_proba <- function(x_train, y_train, x_test, k) {
  xx <- rowSums(x_test^2)
  yy <- rowSums(x_train^2)
  d2 <- outer(xx, yy, "+") - 2 * tcrossprod(x_test, x_train)
  d2[d2 < 0] <- 0
  apply(d2, 1, function(d_row) {
    nn <- order(d_row)[seq_len(min(k, length(d_row)))]
    mean(y_train[nn])
  })
}

#' @keywords internal
#' @noRd
.auc <- function(labels, scores) {
  pos <- scores[labels == 1L]
  neg <- scores[labels == 0L]
  if (length(pos) == 0L || length(neg) == 0L) return(0.5)
  r <- rank(c(pos, neg), ties.method = "average")
  r_pos <- r[seq_along(pos)]
  (sum(r_pos) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}
