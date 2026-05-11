#' Compute a posterior correlation or covariance matrix
#'
#' Thin wrapper around [stats::cor()] / [stats::cov()] that:
#'
#' * accepts the same draws-like inputs as [compress_posterior()] (matrix,
#'   data.frame, `posterior::draws_*`, 3-D iter x chain x var array);
#' * keeps parameter names on rows/cols;
#' * stores a small set of attributes (`type`, `method`, `n_draws`,
#'   `n_params`) so downstream helpers like [plot_posterior_correlation()]
#'   and [partition_parameters()] can self-describe their plots/output.
#'
#' This is intended as the first step of building a hybrid compression
#' that keeps within-component covariance only for the most strongly
#' correlated parameters.
#'
#' @param draws Draws-like input (see [compress_posterior()]).
#' @param type One of `"correlation"` (default) or `"covariance"`.
#' @param variables Optional character vector of parameters to keep.
#' @param method Correlation method: one of `"pearson"` (default),
#'   `"spearman"`, `"kendall"`. Forwarded to [stats::cor()]. Ignored
#'   when `type = "covariance"`.
#' @param use Handling of missing values; forwarded to
#'   [stats::cor()] / [stats::cov()]. Default `"everything"`.
#' @param drop_zero_var Logical; if `TRUE` (default), parameters with
#'   zero posterior SD are dropped before computing correlations (else
#'   `cor()` would fill rows/cols with `NA`).
#'
#' @return A symmetric numeric matrix with row/col names equal to the
#'   parameter names, and attributes:
#'   \describe{
#'     \item{`type`}{`"correlation"` or `"covariance"`}
#'     \item{`method`}{correlation method (or `NA_character_` for cov)}
#'     \item{`n_draws`}{number of MCMC draws used}
#'     \item{`n_params`}{number of parameters in the matrix}
#'     \item{`dropped_zero_var`}{names of parameters dropped, if any}
#'   }
#'
#' @seealso [plot_posterior_correlation()], [partition_parameters()].
#' @export
#' @examples
#' set.seed(1)
#' draws <- matrix(rnorm(500 * 4), ncol = 4,
#'                 dimnames = list(NULL, c("a", "b", "c", "d")))
#' draws[, "b"] <- draws[, "a"] + 0.05 * rnorm(500)  # induce strong corr
#' cm <- posterior_correlation(draws)
#' round(cm, 2)
posterior_correlation <- function(
    draws,
    type = c("correlation", "covariance"),
    variables = NULL,
    method = c("pearson", "spearman", "kendall"),
    use = "everything",
    drop_zero_var = TRUE) {
  type   <- match.arg(type)
  method <- match.arg(method)

  mat <- .as_draws_matrix(draws, variables = variables)

  dropped <- character(0L)
  if (isTRUE(drop_zero_var)) {
    sds <- apply(mat, 2L, stats::sd)
    keep <- is.finite(sds) & sds > 0
    if (any(!keep)) {
      dropped <- colnames(mat)[!keep]
      mat <- mat[, keep, drop = FALSE]
    }
  }

  if (ncol(mat) < 2L) {
    stop(
      "Need at least 2 parameters with non-zero variance to compute a ",
      "correlation/covariance matrix; got ", ncol(mat), ".",
      call. = FALSE
    )
  }

  out <- if (type == "correlation") {
    stats::cor(mat, method = method, use = use)
  } else {
    stats::cov(mat, use = use)
  }

  attr(out, "type")             <- type
  attr(out, "method")           <- if (type == "correlation") method else NA_character_
  attr(out, "n_draws")          <- nrow(mat)
  attr(out, "n_params")         <- ncol(mat)
  attr(out, "dropped_zero_var") <- dropped
  class(out) <- c("posterior_correlation", class(out))
  out
}


#' Diagnostic plots of a posterior correlation matrix
#'
#' Produces lightweight summary plots that help pick a threshold (or
#' top-N) for [partition_parameters()]. Three views are available:
#'
#' \describe{
#'   \item{`"heatmap"`}{tile plot of the correlation matrix; ordered by
#'     hierarchical clustering on `1 - |corr|` for visual block structure;}
#'   \item{`"hist"`}{histogram of off-diagonal values (or their absolute
#'     value if `absolute = TRUE`), useful for choosing a threshold;}
#'   \item{`"max_abs"`}{sorted curve of `max_{j != i} |corr(i, j)|` per
#'     parameter, with optional reference lines at proposed thresholds.}
#' }
#'
#' Requires `ggplot2` (Suggests). Returns `ggplot` objects rather than
#' drawing, so callers can compose / save them.
#'
#' @param x A matrix from [posterior_correlation()] (or any symmetric
#'   numeric matrix; rows/cols are treated as parameters).
#' @param which One or more of `"heatmap"`, `"hist"`, `"max_abs"`. The
#'   default returns all three.
#' @param absolute Logical; if `TRUE` (default for `"hist"` and
#'   `"max_abs"`) summarise on the absolute scale, since both signs of
#'   correlation are equally informative for partitioning.
#' @param threshold Optional numeric reference threshold(s) drawn on the
#'   `"hist"` and `"max_abs"` panels (e.g. `0.3`).
#' @param max_params_heatmap Integer; if the matrix has more parameters
#'   than this, the heatmap is restricted to the top-`max_params_heatmap`
#'   parameters by `max_{j != i} |corr(i, j)|`. Avoids unreadable
#'   1000x1000 tiles. Default `200`.
#' @param dendrogram One of `"none"` (default), `"both"`, `"row"`, or
#'   `"column"`. Adds dendrograms next to the heatmap (return value becomes
#'   a `patchwork` object when not `"none"`). Requires the `ggdendro` and
#'   `patchwork` packages (Suggests); silently falls back to a plain heatmap
#'   if either is missing.
#' @param partition Optional [partition_parameters_clusters()] result
#'   (**`simple_output = TRUE`**: plain list of cluster vectors; or
#'   **`simple_output = FALSE`**: a `poco_partition_clusters` object), or any
#'   list with the same fields as the full object (`cluster_id`, `blocks`,
#'   `remainder`, ...). For **simple** lists (no `cluster_id`), cluster
#'   outlines are drawn by mapping block membership to the **visible**
#'   heatmap columns (after any `max_params_heatmap` subset and dendrogram
#'   reorder). Parameters not in any block are treated as remainder for the
#'   annotation subtitle.
#' @param show_cluster_labels Logical or `NULL` (default). If `NULL`, on-plot
#'   cluster **text** is omitted when there are many clusters or the heatmap
#'   is dense (few parameters per block); outlines are always drawn. Use
#'   `TRUE` / `FALSE` to override.
#' @param show_cluster_legend Logical or `NULL` (default). If `NULL`, the
#'   cluster **colour legend** follows similar heuristics (slightly more
#'   permissive than labels). Use `TRUE` / `FALSE` to override.
#' @param ... Currently unused.
#'
#' @return A named list of `ggplot` objects (one entry per `which`).
#'   When called for a single panel, returns the single `ggplot`.
#'
#' @seealso [posterior_correlation()], [partition_parameters()].
#' @export
plot_posterior_correlation <- function(
    x,
    which = c("heatmap", "hist", "max_abs"),
    absolute = TRUE,
    threshold = NULL,
    max_params_heatmap = 200L,
    dendrogram = c("none", "both", "row", "column"),
    partition = NULL,
    show_cluster_labels = NULL,
    show_cluster_legend = NULL,
    ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for plot_posterior_correlation().\n",
      "  install.packages('ggplot2')",
      call. = FALSE
    )
  }

  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x)) {
    stop("`x` must be a square numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(x))) {
    rownames(x) <- colnames(x) <- paste0("V", seq_len(nrow(x)))
  } else if (is.null(colnames(x))) {
    colnames(x) <- rownames(x)
  }

  which      <- match.arg(which, several.ok = TRUE)
  dendrogram <- match.arg(dendrogram)

  type <- attr(x, "type") %||% "correlation"

  abs_offdiag <- .abs_offdiag_per_param(x)

  plots <- list()

  if ("heatmap" %in% which) {
    plots$heatmap <- .plot_corr_heatmap(
      x,
      type                 = type,
      max_params_heatmap   = max_params_heatmap,
      abs_offdiag          = abs_offdiag,
      dendrogram           = dendrogram,
      partition            = partition,
      show_cluster_labels  = show_cluster_labels,
      show_cluster_legend  = show_cluster_legend
    )
  }

  if ("hist" %in% which) {
    plots$hist <- .plot_corr_hist(
      x,
      type      = type,
      absolute  = absolute,
      threshold = threshold
    )
  }

  if ("max_abs" %in% which) {
    plots$max_abs <- .plot_corr_max_abs(
      abs_offdiag,
      type      = type,
      threshold = threshold
    )
  }

  if (length(plots) == 1L) plots[[1L]] else plots
}


#' Partition parameters into a "correlated" block A and a weaker block B
#'
#' Given a correlation matrix (e.g. from [posterior_correlation()]),
#' classify each parameter as belonging to **block A** (strongly
#' correlated with at least one other parameter) or **block B**
#' (everything else) using one of several explicit rules.
#'
#' The intended use is in a hybrid compression: fit a richer
#' (non-diagonal) mixture on **block A**, and keep a cheap
#' diagonal/independent approximation on **block B**.
#'
#' @param x A matrix from [posterior_correlation()] (or any square
#'   correlation matrix).
#' @param rule One of:
#'   \describe{
#'     \item{`"threshold"`}{block A = parameters whose
#'       `max_{j != i} |corr(i, j)| >= threshold`. Requires `threshold`.}
#'     \item{`"top_n"`}{block A = the `n` parameters with the largest
#'       max-off-diagonal `|corr|`. Requires `n`.}
#'     \item{`"top_prop"`}{block A = the top `prop` (e.g. `0.1` for the
#'       top 10%) of parameters by max-off-diagonal `|corr|`. Requires
#'       `prop` in `(0, 1]`.}
#'     \item{`"min_degree"`}{block A = parameters with at least `min_degree`
#'       neighbours whose `|corr| >= threshold` (graph-degree rule).
#'       Requires both `threshold` and `min_degree`.}
#'   }
#' @param threshold Numeric; required for `"threshold"` and
#'   `"min_degree"`.
#' @param n Integer; required for `"top_n"`.
#' @param prop Numeric in `(0, 1]`; required for `"top_prop"`.
#' @param min_degree Integer; required for `"min_degree"`.
#'
#' @return A list of class `"poco_partition"` with components:
#'   \describe{
#'     \item{`block_a`}{character vector of parameter names assigned to A;}
#'     \item{`block_b`}{character vector of parameter names assigned to B;}
#'     \item{`abs_max_offdiag`}{named numeric vector of per-parameter
#'       max-off-diagonal `|corr|` (same names as `x`'s rows);}
#'     \item{`rule`}{the rule used;}
#'     \item{`params`}{a named list of the rule's parameters
#'       (`threshold`, `n`, `prop`, `min_degree`).}
#'   }
#'
#' @seealso [posterior_correlation()], [plot_posterior_correlation()].
#' @export
#' @examples
#' set.seed(2)
#' draws <- matrix(rnorm(500 * 5), ncol = 5,
#'                 dimnames = list(NULL, letters[1:5]))
#' draws[, "b"] <- draws[, "a"] + 0.05 * rnorm(500)
#' draws[, "c"] <- draws[, "a"] - 0.05 * rnorm(500)
#' cm <- posterior_correlation(draws)
#' partition_parameters(cm, rule = "threshold", threshold = 0.3)
partition_parameters <- function(
    x,
    rule = c("threshold", "top_n", "top_prop", "min_degree"),
    threshold = NULL,
    n = NULL,
    prop = NULL,
    min_degree = NULL) {
  rule <- match.arg(rule)

  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x)) {
    stop("`x` must be a square numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(x))) {
    rownames(x) <- colnames(x) <- paste0("V", seq_len(nrow(x)))
  }

  abs_max_offdiag <- .abs_offdiag_per_param(x)
  param_names     <- names(abs_max_offdiag)

  block_a <- switch(
    rule,
    threshold = {
      if (is.null(threshold)) {
        stop("`threshold` is required for rule = 'threshold'.", call. = FALSE)
      }
      param_names[abs_max_offdiag >= threshold]
    },
    top_n = {
      if (is.null(n)) {
        stop("`n` is required for rule = 'top_n'.", call. = FALSE)
      }
      n <- min(as.integer(n), length(param_names))
      ord <- order(abs_max_offdiag, decreasing = TRUE)
      param_names[ord[seq_len(n)]]
    },
    top_prop = {
      if (is.null(prop) || prop <= 0 || prop > 1) {
        stop("`prop` must be in (0, 1].", call. = FALSE)
      }
      n_pick <- max(1L, ceiling(prop * length(param_names)))
      ord <- order(abs_max_offdiag, decreasing = TRUE)
      param_names[ord[seq_len(n_pick)]]
    },
    min_degree = {
      if (is.null(threshold) || is.null(min_degree)) {
        stop(
          "`threshold` and `min_degree` are required for rule = 'min_degree'.",
          call. = FALSE
        )
      }
      adj <- abs(x) >= threshold
      diag(adj) <- FALSE
      deg <- rowSums(adj)
      param_names[deg >= min_degree]
    }
  )

  block_b <- setdiff(param_names, block_a)

  out <- list(
    block_a         = block_a,
    block_b         = block_b,
    abs_max_offdiag = abs_max_offdiag,
    rule            = rule,
    params          = list(
      threshold  = threshold,
      n          = n,
      prop       = prop,
      min_degree = min_degree
    )
  )
  class(out) <- c("poco_partition", "list")
  out
}


#' Partition parameters into multiple correlated clusters + a remainder
#'
#' Hierarchical-clustering counterpart to [partition_parameters()]. Builds
#' clusters from the same `1 - |corr|` distance used by the heatmap, cuts
#' the dendrogram by either a correlation **`threshold`** (cluster members
#' have average linkage `|corr| >= threshold`) **or** by a target number
#' of clusters **`k`**, and then keeps only clusters whose size is at
#' least **`min_size`**. Parameters not in any kept cluster are collected
#' in **`remainder`** (weakly correlated / small groups / singletons).
#'
#' Intended use: feed each cluster (block A1, A2, ...) into a richer
#' (non-diagonal) `mclust` family, while keeping **`remainder`** on a cheap
#' diagonal/independent approximation.
#'
#' @param x A correlation matrix (e.g. from [posterior_correlation()]).
#' @param threshold Numeric in `(0, 1]`; cut the dendrogram at height
#'   `1 - threshold`. Mutually exclusive with `k`.
#' @param k Integer; cut the dendrogram into exactly `k` clusters.
#'   Mutually exclusive with `threshold`.
#' @param min_size Minimum cluster size for a group to be kept as its own
#'   block. Either a **positive integer** (count), or a **numeric proportion
#'   strictly between 0 and 1**, in which case it is converted with
#'   `ceiling(min_size * n_parameters)` and an informative message is printed.
#'   For count-style inputs, negative values are converted to their absolute
#'   value and non-integer values are rounded up with `ceiling()`, both with
#'   warnings. The effective minimum is then clamped to at least `2`; if
#'   conversion or direct input yields `< 2`, it is reset to `2` with a warning.
#'   Default `10` (integer count).
#' @param linkage Linkage method for [stats::hclust()]. Default
#'   `"average"` (matches what the heatmap uses to reorder).
#' @param force_remainder Optional; parameters forced into **`remainder`**
#'   (never assigned to **`blocks`**), regardless of clustering. `NULL`
#'   (default) keeps the partition as cut from the tree. Otherwise supply
#'   a tidy-select expression (interpreted like [dplyr::select()] column
#'   semantics), e.g. `tidyselect::starts_with("cor_")`,
#'   `tidyselect::contains("gamma")`, or a `c()` / character vector of names.
#'   Names must appear in `rownames(x)`.
#' @param simple_output Logical. If `TRUE` (default), return a **plain**
#'   [base::list()] of character vectors only: one element per kept cluster
#'   (names `cluster_1`, `cluster_2`, ...), same shape as a manual
#'   **`partition`** for [compress_posterior()] / [compress_brmsfit()]. All
#'   **cluster-related metadata is omitted** (`cluster_id`, `hclust`,
#'   `params`, and there is no `remainder` component). Parameters **not**
#'   appearing in any list element are still the **remainder** for
#'   blockwise compression: downstream code treats them as unclustered. If
#'   `FALSE`, return the full **`poco_partition_clusters`** object below
#'   (for [plot_posterior_correlation()] cluster outlines, printing, and
#'   inspection of `remainder` explicitly).
#'
#' @return If **`simple_output`** is `TRUE`: a plain named `list` of
#'   character vectors (`cluster_1`, `cluster_2`, ...), or **`list()`** when
#'   no cluster passes **`min_size`** (all parameters are remainder for
#'   [compress_posterior()] / [compress_brmsfit()]). If `FALSE`, a list
#'   of class `"poco_partition_clusters"` with components:
#'   \describe{
#'     \item{`blocks`}{named list of character vectors (`cluster_1`,
#'       `cluster_2`, ...), sorted by descending size;}
#'     \item{`remainder`}{character vector of parameters not in any kept cluster;}
#'     \item{`cluster_id`}{named integer vector mapping every parameter
#'       to its kept-cluster id (1..K) or `NA` if it landed in `remainder`;}
#'     \item{`hclust`}{the [stats::hclust()] object used for cutting and
#'       for reordering (re-used by [plot_posterior_correlation()]);}
#'     \item{`params`}{the rule's parameters (`threshold`, `k`, `min_size`
#'       effective count, optional `min_size_proportion` if a proportion was
#'       given, `linkage`, optional `force_remainder` if any names were forced).}
#'   }
#'
#' @seealso [posterior_correlation()], [partition_parameters()],
#'   [plot_posterior_correlation()].
#' @importFrom tidyselect eval_select
#' @export
#' @examples
#' set.seed(3)
#' n <- 600L
#' a <- matrix(rnorm(n), n, 3)
#' a[, 2] <- a[, 1] + 0.05 * rnorm(n)
#' a[, 3] <- a[, 1] - 0.05 * rnorm(n)
#' b <- matrix(rnorm(n * 4L), n, 4L)
#' draws <- cbind(a, b)
#' colnames(draws) <- c("a1", "a2", "a3", "b1", "b2", "b3", "b4")
#' cm <- posterior_correlation(draws)
#' # Default: plain list of clusters (pass-through to compress_* `partition`)
#' partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L)
#' # Full metadata (heatmap outlines, `$remainder`, `$hclust`, ...)
#' partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L,
#'   simple_output = FALSE)
#' # Proportion in (0, 1): minimum size = ceiling(0.25 * ncol(cm))
#' partition_parameters_clusters(cm, threshold = 0.5, min_size = 0.25)
#' # Keep cor_* / residual correlations in the cheap remainder block:
#' partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L,
#'   force_remainder = tidyselect::starts_with("cor_"))
#' partition_parameters_clusters(cm, threshold = 0.5, min_size = 2L,
#'   force_remainder = c("a1", "a2"))
partition_parameters_clusters <- function(
    x,
    threshold = NULL,
    k         = NULL,
    min_size  = 10L,
    linkage   = "average",
    force_remainder = NULL,
    simple_output = TRUE) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x)) {
    stop("`x` must be a square numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(x))) {
    rownames(x) <- colnames(x) <- paste0("V", seq_len(nrow(x)))
  }
  if (is.null(threshold) && is.null(k)) {
    stop("Provide either `threshold` or `k`.", call. = FALSE)
  }
  if (!is.null(threshold) && !is.null(k)) {
    stop("Provide only one of `threshold` or `k`.", call. = FALSE)
  }

  n_params <- nrow(x)
  min_size_in <- min_size
  min_size_prop <- NULL
  if (length(min_size_in) != 1L) {
    stop("`min_size` must be length 1.", call. = FALSE)
  }
  if (is.numeric(min_size_in) && !is.na(min_size_in) &&
      min_size_in > 0 && min_size_in < 1) {
    min_size_prop <- as.numeric(min_size_in)
    min_size <- as.integer(ceiling(min_size_prop * n_params))
    .message(
      "`min_size` = ", min_size_prop,
      " is a proportion in (0, 1): using minimum cluster size ",
      min_size, " = ceiling(", min_size_prop, " * ", n_params, ")."
    )
  } else {
    min_size_num <- suppressWarnings(as.numeric(min_size_in))
    if (length(min_size_num) != 1L || is.na(min_size_num) ||
        !is.finite(min_size_num)) {
      stop(
        "`min_size` must be numeric: either a count-like value or a ",
        "proportion in (0, 1).",
        call. = FALSE
      )
    }
    if (min_size_num < 0) {
      warning(
        "`min_size` is negative (", min_size_num,
        "); using absolute value ", abs(min_size_num), ".",
        call. = FALSE
      )
      min_size_num <- abs(min_size_num)
    }
    min_size_num_up <- ceiling(min_size_num)
    if (!isTRUE(all.equal(min_size_num, min_size_num_up))) {
      warning(
        "`min_size` should be an integer count; using ceiling(",
        min_size_num, ") = ", min_size_num_up, ".",
        call. = FALSE
      )
    }
    min_size <- as.integer(min_size_num_up)
  }
  if (min_size < 2L) {
    warning(
      "`min_size` resolved to ", min_size,
      ", but the minimum allowed effective cluster size is 2; using 2.",
      call. = FALSE
    )
    min_size <- 2L
  }

  d  <- stats::as.dist(1 - abs(x))
  hc <- stats::hclust(d, method = linkage)

  raw_id <- if (!is.null(k)) {
    stats::cutree(hc, k = as.integer(k))
  } else {
    stats::cutree(hc, h = 1 - threshold)
  }
  names(raw_id) <- rownames(x)

  sizes    <- table(raw_id)
  keep_ids <- as.integer(names(sizes)[sizes >= min_size])
  if (length(keep_ids) == 0L) {
    blocks <- list()
    remainder <- rownames(x)
    new_id <- rep(NA_integer_, length(raw_id))
    names(new_id) <- names(raw_id)
  } else {
    keep_ids <- keep_ids[
      order(as.integer(sizes[as.character(keep_ids)]), decreasing = TRUE)
    ]
    blocks <- lapply(keep_ids, function(id) {
      names(raw_id)[raw_id == id]
    })
    names(blocks) <- paste0("cluster_", seq_along(blocks))

    remap <- stats::setNames(
      seq_along(keep_ids),
      as.character(keep_ids)
    )
    new_id <- rep(NA_integer_, length(raw_id))
    matched <- as.character(raw_id) %in% names(remap)
    new_id[matched] <- as.integer(remap[as.character(raw_id[matched])])
    names(new_id) <- names(raw_id)

    remainder <- names(raw_id)[is.na(new_id)]
  }

  force_remainder_applied <- NULL
  force_quo <- rlang::enquo(force_remainder)
  if (!rlang::quo_is_null(force_quo)) {
    sel_df <- as.data.frame(
      matrix(
        ncol = n_params,
        nrow = 0L,
        dimnames = list(NULL, rownames(x))
      ),
      check.names = FALSE
    )
    forced_pos <- tidyselect::eval_select(
      force_quo,
      data = sel_df,
      allow_rename = FALSE
    )
    forced_names <- names(forced_pos)
    if (length(forced_names) > 0L) {
      blocks <- lapply(blocks, function(b) setdiff(b, forced_names))
      blocks <- blocks[lengths(blocks) > 0L]
      if (length(blocks)) {
        names(blocks) <- paste0("cluster_", seq_along(blocks))
      }
      new_id <- rep(NA_integer_, n_params)
      names(new_id) <- rownames(x)
      if (length(blocks)) {
        for (i in seq_along(blocks)) {
          new_id[blocks[[i]]] <- as.integer(i)
        }
      }
      remainder <- names(new_id)[is.na(new_id)]
      force_remainder_applied <- forced_names
    }
  }

  params <- list(
    threshold = threshold,
    k         = k,
    min_size  = min_size,
    linkage   = linkage
  )
  if (!is.null(min_size_prop)) {
    params$min_size_proportion <- min_size_prop
  }
  if (!is.null(force_remainder_applied)) {
    params$force_remainder <- force_remainder_applied
  }

  if (isTRUE(simple_output)) {
    out <- blocks
    class(out) <- "list"
    return(out)
  }

  out <- list(
    blocks     = blocks,
    remainder  = remainder,
    cluster_id = new_id,
    hclust     = hc,
    params     = params
  )
  class(out) <- c("poco_partition_clusters", "list")
  out
}


#' @export
print.poco_partition_clusters <- function(x, ...) {
  n_total <- length(x$cluster_id)
  rem <- .partition_remainder(x)
  rp <- x$params[!vapply(x$params, is.null, logical(1L))]
  rp_str <- paste(names(rp), unlist(rp), sep = " = ", collapse = ", ")
  cat("<poco_partition_clusters>\n")
  cat("  parameters : ", rp_str, "\n", sep = "")
  cat(
    "  kept       : ", length(x$blocks),
    " cluster(s); ", n_total - length(rem),
    " of ", n_total, " params",
    sprintf(" (%.1f%%)\n", 100 * (n_total - length(rem)) / max(n_total, 1L)),
    sep = ""
  )
  cat(
    "  remainder  : ", length(rem), " / ", n_total,
    sprintf(" (%.1f%%)\n", 100 * length(rem) / max(n_total, 1L)),
    sep = ""
  )
  if (length(x$blocks)) {
    sizes <- vapply(x$blocks, length, integer(1L))
    show <- min(10L, length(sizes))
    cat("  block sizes:\n")
    for (i in seq_len(show)) {
      cat(sprintf("    %s : %d\n", names(sizes)[i], sizes[[i]]))
    }
    if (length(sizes) > show) {
      cat("    ... (", length(sizes) - show, " more)\n", sep = "")
    }
  }
  invisible(x)
}


#' @export
print.poco_partition <- function(x, ...) {
  n_total <- length(x$block_a) + length(x$block_b)
  cat("<poco_partition>\n")
  cat("  rule       : ", x$rule, "\n", sep = "")
  rp <- x$params[!vapply(x$params, is.null, logical(1L))]
  if (length(rp)) {
    rp_str <- paste(names(rp), unlist(rp), sep = " = ", collapse = ", ")
    cat("  parameters : ", rp_str, "\n", sep = "")
  }
  cat(
    "  block A    : ", length(x$block_a), " / ", n_total,
    sprintf(" (%.1f%%)\n", 100 * length(x$block_a) / max(n_total, 1L)),
    sep = ""
  )
  cat(
    "  block B    : ", length(x$block_b), " / ", n_total,
    sprintf(" (%.1f%%)\n", 100 * length(x$block_b) / max(n_total, 1L)),
    sep = ""
  )
  if (length(x$block_a)) {
    head_n <- min(5L, length(x$block_a))
    cat(
      "  A head     : ",
      paste(utils::head(x$block_a, head_n), collapse = ", "),
      if (length(x$block_a) > head_n) ", ..." else "",
      "\n",
      sep = ""
    )
  }
  invisible(x)
}


# ----------------------------- internals ------------------------------- #

#' Per-parameter max off-diagonal |correlation|
#' @keywords internal
#' @noRd
.abs_offdiag_per_param <- function(x) {
  m <- abs(x)
  diag(m) <- NA_real_
  out <- suppressWarnings(apply(m, 2L, max, na.rm = TRUE))
  out[!is.finite(out)] <- 0
  if (is.null(names(out))) names(out) <- colnames(x)
  out
}


#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Parameters not assigned to any kept cluster (`remainder`).
#' @keywords internal
#' @noRd
.partition_remainder <- function(partition) {
  if (is.null(partition)) {
    return(character())
  }
  if (!is.null(partition$remainder)) {
    return(partition$remainder)
  }
  character()
}


#' Named vector of outline/label colours for cluster block names.
#' @keywords internal
#' @noRd
.cluster_outline_palette <- function(block_names) {
  k <- length(block_names)
  if (k < 1L) {
    return(stats::setNames(character(), character()))
  }
  cols <- grDevices::palette.colors(
    n = k,
    palette = "Classic Tableau",
    recycle = TRUE
  )
  stats::setNames(as.character(cols), block_names)
}


#' Heuristics for cluster text labels vs colour legend on correlation heatmaps.
#'
#' @keywords internal
#' @noRd
.auto_cluster_heatmap_annotation <- function(n_blocks, n_params) {
  if (n_blocks < 1L) {
    return(list(
      show_cluster_labels = FALSE,
      show_cluster_legend = FALSE
    ))
  }
  avg_block <- n_params / n_blocks
  show_labels <-
    n_blocks <= 6L && n_params <= 350L && avg_block >= 5
  show_legend <- n_blocks <= 12L && n_params <= 400L
  if (n_blocks > 10L) {
    show_legend <- FALSE
  }
  list(
    show_cluster_labels = show_labels,
    show_cluster_legend = show_legend
  )
}


#' One row per distinct **displayed** cluster name: keeps the contiguous run
#' with largest tile area. Using the label string avoids duplicate on-plot text
#' when `cluster_id` mixes types (e.g. integer vs character) or other edge
#' cases that still map to the same `cluster_k` name.
#'
#' @keywords internal
#' @noRd
.pick_cluster_label_row_per_id <- function(rect_df) {
  if (!nrow(rect_df)) {
    return(rect_df)
  }
  area <- (rect_df$xmax - rect_df$xmin) * (rect_df$ymax - rect_df$ymin)
  key <- if ("cluster_lab" %in% names(rect_df)) {
    as.character(rect_df$cluster_lab)
  } else {
    as.character(rect_df$cluster_id)
  }
  uid <- unique(key)
  pick <- vapply(uid, function(nm) {
    idx <- which(key == nm)
    idx[which.max(area[idx])]
  }, integer(1L))
  rect_df[sort(unname(pick)), , drop = FALSE]
}


#' Build `cluster_id` / `blocks` / `remainder` for heatmap annotation when
#' `partition` is not a `poco_partition_clusters` object and has no
#' `cluster_id` (e.g. [partition_parameters_clusters()] with default
#' **`simple_output = TRUE`**, or any plain `list()` of character vectors
#' suitable for `compress_*` `partition`). Visible columns only. Lists whose
#' elements are not all non-empty character vectors are returned unchanged.
#'
#' @keywords internal
#' @noRd
.enrich_partition_for_heatmap_plot <- function(partition, visible, linkage) {
  if (is.null(partition)) {
    return(NULL)
  }
  if (inherits(partition, "poco_partition_clusters")) {
    return(partition)
  }
  if (!is.null(partition[["cluster_id"]])) {
    return(partition)
  }
  if (!is.list(partition) || !length(partition)) {
    return(partition)
  }
  bl <- list()
  for (i in seq_along(partition)) {
    el <- partition[[i]]
    if (!is.character(el) || !length(el)) {
      return(partition)
    }
    hit <- intersect(el, visible)
    hit <- unique(hit[nzchar(hit)])
    if (length(hit)) {
      bl[[length(bl) + 1L]] <- hit
    }
  }
  if (!length(bl)) {
    return(partition)
  }
  names(bl) <- paste0("cluster_", seq_along(bl))
  rem <- setdiff(visible, unlist(bl, use.names = FALSE))
  cid <- rep(NA_integer_, length(visible))
  names(cid) <- visible
  for (k in seq_along(bl)) {
    cid[bl[[k]]] <- k
  }
  list(
    blocks     = bl,
    remainder  = rem,
    cluster_id = cid,
    params     = list(linkage = linkage)
  )
}


#' @keywords internal
#' @noRd
.plot_corr_heatmap <- function(x, type, max_params_heatmap, abs_offdiag,
                               dendrogram = "both",
                               partition  = NULL,
                               show_cluster_labels = NULL,
                               show_cluster_legend = NULL) {
  full_n <- ncol(x)
  if (ncol(x) > max_params_heatmap) {
    keep <- names(sort(abs_offdiag, decreasing = TRUE))[seq_len(max_params_heatmap)]
    x <- x[keep, keep, drop = FALSE]
    subtitle <- paste0(
      "Showing top ", max_params_heatmap,
      " of ", full_n,
      " parameters by max |", type, "|"
    )
  } else {
    subtitle <- paste0(ncol(x), " parameters")
  }

  linkage <- if (!is.null(partition) &&
                 !is.null(partition$params$linkage)) {
    partition$params$linkage
  } else {
    "average"
  }

  hc <- tryCatch(
    stats::hclust(stats::as.dist(1 - abs(x)), method = linkage),
    error = function(e) NULL
  )
  ord <- if (!is.null(hc)) hc$order else seq_len(ncol(x))
  x <- x[ord, ord, drop = FALSE]

  long <- data.frame(
    row   = factor(rep(rownames(x), times = ncol(x)), levels = rownames(x)),
    col   = factor(rep(colnames(x), each  = nrow(x)), levels = colnames(x)),
    value = as.vector(x)
  )

  fill_lim <- if (type == "correlation") c(-1, 1) else NULL

  p <- ggplot2::ggplot(long, ggplot2::aes(.data$col, .data$row, fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low      = "blue",
      mid      = "white",
      high     = "red",
      midpoint = 0,
      limits   = fill_lim,
      name     = type
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::labs(
      title    = paste0("Posterior ", type, " heatmap"),
      subtitle = subtitle,
      x        = NULL,
      y        = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank(),
      legend.position = "right"
    )

  partition_plot <- .enrich_partition_for_heatmap_plot(
    partition,
    colnames(x),
    linkage
  )

  if (!is.null(partition_plot) && !is.null(partition_plot$cluster_id)) {
    visible <- colnames(x)
    cid <- partition_plot$cluster_id[visible]
    runs <- .contiguous_runs(cid)
    n_blocks <- length(partition_plot$blocks)
    auto_ann <- .auto_cluster_heatmap_annotation(n_blocks, length(visible))
    if (is.null(show_cluster_labels)) {
      show_cluster_labels <- auto_ann$show_cluster_labels
    }
    if (is.null(show_cluster_legend)) {
      show_cluster_legend <- auto_ann$show_cluster_legend
    }
    if (nrow(runs) > 0L) {
      block_names <- names(partition_plot$blocks)
      pal <- .cluster_outline_palette(block_names)
      lab <- ifelse(
        runs$cluster >= 1L & runs$cluster <= length(block_names),
        block_names[runs$cluster],
        paste0("C", runs$cluster)
      )
      rect_df <- data.frame(
        xmin        = runs$start - 0.5,
        xmax        = runs$end   + 0.5,
        ymin        = runs$start - 0.5,
        ymax        = runs$end   + 0.5,
        cluster_id  = as.integer(runs$cluster),
        cluster_lab = lab,
        stringsAsFactors = FALSE
      )
      rect_df$outline_col <- pal[rect_df$cluster_lab]
      rect_df$outline_col[
        is.na(rect_df$outline_col) | !nzchar(rect_df$outline_col)
      ] <- "#333333"
      rect_df$x_text <- (rect_df$xmin + rect_df$xmax) / 2
      rect_df$y_text <- (rect_df$ymin + rect_df$ymax) / 2

      label_df <- .pick_cluster_label_row_per_id(rect_df)

      rect_df$cluster_lab <- factor(
        rect_df$cluster_lab,
        levels = names(pal)
      )
      label_df$cluster_lab <- factor(
        as.character(label_df$cluster_lab),
        levels = names(pal)
      )

      colour_scale <- ggplot2::scale_colour_manual(
        name     = "Cluster",
        values   = pal,
        drop     = FALSE,
        na.value = "#333333",
        guide    = if (isTRUE(show_cluster_legend)) {
          ggplot2::guide_legend(order = 1)
        } else {
          "none"
        }
      )

      p <- p +
        ggplot2::geom_rect(
          data = rect_df,
          ggplot2::aes(
            xmin = .data$xmin, xmax = .data$xmax,
            ymin = .data$ymin, ymax = .data$ymax,
            colour = .data$cluster_lab
          ),
          inherit.aes = FALSE,
          fill        = NA,
          linewidth   = 1.35,
          show.legend = isTRUE(show_cluster_legend)
        ) +
        colour_scale

      if (isTRUE(show_cluster_labels)) {
        p <- p +
          ggplot2::geom_label(
            data = label_df,
            ggplot2::aes(
              x      = .data$x_text,
              y      = .data$y_text,
              label  = .data$cluster_lab,
              colour = .data$cluster_lab
            ),
            inherit.aes = FALSE,
            fill         = ggplot2::alpha("white", 0.92),
            label.size   = 0.45,
            size         = 3.2,
            fontface     = "bold",
            show.legend  = FALSE
          )
      }
    }
    n_rem <- length(.partition_remainder(partition_plot))
    if (length(stats::na.omit(as.integer(cid))) > 0L) {
      p <- p +
        ggplot2::labs(
          subtitle = paste0(
            subtitle,
            "  |  ", n_blocks, " cluster(s) marked; ",
            n_rem, " in remainder"
          )
        )
    }
  }

  if (dendrogram != "none" && !is.null(hc)) {
    seg <- .dendro_segments(hc)
    if (!is.null(seg)) {
      p <- .compose_heatmap_with_dendro(p, seg, dendrogram)
    }
  }

  p
}


#' Contiguous runs of a (possibly NA-bearing) integer cluster id vector.
#'
#' Returns a data.frame with columns `cluster`, `start`, `end` for runs
#' where the cluster id is not NA.
#' @keywords internal
#' @noRd
.contiguous_runs <- function(cid) {
  cid <- as.integer(cid)
  n <- length(cid)
  if (n == 0L) {
    return(data.frame(cluster = integer(), start = integer(), end = integer()))
  }
  marker <- ifelse(is.na(cid), -1L, cid)
  rl <- rle(marker)
  ends   <- cumsum(rl$lengths)
  starts <- ends - rl$lengths + 1L
  df <- data.frame(cluster = rl$values, start = starts, end = ends)
  df[df$cluster > 0L, , drop = FALSE]
}


#' Get ggdendro segments for a hclust, or NULL if ggdendro missing.
#' @keywords internal
#' @noRd
.dendro_segments <- function(hc) {
  if (!requireNamespace("ggdendro", quietly = TRUE)) {
    return(NULL)
  }
  ggdendro::dendro_data(hc, type = "rectangle")$segments
}


#' Compose a heatmap ggplot with dendrogram(s) on top/left via patchwork.
#' Returns the original plot unchanged if patchwork is missing.
#' @keywords internal
#' @noRd
.compose_heatmap_with_dendro <- function(heatmap_plot, segments, where) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    return(heatmap_plot)
  }

  blank <- ggplot2::ggplot() + ggplot2::theme_void()

  top <- ggplot2::ggplot(
    segments,
    ggplot2::aes(
      x = .data$x, y = .data$y,
      xend = .data$xend, yend = .data$yend
    )
  ) +
    ggplot2::geom_segment(linewidth = 0.3) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

  left <- ggplot2::ggplot(
    segments,
    ggplot2::aes(
      x = -.data$y, y = .data$x,
      xend = -.data$yend, yend = .data$xend
    )
  ) +
    ggplot2::geom_segment(linewidth = 0.3) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0))) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

  switch(
    where,
    column = patchwork::wrap_plots(
      top, heatmap_plot,
      ncol = 1L, heights = c(0.18, 1)
    ),
    row = patchwork::wrap_plots(
      left, heatmap_plot,
      nrow = 1L, widths = c(0.18, 1)
    ),
    both = patchwork::wrap_plots(
      blank, top, left, heatmap_plot,
      ncol = 2L, nrow = 2L,
      widths  = c(0.18, 1),
      heights = c(0.18, 1)
    ),
    heatmap_plot
  )
}


#' @keywords internal
#' @noRd
.plot_corr_hist <- function(x, type, absolute, threshold) {
  vals <- x[upper.tri(x)]
  if (isTRUE(absolute)) vals <- abs(vals)
  df <- data.frame(value = vals)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$value)) +
    ggplot2::geom_histogram(bins = 60, fill = "steelblue",
                            color = "white", alpha = 0.85) +
    ggplot2::labs(
      title    = paste0(
        "Off-diagonal ", type,
        if (isTRUE(absolute)) " (|.|)" else ""
      ),
      subtitle = paste(length(vals), "unique pairs"),
      x        = if (isTRUE(absolute)) paste0("|", type, "|") else type,
      y        = "count"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(threshold)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = threshold,
        linetype   = "dashed",
        color      = "red"
      )
  }
  p
}


#' @keywords internal
#' @noRd
.plot_corr_max_abs <- function(abs_offdiag, type, threshold) {
  ord <- order(abs_offdiag, decreasing = TRUE)
  df <- data.frame(
    rank  = seq_along(abs_offdiag),
    value = abs_offdiag[ord],
    name  = factor(names(abs_offdiag)[ord], levels = names(abs_offdiag)[ord])
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$rank, .data$value)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::labs(
      title    = paste0("Max off-diagonal |", type, "| per parameter"),
      subtitle = "Sorted descending; pick a threshold or top-N for partition_parameters()",
      x        = "Parameter rank",
      y        = paste0("max_{j != i} |", type, "(i, j)|")
    ) +
    ggplot2::theme_minimal()

  if (!is.null(threshold)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = threshold,
        linetype   = "dashed",
        color      = "red"
      )
  }
  p
}
