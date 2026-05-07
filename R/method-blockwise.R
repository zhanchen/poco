#' Blockwise (per-cluster) posterior compression
#'
#' Internal worker behind the `partition` argument of [compress_posterior()]
#' and [compress_brmsfit()].
#'
#' Each cluster of `partition$blocks` is compressed independently with a
#' typically richer covariance family (`cluster_model_name`, default
#' auto-selected by mclust/BIC), and the leftover `partition$remainder`
#' parameters are compressed together with a cheaper diagonal family
#' (`remainder_model_name`, default `"VVI"`).
#'
#' The result is wrapped in a `posterior_compressed_blockwise` object
#' that [sample_posterior()] and [density_posterior()] dispatch on:
#' samples from each block are drawn independently and then re-aligned
#' with the original parameter order, so [reconstruct_brmsfit()] keeps
#' working unchanged.
#'
#' This implies a **between-block independence** approximation: posterior
#' correlations *within* each cluster are preserved by the richer
#' mixture, but correlations *between* blocks (cluster-to-cluster, and
#' cluster-to-remainder) are dropped.
#'
#' @keywords internal
#' @noRd
.compress_blockwise <- function(
    draws_mat,
    partition,
    method               = "mclust",
    n_components         = 3L,
    cluster_model_name   = NULL,
    remainder_model_name = NULL,
    verbose              = FALSE,
    ...) {
  param_names <- colnames(draws_mat)
  if (is.null(param_names)) {
    stop("Blockwise compression requires `draws_mat` to have column names.",
         call. = FALSE)
  }

  if (!inherits(partition, "poco_partition_clusters")) {
    stop(
      "`partition` must be a `poco_partition_clusters` object ",
      "(see `partition_parameters_clusters()`).",
      call. = FALSE
    )
  }

  blocks    <- partition$blocks %||% list()
  remainder <- .partition_remainder(partition)

  declared <- c(unlist(blocks, use.names = FALSE), remainder)
  missing_in_draws <- setdiff(declared, param_names)
  if (length(missing_in_draws)) {
    stop(
      "These partition parameters are not present in `draws_mat`:\n  ",
      paste(utils::head(missing_in_draws, 5L), collapse = ", "),
      if (length(missing_in_draws) > 5L) ", ..." else "",
      call. = FALSE
    )
  }
  not_in_partition <- setdiff(param_names, declared)
  if (length(not_in_partition)) {
    .message(
      "Blockwise compression: ", length(not_in_partition),
      " parameter(s) in `draws_mat` were not in the partition; ",
      "treating them as additional remainder."
    )
    remainder <- c(remainder, not_in_partition)
  }

  block_membership <- rep(NA_character_, length(param_names))
  names(block_membership) <- param_names

  fit_block <- function(name, members, model_name) {
    if (length(members) == 0L) return(NULL)
    block_draws <- draws_mat[, members, drop = FALSE]
    if (verbose) {
      .message(
        "blockwise: compressing '", name, "' (",
        length(members), " params)..."
      )
    }
    comp <- compress_posterior(
      block_draws,
      method       = method,
      n_components = n_components,
      model_name   = model_name,
      verbose      = verbose,
      partition    = NULL,
      ...
    )
    block_membership[members] <<- name
    comp
  }

  # If the caller did not specify model families, enforce the intended
  # constraint: clusters use non-diagonal covariance families (BIC within
  # that set), remainder uses diagonal families (BIC within that set).
  if (identical(method, "mclust")) {
    if (is.null(cluster_model_name)) {
      # Ellipsoidal (non-diagonal) family names. mclust will drop
      # unsupported ones as needed; our own mclust wrapper further
      # restricts when n <= d within each block.
      cluster_model_name <- c(
        "EEE", "VEE", "EVE", "VVE",
        "EEV", "VEV", "EVV", "VVV"
      )
    }
    if (is.null(remainder_model_name)) {
      # Diagonal/spherical families only (fast + identifiable in high-d).
      remainder_model_name <- c("EII", "VII", "EEI", "EVI", "VEI", "VVI")
    }
  }

  cluster_results <- lapply(
    names(blocks),
    function(nm) fit_block(nm, blocks[[nm]], cluster_model_name)
  )
  names(cluster_results) <- names(blocks)
  cluster_results <- cluster_results[!vapply(cluster_results, is.null, logical(1))]

  remainder_result <- fit_block("remainder", remainder, remainder_model_name)

  blocks_out <- cluster_results
  if (!is.null(remainder_result)) {
    blocks_out$remainder <- remainder_result
  }

  if (length(blocks_out) == 0L) {
    stop("Blockwise compression produced no blocks (empty partition?).",
         call. = FALSE)
  }

  out <- list(
    method              = "blockwise",
    base_method         = method,
    n_params            = length(param_names),
    n_components        = n_components,
    n_draws             = nrow(draws_mat),
    param_names         = param_names,
    block_membership    = block_membership,
    blocks              = blocks_out,
    cluster_block_names = setdiff(names(blocks_out), "remainder"),
    has_remainder       = !is.null(remainder_result)
  )
  class(out) <- c(
    "posterior_compressed_blockwise",
    "posterior_compressed",
    "list"
  )
  out
}


#' @keywords internal
#' @noRd
.sample_blockwise <- function(comp, n_draws = NULL, verbose = FALSE) {
  if (is.null(n_draws)) n_draws <- comp$n_draws
  n_draws <- as.integer(n_draws)
  if (length(n_draws) != 1L || is.na(n_draws) || n_draws <= 0L) {
    stop("`n_draws` must be a single positive integer.", call. = FALSE)
  }

  out <- matrix(
    NA_real_,
    nrow = n_draws,
    ncol = length(comp$param_names),
    dimnames = list(NULL, comp$param_names)
  )

  for (nm in names(comp$blocks)) {
    block <- comp$blocks[[nm]]
    if (isTRUE(verbose)) {
      .message(
        "blockwise sample: '", nm, "' (",
        length(block$param_names), " params, n_draws = ", n_draws, ")"
      )
    }
    block_samples <- sample_posterior(block, n_draws = n_draws)
    if (!identical(colnames(block_samples), block$param_names)) {
      block_samples <- block_samples[, block$param_names, drop = FALSE]
    }
    out[, block$param_names] <- block_samples
  }

  if (anyNA(out)) {
    missing_cols <- comp$param_names[
      apply(out, 2L, function(col) anyNA(col))
    ]
    stop(
      "Blockwise sampler left ", length(missing_cols),
      " parameter(s) unfilled (likely a partition / block mismatch): ",
      paste(utils::head(missing_cols, 5L), collapse = ", "),
      if (length(missing_cols) > 5L) ", ..." else "",
      call. = FALSE
    )
  }

  out
}


#' Density evaluation for a blockwise compressed posterior.
#'
#' Assumes the blocks are independent in the approximation, so
#' `log p(x) = sum_block log p_block(x_block)`.
#'
#' @keywords internal
#' @noRd
.density_blockwise <- function(comp, x, log = FALSE) {
  param_names <- comp$param_names
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    if (is.vector(x) && is.numeric(x)) {
      if (length(x) != length(param_names)) {
        stop(
          "x has length ", length(x), " but the compressed posterior has ",
          length(param_names), " parameters.",
          call. = FALSE
        )
      }
      x <- matrix(x, nrow = 1L, dimnames = list(NULL, param_names))
    } else {
      x <- as.matrix(x)
    }
  }

  if (is.null(colnames(x))) {
    if (ncol(x) != length(param_names)) {
      stop(
        "x has ", ncol(x), " columns but the compressed posterior has ",
        length(param_names), " parameters and no column names to align by.",
        call. = FALSE
      )
    }
    colnames(x) <- param_names
  } else {
    missing_pn <- setdiff(param_names, colnames(x))
    if (length(missing_pn)) {
      stop(
        "x is missing ", length(missing_pn),
        " parameter column(s) needed by the compressed posterior: ",
        paste(utils::head(missing_pn, 5L), collapse = ", "),
        if (length(missing_pn) > 5L) ", ..." else "",
        call. = FALSE
      )
    }
  }

  log_d <- numeric(nrow(x))
  for (block in comp$blocks) {
    sub <- x[, block$param_names, drop = FALSE]
    log_d <- log_d + density_posterior(block, sub, log = TRUE)
  }
  if (log) log_d else exp(log_d)
}


#' Extract a single block from a blockwise compressed posterior
#'
#' Convenience accessor: returns the per-block `posterior_compressed`
#' object for a named cluster (`"cluster_1"`, ...) or the literal name
#' `"remainder"`. Useful when you want to inspect or re-sample one
#' block in isolation (e.g. for diagnostics).
#'
#' @param x A `posterior_compressed_blockwise` object.
#' @param name Character; block name (one of `names(x$blocks)`).
#'
#' @return The matching `posterior_compressed` block, or an error if the
#'   name is not found.
#' @export
get_compressed_block <- function(x, name) {
  if (!inherits(x, "posterior_compressed_blockwise")) {
    stop(
      "`x` must be a posterior_compressed_blockwise object.",
      call. = FALSE
    )
  }
  if (!name %in% names(x$blocks)) {
    stop(
      "No block named '", name, "'. Available blocks: ",
      paste(names(x$blocks), collapse = ", "),
      call. = FALSE
    )
  }
  x$blocks[[name]]
}


#' @export
summary.posterior_compressed_blockwise <- function(object, ...) {
  blocks <- object$blocks
  bdf <- data.frame(
    block          = names(blocks),
    role           = ifelse(
      names(blocks) == "remainder", "remainder", "cluster"
    ),
    n_params       = vapply(blocks, function(b) {
      as.integer(b$n_params %||% length(b$param_names))
    }, integer(1L)),
    n_components   = vapply(blocks, function(b) {
      as.integer(b$n_components %||% NA_integer_)
    }, integer(1L)),
    model_name     = vapply(blocks, function(b) {
      as.character(b$model_name %||% NA_character_)
    }, character(1L)),
    covariance     = vapply(blocks, function(b) {
      as.character(b$covariance_type %||% NA_character_)
    }, character(1L)),
    stringsAsFactors = FALSE,
    row.names      = NULL
  )

  out <- list(
    method      = object$method,
    base_method = object$base_method,
    n_params    = object$n_params,
    n_draws     = object$n_draws,
    n_blocks    = length(blocks),
    n_clusters  = length(object$cluster_block_names),
    has_remainder = isTRUE(object$has_remainder),
    object_size = utils::object.size(object),
    blocks      = bdf
  )
  class(out) <- c(
    "summary.posterior_compressed_blockwise",
    "summary.posterior_compressed"
  )
  out
}


#' @export
print.summary.posterior_compressed_blockwise <- function(x, ...) {
  cat("Compressed posterior (blockwise)\n")
  cat("  base method   : ", x$base_method, "\n", sep = "")
  cat("  parameters    : ", x$n_params, "\n", sep = "")
  cat("  original draws: ", x$n_draws, "\n", sep = "")
  cat("  n_blocks      : ", x$n_blocks,
      "  (", x$n_clusters, " cluster(s)",
      if (x$has_remainder) " + remainder)" else ")",
      "\n", sep = "")
  cat("  in-memory size: ", format(x$object_size, units = "auto"), "\n",
      sep = "")
  cat("  block summary :\n")
  print(x$blocks, row.names = FALSE)
  invisible(x)
}


#' @export
print.posterior_compressed_blockwise <- function(x, ...) {
  cat("<posterior_compressed: blockwise>\n")
  cat(" base method     : ", x$base_method, "\n", sep = "")
  cat(" parameters      : ", x$n_params, "\n", sep = "")
  cat(" original draws  : ", x$n_draws, "\n", sep = "")
  cat(" cluster blocks  : ", length(x$cluster_block_names), "\n", sep = "")
  cat(" has remainder   : ", x$has_remainder, "\n", sep = "")
  cat(" block summary   :\n")
  for (nm in names(x$blocks)) {
    block <- x$blocks[[nm]]
    model <- if (!is.null(block$model_name)) {
      paste0("model=", block$model_name)
    } else {
      paste0("method=", block$method)
    }
    G <- if (!is.null(block$n_components)) {
      paste0(", G=", block$n_components)
    } else {
      ""
    }
    cat(sprintf(
      "    %-20s n=%-5d  %s%s\n",
      nm, block$n_params, model, G
    ))
  }
  invisible(x)
}
