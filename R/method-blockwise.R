#' Declare manual correlation clusters for blockwise compression
#'
#' Captures block specifications with [rlang::enquos()] so you can mix
#' literal parameter names and **tidy-select** expressions (same semantics as
#' column selection, e.g. [tidyselect::starts_with()], [tidyselect::contains()]).
#' Parameters matched by no block are placed in **`remainder`** when you pass
#' the result as `partition` to [compress_posterior()] or [compress_brmsfit()].
#'
#' For **character-only** blocks you can instead pass a plain `list()` of
#' character vectors; that form is coerced automatically.
#'
#' **Tidy-select** helpers must be captured: use `partition_blocks(...)` or,
#' in a plain list, a **formula** on the RHS such as
#' `~ tidyselect::starts_with("b_")`. Writing `list(tidyselect::starts_with("b_"))`
#' does **not** work, because R evaluates `starts_with()` before compression
#' and the result is not a column selection expression.
#'
#' Overlap is resolved **in order**: a parameter is assigned to the first
#' block that matches it; later blocks do not steal names.
#'
#' **Singleton blocks** (only one parameter after resolution) are not used as
#' correlation clusters: those parameters are moved to **`remainder`** with a
#' [warning()], because within-block correlations are undefined for a single
#' column.
#'
#' @param ... One or more block specifications: a character vector of names,
#'   a `list()` of character vectors, a **formula** with tidy-select on the
#'   RHS (e.g. `~ tidyselect::starts_with("b_")`), or any expression
#'   understood by [tidyselect::eval_select()] (including bare
#'   `tidyselect::starts_with("x")` when captured here).
#'
#' @return An object of class `"poco_partition_block_spec"` (a list of
#'   quosures). Normally passed directly to `compress_*(..., partition = )`;
#'   resolution against actual parameter names happens during compression.
#'
#' @seealso [compress_posterior()], [compress_brmsfit()],
#'   [partition_parameters_clusters()] (default **`simple_output = TRUE`**
#'   plain list, or **`simple_output = FALSE`** for the full object).
#' @export
#' @examples
#' \dontrun{
#' draws <- matrix(rnorm(100 * 5), ncol = 5)
#' colnames(draws) <- c("a1", "b_1", "b_2", "dataset_id_x", "z")
#' compress_posterior(
#'   draws,
#'   method = "mclust",
#'   n_components = 2L,
#'   partition = partition_blocks(
#'     c("a1", "z"),
#'     tidyselect::starts_with("b_"),
#'     tidyselect::contains("dataset_id")
#'   )
#' )
#' # Plain list of name vectors only:
#' compress_posterior(draws, method = "mclust", partition = list(c("a1"), c("b_1", "b_2")))
#' }
partition_blocks <- function(...) {
  structure(rlang::enquos(...), class = c("poco_partition_block_spec", "list"))
}


#' @keywords internal
#' @noRd
.partition_sel_df <- function(param_names) {
  as.data.frame(
    matrix(
      ncol = length(param_names),
      nrow = 0L,
      dimnames = list(NULL, param_names)
    ),
    check.names = FALSE
  )
}


#' @keywords internal
#' @noRd
.block_partition_looks_reserved <- function(x) {
  nm <- names(x)
  if (!length(nm) || !all(nzchar(nm))) {
    return(FALSE)
  }
  any(
    c("blocks", "cluster_id", "remainder", "hclust", "params") %in% nm
  )
}


#' @keywords internal
#' @noRd
.resolve_partition_block_spec <- function(spec, param_names, sel_df) {
  if (rlang::is_quosure(spec)) {
    return(names(tidyselect::eval_select(spec, data = sel_df, allow_rename = FALSE)))
  }
  if (inherits(spec, "formula")) {
    rq <- rlang::new_quosure(
      rlang::f_rhs(spec),
      env = rlang::f_env(spec)
    )
    return(names(tidyselect::eval_select(rq, data = sel_df, allow_rename = FALSE)))
  }
  if (is.character(spec)) {
    bad <- setdiff(spec, param_names)
    if (length(bad)) {
      warning(
        "Ignoring ", length(bad), " partition name(s) not found in draws: ",
        paste(utils::head(bad, 10L), collapse = ", "),
        if (length(bad) > 10L) ", ..." else "",
        call. = FALSE
      )
    }
    return(spec[spec %in% param_names])
  }
  if (is.list(spec) && length(spec)) {
    all_chr <- all(vapply(spec, is.character, logical(1L)))
    if (all_chr) {
      v <- unlist(spec, use.names = FALSE)
      bad <- setdiff(v, param_names)
      if (length(bad)) {
        warning(
          "Ignoring ", length(bad), " partition name(s) not found in draws: ",
          paste(utils::head(bad, 10L), collapse = ", "),
          if (length(bad) > 10L) ", ..." else "",
          call. = FALSE
        )
      }
      return(v[v %in% param_names])
    }
  }
  stop(
    "Invalid manual `partition` block (element ",
    "must be a character vector, list of character vectors, formula, ",
    "or quosure; use partition_blocks() for bare tidy-select helpers).",
    call. = FALSE
  )
}


#' @keywords internal
#' @noRd
.partition_from_block_list <- function(block_specs, param_names) {
  if (!length(param_names)) {
    stop("`param_names` is empty.", call. = FALSE)
  }
  sel_df <- .partition_sel_df(param_names)
  assigned <- character()
  blocks <- list()
  for (i in seq_along(block_specs)) {
    spec <- block_specs[[i]]
    nm <- .resolve_partition_block_spec(spec, param_names, sel_df)
    nm <- unique(nm[!is.na(nm) & nzchar(nm)])
    nm <- nm[nm %in% param_names]
    nm <- setdiff(nm, assigned)
    if (!length(nm)) {
      next
    }
    assigned <- c(assigned, nm)
    blocks[[paste0("cluster_", length(blocks) + 1L)]] <- nm
  }
  remainder <- setdiff(param_names, assigned)
  cluster_id <- rep(NA_integer_, length(param_names))
  names(cluster_id) <- param_names
  for (k in seq_along(blocks)) {
    cluster_id[blocks[[k]]] <- as.integer(k)
  }
  out <- list(
    blocks     = blocks,
    remainder  = remainder,
    cluster_id = cluster_id,
    hclust     = NULL,
    params     = list(source = "partition_manual")
  )
  class(out) <- c("poco_partition_clusters", "list")
  out
}


#' @keywords internal
#' @noRd
.coerce_blockwise_partition <- function(partition, param_names) {
  if (is.null(partition)) {
    return(NULL)
  }
  if (inherits(partition, "poco_partition_clusters")) {
    return(partition)
  }
  if (inherits(partition, "poco_partition_block_spec")) {
    specs <- unclass(partition)
    return(.partition_from_block_list(specs, param_names))
  }
  # Includes `list()` (length 0): all parameters go to remainder (e.g.
  # `partition_parameters_clusters(..., simple_output = TRUE)` when no cluster
  # passes `min_size`).
  if (is.list(partition) && !.block_partition_looks_reserved(partition)) {
    return(.partition_from_block_list(partition, param_names))
  }
  partition
}


#' Drop partition names absent from draws; add draw columns missing from the
#' partition to remainder. Emits [warning()] for both cases.
#'
#' @keywords internal
#' @noRd
.align_partition_to_draws <- function(partition, param_names) {
  if (!inherits(partition, "poco_partition_clusters")) {
    return(partition)
  }
  blocks <- partition$blocks %||% list()
  remainder_orig <- .partition_remainder(partition)
  in_blocks_orig <- unlist(blocks, use.names = FALSE)
  all_declared <- unique(c(in_blocks_orig, remainder_orig))
  bad <- setdiff(all_declared, param_names)
  if (length(bad)) {
    warning(
      "Ignoring ", length(bad), " partition name(s) not present in draws: ",
      paste(utils::head(bad, 10L), collapse = ", "),
      if (length(bad) > 10L) ", ..." else "",
      call. = FALSE
    )
  }

  new_blocks <- list()
  for (i in seq_along(blocks)) {
    nm <- blocks[[i]]
    nm <- intersect(nm, param_names)
    nm <- unique(nm[nzchar(nm)])
    if (length(nm)) {
      new_blocks[[length(new_blocks) + 1L]] <- nm
    }
  }
  if (length(new_blocks)) {
    names(new_blocks) <- paste0("cluster_", seq_along(new_blocks))
  } else {
    new_blocks <- list()
  }

  assigned <- unlist(new_blocks, use.names = FALSE)
  remainder <- intersect(remainder_orig, param_names)
  remainder <- unique(remainder[nzchar(remainder)])
  remainder <- setdiff(remainder, assigned)

  missing_fit <- setdiff(param_names, c(assigned, remainder))
  if (length(missing_fit)) {
    warning(
      length(missing_fit),
      " draw column(s) were not listed in any cluster or remainder; ",
      "assigning to remainder: ",
      paste(utils::head(missing_fit, 10L), collapse = ", "),
      if (length(missing_fit) > 10L) ", ..." else "",
      call. = FALSE
    )
    remainder <- c(remainder, missing_fit)
  }

  cluster_id <- rep(NA_integer_, length(param_names))
  names(cluster_id) <- param_names
  for (k in seq_along(new_blocks)) {
    cluster_id[new_blocks[[k]]] <- as.integer(k)
  }

  partition$blocks     <- new_blocks
  partition$remainder  <- remainder
  partition$cluster_id <- cluster_id
  partition
}


#' Move singleton-parameter cluster blocks into remainder (blockwise only).
#'
#' A cluster block with a single parameter cannot support within-block
#' correlations; such blocks are dropped and those parameters are appended to
#' `remainder`. Emits one [warning()] listing the moved parameters.
#'
#' @keywords internal
#' @noRd
.partition_demote_singleton_cluster_blocks <- function(partition, param_names) {
  if (!inherits(partition, "poco_partition_clusters")) {
    return(partition)
  }
  blocks <- partition$blocks %||% list()
  if (!length(blocks)) {
    return(partition)
  }

  singletons <- character()
  kept <- list()
  for (i in seq_along(blocks)) {
    nm <- blocks[[i]]
    if (length(nm) > 1L) {
      kept[[length(kept) + 1L]] <- nm
    } else if (length(nm) == 1L) {
      singletons <- c(singletons, nm)
    }
  }

  if (!length(singletons)) {
    return(partition)
  }

  singletons <- unique(singletons[singletons %in% param_names])
  rem_prev <- .partition_remainder(partition)
  remainder <- unique(c(intersect(rem_prev, param_names), singletons))

  if (length(kept)) {
    names(kept) <- paste0("cluster_", seq_along(kept))
  }

  warning(
    "Blockwise partition: ",
    length(singletons),
    " cluster block(s) contained only one parameter ",
    "(within-block correlations are undefined); ",
    "those parameters were moved to remainder: ",
    paste(utils::head(singletons, 15L), collapse = ", "),
    if (length(singletons) > 15L) ", ..." else "",
    call. = FALSE
  )

  cluster_id <- rep(NA_integer_, length(param_names))
  names(cluster_id) <- param_names
  for (k in seq_along(kept)) {
    cluster_id[kept[[k]]] <- as.integer(k)
  }

  partition$blocks <- kept
  partition$remainder <- remainder
  partition$cluster_id <- cluster_id
  partition
}


#' Upper bound on parallel workers from cluster-block count and CPUs.
#'
#' @keywords internal
#' @noRd
.bp_cluster_worker_cap <- function(n_cluster_blocks) {
  nb <- as.integer(n_cluster_blocks)
  if (length(nb) != 1L || is.na(nb) || nb < 1L) {
    return(1L)
  }
  dc <- parallel::detectCores()
  if (is.na(dc) || dc < 2L) {
    dc <- 2L
  }
  min(nb, max(1L, as.integer(dc) - 1L))
}


#' Build `SnowParam` / `MulticoreParam` for cluster-block `bplapply()`.
#'
#' @return A `BiocParallelParam`, or `NULL` if **BiocParallel** is not installed.
#'
#' @keywords internal
#' @noRd
.bp_new_cluster_parallel_param <- function(workers, verbose = FALSE) {
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    return(NULL)
  }
  workers <- max(1L, as.integer(workers))
  .bp_mk_param <- function(constructor, w, show_bar) {
    args <- list(workers = as.integer(w))
    if (isTRUE(show_bar)) {
      args$progressbar <- TRUE
    }
    tryCatch(
      do.call(constructor, args),
      error = function(e) do.call(constructor, list(workers = args$workers))
    )
  }
  if (identical(.Platform$OS.type, "windows")) {
    .bp_mk_param(BiocParallel::SnowParam, workers, verbose)
  } else {
    .bp_mk_param(BiocParallel::MulticoreParam, workers, verbose)
  }
}


#' Resolve parallel backend for cluster-block compression.
#'
#' @param cluster_BPPARAM `NULL` (auto), `FALSE` (sequential), a single
#'   non-negative numeric worker count (e.g. `2` or `2L`, coerced with
#'   [ceiling()] when non-integer), or a `BiocParallelParam`.
#' @param n_cluster_blocks Number of cluster blocks (not counting remainder).
#'
#' @return `NULL` for sequential [lapply()], else a `BiocParallelParam` for
#'   [BiocParallel::bplapply()].
#'
#' @param verbose Passed to auto-built `SnowParam` / `MulticoreParam` as
#'   `progressbar = verbose` when supported (BiocParallel progress text).
#'
#' @keywords internal
#' @noRd
.resolve_cluster_BPPARAM <- function(
    cluster_BPPARAM,
    n_cluster_blocks,
    verbose = FALSE) {
  if (isFALSE(cluster_BPPARAM)) {
    return(NULL)
  }

  if (!is.null(cluster_BPPARAM)) {
    if (methods::is(cluster_BPPARAM, "BiocParallelParam")) {
      if (!requireNamespace("BiocParallel", quietly = TRUE)) {
        stop(
          "Install BiocParallel to use a custom parallel backend ",
          "(BiocManager::install(\"BiocParallel\")).",
          call. = FALSE
        )
      }
      return(cluster_BPPARAM)
    }

    if (length(cluster_BPPARAM) == 1L &&
        is.numeric(cluster_BPPARAM) &&
        !is.logical(cluster_BPPARAM)) {
      raw_disp <- cluster_BPPARAM
      raw <- suppressWarnings(as.numeric(cluster_BPPARAM))
      if (is.na(raw) || !is.finite(raw)) {
        .message(
          "cluster_BPPARAM = ", raw_disp,
          " is not a finite number; using sequential cluster compression."
        )
        return(NULL)
      }
      if (raw < 2) {
        .message(
          "cluster_BPPARAM = ", raw_disp,
          " requests fewer than 2 parallel workers; using sequential cluster ",
          "compression (same as cluster_BPPARAM = FALSE)."
        )
        return(NULL)
      }

      req <- as.integer(ceiling(raw))
      cap <- .bp_cluster_worker_cap(n_cluster_blocks)
      w_final <- min(req, cap)

      kind <- if (identical(.Platform$OS.type, "windows")) {
        "BiocParallel::SnowParam"
      } else {
        "BiocParallel::MulticoreParam"
      }

      if (w_final < 2L) {
        .message(
          "cluster_BPPARAM = ", raw_disp,
          " (ceiling = ", req, " worker",
          if (!identical(req, 1L)) "s" else "",
          "): fewer than 2 parallel workers are available (cluster blocks = ",
          n_cluster_blocks, ", worker cap = ", cap,
          "). Using sequential cluster compression."
        )
        return(NULL)
      }

      parts <- character()
      if (abs(raw - as.double(req)) > 1e-12) {
        parts <- c(
          parts,
          paste0("ceiling(", raw_disp, ") = ", req, " worker",
                 if (!identical(req, 1L)) "s" else "")
        )
      }
      if (w_final < req) {
        parts <- c(
          parts,
          paste0("capped ", req, " to ", w_final, " (max ", cap,
                 " from cluster blocks / CPUs)")
        )
      }

      bp <- .bp_new_cluster_parallel_param(w_final, verbose)
      if (is.null(bp)) {
        message(
          "BiocParallel is not installed; compressing cluster blocks sequentially. ",
          "Install with BiocManager::install(\"BiocParallel\") to enable parallel ",
          "cluster compression."
        )
        return(NULL)
      }

      msg_core <- paste0(
        "cluster_BPPARAM = ", raw_disp, ": using ", w_final, " worker",
        if (!identical(w_final, 1L)) "s" else "",
        " via ", kind, "."
      )
      if (length(parts)) {
        .message(msg_core, " (", paste(parts, collapse = "; "), ")")
      } else {
        .message(msg_core)
      }
      return(bp)
    }

    stop(
      "`cluster_BPPARAM` must be NULL (auto), FALSE (sequential), ",
      "a single numeric worker count (e.g. 2 or 2L), ",
      "or a BiocParallelParam object.",
      call. = FALSE
    )
  }

  if (n_cluster_blocks < 2L) {
    return(NULL)
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    message(
      "BiocParallel is not installed; compressing cluster blocks sequentially. ",
      "Install with BiocManager::install(\"BiocParallel\") to enable parallel ",
      "cluster compression (worker count is chosen automatically from the ",
      "number of clusters and available CPUs when cluster_BPPARAM = NULL)."
    )
    return(NULL)
  }
  w <- .bp_cluster_worker_cap(n_cluster_blocks)
  if (w < 2L) {
    return(NULL)
  }
  .bp_new_cluster_parallel_param(w, verbose)
}


#' Copy a BiocParallelParam with `progressbar = TRUE` when possible (does not
#' mutate the caller's object).
#'
#' @keywords internal
#' @noRd
.bp_enable_progress_copy <- function(bp, enable) {
  if (!isTRUE(enable) || is.null(bp)) {
    return(bp)
  }
  if (!methods::is(bp, "BiocParallelParam")) {
    return(bp)
  }
  if (exists("bpprogressbar", envir = asNamespace("BiocParallel"), inherits = FALSE)) {
    already <- tryCatch(
      isTRUE(BiocParallel::bpprogressbar(bp)),
      error = function(e) FALSE
    )
    if (already) {
      return(bp)
    }
    bp2 <- bp
    ok <- tryCatch({
      BiocParallel::bpprogressbar(bp2) <- TRUE
      TRUE
    }, error = function(e) FALSE)
    if (ok) {
      return(bp2)
    }
  }
  sn <- tryCatch(methods::slotNames(bp), error = function(e) character())
  if (!("progressbar" %in% sn)) {
    return(bp)
  }
  if (isTRUE(methods::slot(bp, "progressbar"))) {
    return(bp)
  }
  bp2 <- bp
  ok <- tryCatch({
    methods::slot(bp2, "progressbar") <- TRUE
    TRUE
  }, error = function(e) FALSE)
  if (ok) bp2 else bp
}


#' Integer worker count for messaging.
#'
#' @keywords internal
#' @noRd
.bp_workers_n <- function(bp) {
  if (is.null(bp) || !requireNamespace("BiocParallel", quietly = TRUE)) {
    return(NA_integer_)
  }
  tryCatch(
    as.integer(BiocParallel::bpworkers(bp)),
    error = function(e) NA_integer_
  )
}


#' Whether BiocParallel progress bar is enabled on a backend object.
#'
#' @keywords internal
#' @noRd
.bp_progress_enabled <- function(bp) {
  if (is.null(bp) || !requireNamespace("BiocParallel", quietly = TRUE)) {
    return(FALSE)
  }
  if (exists("bpprogressbar", envir = asNamespace("BiocParallel"), inherits = FALSE)) {
    val <- tryCatch(BiocParallel::bpprogressbar(bp), error = function(e) NULL)
    if (!is.null(val)) {
      return(isTRUE(val))
    }
  }
  sn <- tryCatch(methods::slotNames(bp), error = function(e) character())
  if ("progressbar" %in% sn) {
    return(isTRUE(tryCatch(methods::slot(bp, "progressbar"), error = function(e) FALSE)))
  }
  FALSE
}


#' Blockwise (per-cluster) posterior compression
#'
#' Internal worker behind the `partition` argument of [compress_posterior()]
#' and [compress_brmsfit()].
#'
#' Each cluster of `partition$blocks` is compressed independently with a
#' typically richer covariance family (`cluster_model_name`, default
#' ellipsoidal models with BIC selection). If that fit errors (e.g. block
#' dimension too large), the same block is retried once with spherical /
#' diagonal models when `cluster_model_name` was left at `NULL`. The leftover
#' `partition$remainder` parameters are compressed together with a cheaper
#' diagonal family (`remainder_model_name`, default spherical/diagonal set).
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
    cluster_BPPARAM      = NULL,
    ...) {
  param_names <- colnames(draws_mat)
  if (is.null(param_names)) {
    stop("Blockwise compression requires `draws_mat` to have column names.",
         call. = FALSE)
  }

  partition <- .coerce_blockwise_partition(partition, param_names)

  if (!inherits(partition, "poco_partition_clusters")) {
    stop(
      "`partition` must be coercible: output of ",
      "`partition_parameters_clusters()` (default plain list or full object), ",
      "`partition_blocks()`, or a plain `list()` of character vectors ",
      "(see `?partition_blocks`; use `list()` for an all-remainder block).",
      call. = FALSE
    )
  }

  partition <- .align_partition_to_draws(partition, param_names)
  partition <- .partition_demote_singleton_cluster_blocks(partition, param_names)
  blocks    <- partition$blocks %||% list()
  remainder <- .partition_remainder(partition)

  block_membership <- rep(NA_character_, length(param_names))
  names(block_membership) <- param_names

  extra_args <- list(...)

  # Spherical/diagonal mclust families (BIC among these). Used for remainder
  # by default and as a cluster-block fallback when default ellipsoidal
  # cluster fits fail (e.g. very large blocks).
  mclust_diag_family <- c("EII", "VII", "EEI", "EVI", "VEI", "VVI")
  cluster_mclust_diag_fallback <- FALSE

  # If the caller did not specify model families, enforce the intended
  # constraint: clusters use non-diagonal covariance families (BIC within
  # that set), remainder uses diagonal families (BIC within that set).
  if (identical(method, "mclust")) {
    cluster_mclust_diag_fallback <- is.null(cluster_model_name)
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
      remainder_model_name <- mclust_diag_family
    }
  }

  block_nms <- names(blocks)
  bp <- .resolve_cluster_BPPARAM(
    cluster_BPPARAM,
    length(block_nms),
    verbose = isTRUE(verbose)
  )
  bp <- .bp_enable_progress_copy(bp, isTRUE(verbose))
  parallel_clusters <- !is.null(bp)
  if (isTRUE(verbose) && parallel_clusters) {
    nw <- .bp_workers_n(bp)
    has_pb <- .bp_progress_enabled(bp)
    wpart <- if (!is.na(nw)) {
      paste0(" using ", nw, " worker", if (!identical(nw, 1L)) "s" else "")
    } else {
      ""
    }
    pbpart <- if (has_pb) {
      " Progress bar enabled."
    } else {
      " Progress bar not enabled for this backend object."
    }
    core <- paste0(
      "blockwise: compressing ", length(block_nms),
      " cluster block(s) in parallel with BiocParallel",
      wpart
    )
    .message(paste0(core, ".", pbpart))
  } else if (isTRUE(verbose) && length(block_nms) >= 2L) {
    .message(
      "blockwise: compressing ", length(block_nms),
      " cluster block(s) sequentially."
    )
  }

  compress_one_block <- function(
      name,
      members,
      model_name,
      use_cluster_diag_fallback = FALSE) {
    if (length(members) == 0L) {
      return(NULL)
    }
    block_draws <- draws_mat[, members, drop = FALSE]
    announce <- isTRUE(verbose) &&
      (!parallel_clusters || identical(name, "remainder"))
    if (announce) {
      .message(
        "blockwise: compressing '", name, "' (",
        length(members), " params)..."
      )
    }
    run_compress <- function(mn) {
      do.call(
        compress_posterior,
        c(
          list(
            draws        = block_draws,
            method       = method,
            n_components = n_components,
            model_name   = mn,
            # Keep backend fitting quiet; blockwise progress is reported here.
            verbose      = FALSE,
            partition    = NULL
          ),
          extra_args
        )
      )
    }
    comp <- if (isTRUE(use_cluster_diag_fallback)) {
      tryCatch(
        run_compress(model_name),
        error = function(e) {
          # New line so this does not append to an active BiocParallel progress bar.
          cat("\n", file = stderr())
          .message(
            "blockwise: mclust failed for cluster block \"", name, "\" (",
            length(members), " parameters) using ellipsoidal covariance ",
            "models. Retrying with spherical/diagonal models (BIC within ",
            paste(mclust_diag_family, collapse = ", "), "). ",
            "Original error: ", conditionMessage(e),
            "\n"
          )
          run_compress(mclust_diag_family)
        }
      )
    } else {
      run_compress(model_name)
    }
    list(name = name, members = members, comp = comp)
  }

  if (is.null(bp)) {
    fitted_blocks <- vector("list", length(block_nms))
    names(fitted_blocks) <- block_nms
    for (i in seq_along(block_nms)) {
      nm <- block_nms[[i]]
      fitted_blocks[[nm]] <- compress_one_block(
        nm,
        blocks[[nm]],
        cluster_model_name,
        use_cluster_diag_fallback = cluster_mclust_diag_fallback
      )
      if (isTRUE(verbose)) {
        .message(
          "blockwise: fitted cluster ", i, "/", length(block_nms),
          " (", nm, ")."
        )
      }
    }
  } else {
    fitted_blocks <- BiocParallel::bplapply(
      block_nms,
      function(nm) {
        compress_one_block(
          nm,
          blocks[[nm]],
          cluster_model_name,
          use_cluster_diag_fallback = cluster_mclust_diag_fallback
        )
      },
      BPPARAM = bp
    )
    names(fitted_blocks) <- block_nms
  }

  cluster_results <- list()
  for (nm in block_nms) {
    r <- fitted_blocks[[nm]]
    if (is.null(r)) {
      next
    }
    block_membership[r$members] <- r$name
    cluster_results[[r$name]] <- r$comp
  }

  r_rem <- compress_one_block("remainder", remainder, remainder_model_name)
  remainder_result <- NULL
  if (!is.null(r_rem)) {
    block_membership[r_rem$members] <- r_rem$name
    remainder_result <- r_rem$comp
  }

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
