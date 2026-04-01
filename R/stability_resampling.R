# stability_resampling.R
# Utilities for generating perturbation plans used by stability analysis.

#' Validate counts/metadata inputs for stability analysis.
#' @param counts Integer matrix/data.frame: genes x samples.
#' @param metadata data.frame with sample and condition columns.
#' @param sample_col Metadata column containing sample IDs.
#' @param condition_col Metadata column containing condition labels.
#' @return List with aligned counts matrix and metadata table.
validate_stability_inputs <- function(counts,
                                      metadata,
                                      sample_col = "sample_id",
                                      condition_col = "condition") {
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("counts must be a matrix or data.frame with genes in rows and samples in columns")
  }
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame")
  }
  if (!sample_col %in% colnames(metadata)) {
    stop(sprintf("sample_col '%s' not found in metadata", sample_col))
  }
  if (!condition_col %in% colnames(metadata)) {
    stop(sprintf("condition_col '%s' not found in metadata", condition_col))
  }

  counts_mat <- as.matrix(counts)
  mode(counts_mat) <- "numeric"
  if (is.null(colnames(counts_mat))) {
    stop("counts must have sample IDs in column names")
  }

  metadata2 <- tibble::as_tibble(metadata)
  metadata2[[sample_col]] <- as.character(metadata2[[sample_col]])
  metadata2[[condition_col]] <- as.character(metadata2[[condition_col]])

  common_samples <- intersect(colnames(counts_mat), metadata2[[sample_col]])
  if (length(common_samples) < 4) {
    stop("Insufficient overlap between counts and metadata samples (need >= 4)")
  }

  metadata2 <- dplyr::filter(metadata2, .data[[sample_col]] %in% common_samples)
  metadata2 <- dplyr::arrange(metadata2, match(.data[[sample_col]], common_samples))
  counts_mat <- counts_mat[, metadata2[[sample_col]], drop = FALSE]

  counts_mat <- round(counts_mat)
  mode(counts_mat) <- "integer"

  group_sizes <- metadata2 |>
    dplyr::count(.data[[condition_col]], name = "n")

  if (nrow(group_sizes) != 2) {
    warning(sprintf(
      "Current module is optimized for two-group contrasts; detected %d groups.",
      nrow(group_sizes)
    ))
  }

  list(
    counts = counts_mat,
    metadata = metadata2,
    group_sizes = group_sizes
  )
}

#' Choose default perturbation modes from group sizes.
#' @param metadata Aligned metadata table.
#' @param condition_col Condition column name.
#' @param default_subsample_iter Number of subsampling iterations when enabled.
#' @return List with selected modes and rationale.
choose_stability_modes <- function(metadata,
                                   condition_col = "condition",
                                   default_subsample_iter = 50) {
  sizes <- metadata |>
    dplyr::count(.data[[condition_col]], name = "n")

  min_group_n <- min(sizes$n)
  if (min_group_n <= 4) {
    return(list(
      modes = c("loo"),
      subsample_iterations = 0L,
      rationale = sprintf(
        "At least one group has <= 4 samples (min=%d); defaulting to leave-one-out only.",
        min_group_n
      )
    ))
  }

  list(
    modes = c("loo", "subsample"),
    subsample_iterations = as.integer(default_subsample_iter),
    rationale = sprintf(
      "Group sizes are adequate (min=%d); running both leave-one-out and within-group subsampling.",
      min_group_n
    )
  )
}

#' Generate leave-one-out perturbation plans.
#' @param metadata Aligned metadata table.
#' @param sample_col Sample ID column.
#' @return tibble with one row per LOO iteration.
generate_loo_resamples <- function(metadata, sample_col = "sample_id") {
  samples <- as.character(metadata[[sample_col]])
  purrr::map_dfr(seq_along(samples), function(i) {
    removed <- samples[[i]]
    included <- setdiff(samples, removed)
    tibble::tibble(
      iteration_id = sprintf("loo_%03d", i),
      mode = "loo",
      removed_sample = removed,
      included_samples = list(included),
      excluded_samples = list(removed)
    )
  })
}

#' Generate repeated within-group subsampling plans.
#' @param metadata Aligned metadata table.
#' @param sample_col Sample ID column.
#' @param condition_col Condition label column.
#' @param fraction Fraction to keep per group in each iteration.
#' @param n_iter Number of subsampling iterations.
#' @param seed Random seed for reproducibility.
#' @return tibble with one row per subsampling iteration.
generate_subsample_resamples <- function(metadata,
                                         sample_col = "sample_id",
                                         condition_col = "condition",
                                         fraction = 0.8,
                                         n_iter = 50,
                                         seed = 123L) {
  if (fraction <= 0 || fraction >= 1) {
    stop("fraction must be in (0,1)")
  }
  if (n_iter < 1) {
    stop("n_iter must be >= 1")
  }

  group_df <- metadata |>
    dplyr::group_by(.data[[condition_col]]) |>
    dplyr::summarise(samples = list(as.character(.data[[sample_col]])), .groups = "drop")

  if (any(lengths(group_df$samples) < 3)) {
    warning("At least one group has < 3 samples; subsampling plans may be unstable or skipped.")
  }

  set.seed(seed)
  purrr::map_dfr(seq_len(n_iter), function(i) {
    sampled <- purrr::map(group_df$samples, function(samps) {
      k <- max(2L, floor(length(samps) * fraction))
      k <- min(k, length(samps))
      sort(sample(samps, size = k, replace = FALSE))
    })

    included <- sort(unique(unlist(sampled)))
    all_samples <- as.character(metadata[[sample_col]])
    excluded <- setdiff(all_samples, included)

    tibble::tibble(
      iteration_id = sprintf("subsample_%03d", i),
      mode = "subsample",
      removed_sample = NA_character_,
      included_samples = list(included),
      excluded_samples = list(excluded)
    )
  })
}

#' Build full perturbation plan for stability analysis.
#' @param metadata Aligned metadata table.
#' @param sample_col Sample ID column.
#' @param condition_col Condition column.
#' @param modes Character vector; any of c("loo","subsample").
#' @param subsample_fraction Within-group subsample fraction.
#' @param subsample_iterations Number of subsampling iterations.
#' @param seed Random seed.
#' @return tibble perturbation plan.
build_resample_plan <- function(metadata,
                                sample_col = "sample_id",
                                condition_col = "condition",
                                modes = c("loo", "subsample"),
                                subsample_fraction = 0.8,
                                subsample_iterations = 50,
                                seed = 123L) {
  mode_set <- unique(modes)
  plans <- list()

  if ("loo" %in% mode_set) {
    plans[[length(plans) + 1L]] <- generate_loo_resamples(
      metadata = metadata,
      sample_col = sample_col
    )
  }

  if ("subsample" %in% mode_set && subsample_iterations > 0) {
    plans[[length(plans) + 1L]] <- generate_subsample_resamples(
      metadata = metadata,
      sample_col = sample_col,
      condition_col = condition_col,
      fraction = subsample_fraction,
      n_iter = subsample_iterations,
      seed = seed
    )
  }

  if (length(plans) == 0) {
    stop("No perturbation plans generated; check modes and subsample_iterations settings")
  }

  dplyr::bind_rows(plans)
}
