compute_CNA_metrics_DLP <- function(confusion_df) {
  metrics_df <- confusion_df %>%
    mutate(
      precision = TP / (TP + FP),
      recall = TP / (TP + FN),
      f1 = 2 * precision * recall / (precision + recall)
    )
  
  # Handle NaNs
  metrics_df <- metrics_df %>%
    mutate(
      precision = ifelse(is.nan(precision), 0, precision),
      recall = ifelse(is.nan(recall), 0, recall),
      f1 = ifelse(is.nan(f1), 0, f1)
    )
  
  # Genome-wide metrics
  genome_total <- confusion_df %>%
    group_by(cell_id) %>% 
    summarise(
      TP = sum(TP),
      FP = sum(FP),
      FN = sum(FN),
      chr = "genome"
    ) %>%
    mutate(
      precision = TP / (TP + FP),
      recall = TP / (TP + FN),
      f1 = 2 * precision * recall / (precision + recall)
    ) %>%
    mutate(
      precision = ifelse(is.nan(precision), 0, precision),
      recall = ifelse(is.nan(recall), 0, recall),
      f1 = ifelse(is.nan(f1), 0, f1)
    )
  
  bind_rows(metrics_df, genome_total)
}

compute_confusion_matrix <- function(predicted, truth, delta = 10000, delta_merge = NULL) {
  # Step 0: Sort inputs
  predicted <- sort(predicted)
  truth <- sort(truth)
  
  # Step 1: Recursively merge predicted breakpoints
  if (!is.null(delta_merge) && length(predicted) > 1) {
    changed <- TRUE
    while (changed) {
      diffs <- diff(predicted)
      close_pairs <- which(diffs <= delta_merge)
      if (length(close_pairs) == 0) {
        changed <- FALSE
      } else {
        # Merge first close pair (can be done recursively or all at once, we'll go greedy)
        i <- close_pairs[1]
        new_val <- mean(predicted[i:(i+1)])
        predicted <- c(predicted[-c(i, i+1)], new_val)
        predicted <- sort(predicted)
      }
    }
  }
  
  # Step 2: Matching logic
  matched_truth <- logical(length(truth))
  matched_pred <- logical(length(predicted))
  
  for (i in seq_along(predicted)) {
    diffs <- abs(truth - predicted[i])
    # Find unmatched truth within delta
    candidates <- which(diffs <= delta & !matched_truth)
    if (length(candidates) > 0) {
      best_match <- candidates[which.min(diffs[candidates])]
      matched_truth[best_match] <- TRUE
      matched_pred[i] <- TRUE
    }
  }
  
  TP <- sum(matched_pred)
  FP <- sum(!matched_pred)
  FN <- sum(!matched_truth)
  
  dplyr::tibble(TP = TP, FP = FP, FN = FN)
}
