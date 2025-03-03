# ================================
# evaluate_method_performance.R
# ================================
# This file defines a unified calc_metrics() function to compute:
# - Average Credible Set (CS) Size
# - Coverage (proportion of CSs capturing a causal variant)
# - CS-based FDR and Recall
#
# The function adapts to different model outputs based on the provided method.

calc_metrics <- function(mod, method, X, causal) {
  # Define test.cs based on the method
  if (method == "susie") {
    test.cs <- susie_get_cs(mod, X = X, coverage = 0.95)$cs
  } else if (method == "fineboost") {
    test.cs <- mod$ucos_details$ucos$ucos_index
  } else if (method %in% c("susie_ash", "susie_inf")) {
    test.cs <- mod$sets
  } else {
    stop("Unknown method specified for metric calculation.")
  }

  # Initialize metrics
  cs_size   <- 0
  coverage  <- 0
  cs_fdr    <- 0
  cs_recall <- 0

  if (length(test.cs) > 0) {
    cs_size  <- length(unlist(test.cs)) / length(test.cs)
    coverage <- sum(sapply(test.cs, function(cs) any(causal %in% cs))) / length(test.cs)

    TP_fdr <- sum(sapply(test.cs, function(cs) any(cs %in% causal)))
    FP_fdr <- length(test.cs) - TP_fdr
    cs_fdr <- if ((TP_fdr + FP_fdr) > 0) FP_fdr / (TP_fdr + FP_fdr) else NA

    TP_recall <- sum(causal %in% unlist(test.cs))
    FN_recall <- length(causal) - TP_recall
    cs_recall <- TP_recall / (TP_recall + FN_recall)
  }

  return(list(
    cs_size   = cs_size,
    coverage  = coverage,
    cs_fdr    = cs_fdr,
    cs_recall = cs_recall
  ))
}

# Main function to aggregate metrics from all methods
evaluate_method_performance <- function(susie_out, susie_ash_out, susie_inf_out, fineboost_out, causal, data, precomp = NULL) {
  # data$X is used directly.
  susie_metrics      <- calc_metrics(susie_out,      method = "susie",      X = data$X, causal = causal)
  susie_ash_metrics  <- calc_metrics(susie_ash_out,  method = "susie_ash",  X = data$X, causal = causal)
  susie_inf_metrics  <- calc_metrics(susie_inf_out,  method = "susie_inf",  X = data$X, causal = causal)
  fineboost_metrics  <- calc_metrics(fineboost_out,  method = "fineboost",  X = data$X, causal = causal)

  # Create a summary metrics table
  metrics_table <- data.frame(
    Model     = c("SuSiE", "SuSiE.ash (Marginal)", "SuSiE-inf", "Fineboost"),
    CS_Size   = c(susie_metrics$cs_size,
                  susie_ash_metrics$cs_size,
                  susie_inf_metrics$cs_size,
                  fineboost_metrics$cs_size),
    Coverage  = c(susie_metrics$coverage,
                  susie_ash_metrics$coverage,
                  susie_inf_metrics$coverage,
                  fineboost_metrics$coverage),
    CS_FDR    = c(susie_metrics$cs_fdr,
                  susie_ash_metrics$cs_fdr,
                  susie_inf_metrics$cs_fdr,
                  fineboost_metrics$cs_fdr),
    CS_Recall = c(susie_metrics$cs_recall,
                  susie_ash_metrics$cs_recall,
                  susie_inf_metrics$cs_recall,
                  fineboost_metrics$cs_recall)
  )

  return(list(metrics = metrics_table))
}
