#' @title Calculate log2 fold change
#'
#' @description This function calculates log2(FC), p-values, and adjusted p-values
#' of the data using limma.
#' @param metalyzer_se A Metalyzer object
#' @param categorical A column specifying the two groups
#' @param impute_perc_of_min A numeric value below 1
#' @param impute_NA Logical value whether to impute NA values
#'
#' @return A data frame containing the log2 fold change for each metabolite
#'
#' @import dplyr
#' @import SummarizedExperiment
#' @import limma
#' @importFrom rlang .data
#' @importFrom data.table :=
#' @importFrom stats model.matrix
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer_dataset(file_path = example_mutation_data_xl())
#' metalyzer_se <- filterMetabolites(
#'   metalyzer_se,
#'   drop_metabolites = "Metabolism Indicators"
#' )
#' metalyzer_se <- renameMetaData(
#'   metalyzer_se,
#'   Mutant_Control = "Sample Description"
#' )
#' 
#' metalyzer_se <- calculate_log2FC(
#'   metalyzer_se,
#'   categorical = "Mutant_Control",
#'   impute_perc_of_min = 0.2,
#'   impute_NA = FALSE
#' )
calculate_log2FC <- function(metalyzer_se,
                             categorical,
                             impute_perc_of_min = 0.2,
                             impute_NA = FALSE) {
  
  ## Create a new dataframe to calculate the log2FC
  if (!categorical %in% colnames(metalyzer_se@metadata$aggregated_data)) {
    metalyzer_se <- expand_aggregated_data(metalyzer_se,
                                           meta_data_column = categorical)
  }
  
  metalyzer_se <- data_imputation(metalyzer_se, impute_perc_of_min, impute_NA)
  metalyzer_se <- data_transformation(metalyzer_se)

  aggregated_data <- metalyzer_se@metadata$aggregated_data  %>%
    dplyr::mutate(LogConcentration = log2(.data$Concentration))
    # Prepare abundance data
    feat_data <- dplyr::ungroup(aggregated_data) %>%
      dplyr::select(.data$Metabolite, .data$ID, .data$LogConcentration) %>%
      tidyr::pivot_wider(names_from = 'Metabolite', values_from = 'LogConcentration')

    # Prepare sample metadata of interest
    smp_metadata <- colData(metalyzer_se) %>%
      tibble::as_tibble(rownames = 'ID')

    # Use original column names whose spaces are not replaced with '.'
    colnames(smp_metadata) <- c('ID', colnames(colData(metalyzer_se)))
    smp_metadata <- dplyr::select(smp_metadata, .data$ID, all_of(categorical))

    # Combine abundance data and sample metadata to ensure matched information
    combined_data <- dplyr::left_join(feat_data, smp_metadata, by = 'ID')

    # Retrieve data matrix and sample metadata from combined data to conduct limma
    data_mat <- combined_data[, seq_len(ncol(feat_data))] %>%
      tibble::column_to_rownames('ID') %>%
      t()
    group_vec <- combined_data[, ncol(feat_data)+1, drop = T]

    # Sanity check if specified categorical can split data into two groups
    if (length(unique(group_vec)) != 2) {
      stop("The specified categorical cannot split data into two groups.")
    }
    
    ## Compute log2(FC), p-values, and adjusted p-values using limma
    design <- model.matrix(~ group_vec)
    fit <- limma::lmFit(data_mat, design = design)
    fit <- limma::eBayes(fit)
    log2FCRes <- limma::topTable(fit, coef = 2, number = Inf) %>%
      tibble::rownames_to_column('Metabolite') %>%
      dplyr::select(.data$Metabolite, .data$logFC, .data$P.Value, .data$adj.P.Val) %>%
      dplyr::rename(log2FC = .data$logFC, pval = .data$P.Value, qval = .data$adj.P.Val)
    # Combined all information into a table
    group_info <- combined_data[, c(1, ncol(feat_data)+1)]
    log2FCTab <- dplyr::left_join(aggregated_data, group_info, by = 'ID') %>%
      dplyr::left_join(log2FCRes, by = 'Metabolite') %>%
      dplyr::select(.data$Metabolite, .data$Class, .data$log2FC, .data$pval, .data$qval) %>%
      dplyr::distinct(.data$Metabolite, .keep_all = TRUE)
    metalyzer_se@metadata$log2FC <- log2FCTab

    return(metalyzer_se)
}