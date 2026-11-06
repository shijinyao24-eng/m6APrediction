
#' Encode DNA Sequences into Categorical Features
#'
#' Converts DNA sequences into a data frame in which each nucleotide
#' position is represented as a categorical factor.
#'
#' @param dna_strings A character vector containing DNA sequences
#' of equal length, using A, T, C and G.
#'
#' @return A data frame in which each column represents one nucleotide
#' position and contains categorical factor values.
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Predict m6A Sites in Multiple Sequences
#'
#' Uses a fitted random forest model to predict m6A probabilities
#' and prediction statuses for multiple input sequences.
#'
#' @param ml_fit A fitted random forest classification model.
#' @param feature_df A data frame containing gc_content, RNA_type,
#' RNA_region, exon_length, distance_to_junction,
#' evolutionary_conservation and DNA_5mer.
#' @param positive_threshold A numeric threshold between 0 and 1
#' used to classify predicted m6A probabilities. The default is 0.5.
#'
#' @return A data frame containing the original input features,
#' predicted_m6A_prob and predicted_m6A_status.
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction"))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction"))
#' prediction_multiple(ml_fit, example_df, positive_threshold = 0.6)
#'
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}


#' Predict an Individual m6A Site
#'
#' Uses a fitted random forest model to predict the m6A probability
#' and prediction status for a single input sequence.
#'
#' @param ml_fit A fitted random forest classification model.
#' @param gc_content A numeric value representing the GC content.
#' @param RNA_type A character value specifying the RNA type.
#' Allowed values are mRNA, lincRNA, lncRNA and pseudogene.
#' @param RNA_region A character value specifying the RNA region.
#' Allowed values are CDS, intron, 3'UTR and 5'UTR.
#' @param exon_length A numeric value representing exon length.
#' @param distance_to_junction A numeric value representing
#' the distance to the nearest exon junction.
#' @param evolutionary_conservation A numeric value representing
#' evolutionary conservation.
#' @param DNA_5mer A character string containing a five-nucleotide
#' DNA sequence.
#' @param positive_threshold A numeric threshold between 0 and 1
#' used to classify predicted m6A probabilities. The default is 0.5.
#'
#' @return A named character vector containing predicted_m6A_prob
#' and predicted_m6A_status.
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction"))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8,
#'                   DNA_5mer = "GGACA", positive_threshold = 0.5)
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
