#' Encode nucleotide positions as factors
#' 
#' Split equal-length DNA strings into nucleotide positions.
#' @param dna_strings Character vector of equal-length uppercase DNA strings (A, T, C, G).
#' @return A data frame with one row per string and a factor column per position.
dna_encoding <- function (dna_strings) 
{
    nn <- nchar(dna_strings[1])
    seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
    colnames(seq_m) <- paste0("nt_pos", 1:nn)
    seq_df <- as.data.frame(seq_m)
    seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
    return(seq_df)
}

#' Predict m6A status for multiple sites
#' 
#' Encode nucleotide contexts and obtain Positive-class scores for a feature table.
#' @param ml_fit Fitted course randomForest classification model.
#' @param feature_df Data frame containing gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation and DNA_5mer.
#' @param positive_threshold Numeric cutoff. Only scores strictly above it are Positive.
#' @return Input columns plus numeric predicted_m6A_prob and character predicted_m6A_status.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction", mustWork = TRUE))
#' features <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                  package = "m6APrediction", mustWork = TRUE))
#' result <- prediction_multiple(ml_fit, head(features, 3), positive_threshold = 0.6)
#' print(result[, c("predicted_m6A_prob", "predicted_m6A_status")])
prediction_multiple <- function (ml_fit, feature_df, positive_threshold = 0.5) 
{
    stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", 
        "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
    feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", 
        "pseudogene"))
    feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", 
        "5'UTR"))
    infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), 
        type = "prob")
    feature_df$predicted_m6A_prob <- infered_prob[, "Positive"]
    feature_df$predicted_m6A_status <- "Negative"
    feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
    return(feature_df)
}

#' Predict m6A status for one site
#' 
#' Score one candidate site using the feature values supplied.
#' @param ml_fit Fitted course randomForest classification model.
#' @param gc_content Numeric GC fraction.
#' @param RNA_type RNA category: mRNA, lincRNA, lncRNA or pseudogene.
#' @param RNA_region Region category: CDS, intron, 3'UTR or 5'UTR.
#' @param exon_length Numeric exon length on the course log2 scale.
#' @param distance_to_junction Numeric junction distance on the course log2 scale.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer Uppercase five-base context with A in the central position.
#' @param positive_threshold Numeric cutoff. Only scores strictly above it are Positive.
#' @return Named character vector with predicted_m6A_prob and predicted_m6A_status.
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
#'                   exon_length = 12, distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8, DNA_5mer = "GGACA",
#'                   positive_threshold = 0.5)
prediction_single <- function (ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, 
    DNA_5mer, positive_threshold = 0.5) 
{
    pred_outcome <- prediction_multiple(ml_fit = ml_fit, data.frame(gc_content = gc_content, 
        RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")), 
        RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")), exon_length = exon_length, 
        distance_to_junction = distance_to_junction, evolutionary_conservation = evolutionary_conservation, 
        DNA_5mer = DNA_5mer), positive_threshold = positive_threshold)
    returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
    names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
    return(returned_vector)
}

