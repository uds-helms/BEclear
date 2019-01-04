
#' preprocessBEclear
#'
#' @param data any matrix filled with beta values, column names have to be
#' sample_ids corresponding to the ids listed in "samples", row names have to
#' be gene names.
#' @param samples data frame with two columns, the first column has to contain
#' the sample numbers, the second column has to contain the corresponding batch
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' 
#' @description this methods does some preprocessing steps for the later methods
#' like removing rows containing only missing values
#' 
#' @details 
#' Here we describe the preprocessing steps in the order they are executed:
#' \itemize{
#' \item{Values below 0 or above 1 are set to \code{NA}, as the other methods 
#' expect methylation beta values}
#' \item{columns that only contain \code{NA}s are removed}
#' \item{rows that only contain \code{NA}s are removed}
#' \item{samples that are present in the data, but are not annoted in the samples 
#' are removed. If this is the case with your data-set, please check those samples.}
#' \item{samples that are annoted but not in the data matrix are removed}
#' \item{if there are duplicate sample names in the data matrix, all sample names
#' get replaced through a new unique ID. In this case a \code{\link[data.table]{data.table}}
#' containing the mapping is returned as well}
#' }
#'
#' @return a list containing the processed data and samples and a 
#' \code{\link[data.table]{data.table}} containing mappings from the original
#' sample names to the new ones. If sample names weren't changed this third
#' object is \code{NULL}
#' @export
#'
#' @examples
#' data(BEclearData)
#' res <- preprocessBEclear(ex.data, ex.samples)
preprocessBEclear <- function(data, samples) {
    ## checking if they're are values above 1 or below 0
    if (any(data > 1 | data < 0, na.rm = TRUE)) {
        flog.warn(paste(
            sum(data > 1 | data < 0, na.rm = TRUE),
            "values are above 1 or below 0. Check your data"
        ))
        flog.warn("Replacing them with missing values")
        data[data > 1 | data < 0] <- NA
    }
    
    ## checking if there are columns containing only missing values
    naIndices <- apply(data, 2, function(x) all(is.na(x)))
    if (any(naIndices, na.rm = TRUE)) {
        flog.warn("There are columns, that contain only missing values")
        flog.warn(paste(sum(naIndices), "columns get dropped"))
        data <- data[, !naIndices]
    }
    
    ## checking if there are rows containing only missing values
    naIndices <- apply(data, 1, function(x) all(is.na(x)))
    if (any(naIndices, na.rm = TRUE)) {
        flog.warn("There are rows, that contain only missing values")
        flog.warn(paste(sum(naIndices), "rows get dropped"))
        data <- data[!naIndices, ]
    }
    
    ## checking if there are samples that are not present in the samples matrix
    if (any(!colnames(data) %in% samples$sample_id)) {
        ids <- paste(colnames(data)[!colnames(data) %in% samples$sample_id],
                     collapse = ", "
        )
        flog.warn(paste(
            "The following samples are in the data, but not annotated",
            "in the samples matrix:", ids
        ))
        flog.warn("Dropping those samples now")
        data <- data[, colnames(data) %in% samples$sample_id]
    }
    
    if (any(!samples$sample_id %in% colnames(data))) {
        ids <- paste(samples$sample_id[!samples$sample_id %in% colnames(data)],
                     collapse = ", "
        )
        flog.warn(
            "The following samples are annotated in the sample matrix,",
            "but aren't contained in data matrix:", ids
        )
        flog.warn("Dropping those samples now")
        samples <- samples[sample_id %in% colnames(data)]
    }
    
    samples <- data.table(samples)
    uniqueIDsToSamples <- NULL
    
    ## checking if there are duplicated sample names
    if (any(duplicated(colnames(data)))) {
        flog.warn("Sample names aren't unique")
        flog.warn(paste(
            "Transforming them to unique IDs. List with annotations will",
            "be added to the results"
        ))
        uniqueIDsToSamples <- data.table(
            sample_id = colnames(data),
            unique_id =
                as.character(seq_along(colnames(data)))
        )
        colnames(data) <- uniqueIDsToSamples$unique_id
        samples <- samples[uniqueIDsToSamples, ,
                           on = .(sample_id = sample_id)
                           ][
                               ,
                               .(
                                   sample_id = unique_id,
                                   batch_id = batch_id
                               )
                               ]
    }
    setkey(samples, "batch_id", "sample_id")
    
    return(list(data=data, samples = samples, 
                uniqueIDsToSamples = uniqueIDsToSamples ))
}