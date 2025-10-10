cache <- new.env()
cache$rna <- list()

#' Get a test scRNA-seq dataset
#'
#' Get a test scRNA-seq or CITE-seq dataset with varying levels of processing.
#' This uses caching to avoid recomputation.
#'
#' @param at String specifying the level of processing.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} containing a dataset at the specified level of processing.
#' 
#' For \code{getTestRnaData}, this is a scRNA-seq dataset of the mouse brain,
#' where the main experiment contains RNA counts and the alternative experiments contain ERCC and repeat element counts.
#'
#' For \code{getTestAdtData}, this is a CITE-seq dataset of human PBMCs,
#' where the main experiment contains RNA counts and the alternative experiment contains ADT counts.
#'
#' @author Aaron Lun
#' @examples
#' getTestRnaData.se()
#' getTestAdtData.se()
#'
#' @seealso
#' \code{\link[scRNAseq]{ZeiselBrainData}}, for the original data source for \code{getTestRnaData.se}.
#'
#' \code{\link[scRNAseq]{KotliarovPBMCData}}, for the original data source for \code{getTestAdtData.se}.
#' 
#' @export
#' @name getTestData.se
getTestRnaData.se <- function(at = "start") {
    at <- match.arg(at)

    if (!("start" %in% names(cache$rna))) {
        cache$rna$start <- scRNAseq::ZeiselBrainData()
    }
    if (at == "start") {
        return(cache$rna[[at]])
    }
}

#' @export
#' @rdname getTestData.se
getTestAdtData.se <- function(at = "start") {
    at <- match.arg(at)

    if (!("start" %in% names(cache$adt))) {
        cache$adt$start <- scRNAseq::KotliarovPBMCData()
    }
    if (at == "start") {
        return(cache$adt[[at]])
    }
}
