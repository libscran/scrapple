cache <- new.env()
cache$rna <- list()
cache$adt <- list()
cache$crispr <- list()

#' Get a test scRNA-seq dataset
#'
#' Get a test scRNA-seq or CITE-seq dataset with varying levels of processing.
#' This uses caching to avoid recomputation.
#'
#' @param at String specifying the level of processing.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} containing a dataset at the specified level of processing.
#' 
#' @details
#' For \code{getTestRnaData}, this is a scRNA-seq dataset of the mouse brain,
#' where the main experiment contains RNA counts and the alternative experiments contain ERCC and repeat element counts.
#' This is obtained with \code{fetchDataset("zeisel-brain-2015", "2023-12-14")}.
#'
#' For \code{getTestAdtData}, this is a CITE-seq dataset of human PBMCs,
#' where the main experiment contains RNA counts and the alternative experiment contains ADT counts.
#' This is obtained with \code{fetchDataset("kotliarov-pbmc-2020", "2024-04-18")}.
#'
#' For \code{getTestCrisprData}, this is a Perturb-seq dataset of a pancreatic beta cell line,
#' where the main experiment contains RNA counts and the alternative experiment contains CRISPR guide counts.
#' This is obtained with \code{fetchDataset("cao-pancreas-2025", "2025-10-10", "rqc")}.
#'
#' @author Aaron Lun
#' @examples
#' getTestRnaData.se()
#' getTestAdtData.se()
#' getTestCrisprData.se()
#'
#' @seealso
#' \code{\link[scRNAseq]{fetchDataset}}, used to obtain each dataset.
#' 
#' @export
#' @name getTestData.se
getTestRnaData.se <- function(at = c("start", "qc")) {
    at <- match.arg(at)

    if (!("start" %in% names(cache$rna))) {
        cache$rna$start <- scRNAseq::fetchDataset("zeisel-brain-2015", "2023-12-14", realize.assays=TRUE)
    }
    sce <- cache$rna$start
    if (at == "start") {
        return(sce)
    }

    if (!("qc" %in% names(cache$rna))) {
        sce <- quickRnaQc.se(sce, subsets=list(mito=startsWith(rownames(sce), "mt-")), altexp.proportions=c(ercc="ERCC"))
        sce <- sce[,sce$keep]
        cache$rna$qc <- sce
    }
    if (at == "qc") {
        return(sce)
    }
}

#' @export
#' @rdname getTestData.se
getTestAdtData.se <- function(at = c("start", "qc")) {
    at <- match.arg(at)

    if (!("start" %in% names(cache$adt))) {
        cache$adt$start <- scRNAseq::fetchDataset("kotliarov-pbmc-2020", "2024-04-18", realize.assays=TRUE)
    }
    sce <- cache$adt$start
    if (at == "start") {
        return(sce)
    }

    if (!("qc" %in% names(cache$rna))) {
        sce <- quickRnaQc.se(sce, subsets=list(mito=startsWith(rownames(sce), "MT-")))
        altExp(sce, "ADT") <- quickAdtQc.se(altExp(sce, "ADT"), subsets=list(igg=rowData(altExp(sce, "ADT"))$isotype))
        sce <- sce[,sce$keep & altExp(sce, "ADT")$keep]
        cache$rna$qc <- sce
    }
    sce <- cache$rna$qc
    if (at == "qc") {
        return(sce)
    }
}

#' @export
#' @rdname getTestData.se
getTestCrisprData.se <- function(at = "start") {
    at <- match.arg(at)

    if (!("start" %in% names(cache$crispr))) {
        cache$crispr$start <- scRNAseq::fetchDataset("cao-pancreas-2025", "2025-10-10", "rqc", realize.assays=TRUE)
    }
    if (at == "start") {
        return(cache$crispr[[at]])
    }
}
