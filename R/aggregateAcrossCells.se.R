#' Aggregate expression values across cells
#' 
#' Aggregate expression values across groups of cells for each gene.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param factors,num.threads Arguments to pass to \code{\link[scrapper]{aggregateAcrossCells}}.
#' @param more.aggr.args Named list of additional arguments to pass to \code{\link[scrapper]{aggregateAcrossCells}}.
#' @param assay.type Integer or string specifying the assay of \code{x} to be aggregated.
#' @param output.prefix String containing a prefix to add to the names of the \code{link[SummarizedExperiment]{colData}} columns containing the factor combinations.
#' @param counts.name String containing the name of the \code{\link[SummarizedExperiment]{colData}} column containing the cell count for each factor combination.
#' If \code{NULL}, the cell counts are not reported.
#' @param index.name String containing the name of the \code{\link[SummarizedExperiment]{metadata}} entry containing the combination index for each cell.
#' If \code{NULL}, the indices are not reported.
#' @param include.coldata Logical scalar indicating whether to add the aggregated \code{colData} from \code{x} to the output.
#' @param more.coldata.args Named list of additional arguments to pass to \code{aggregateColData}.
#' Only relevant if \code{include.coldata=TRUE}.
#' @param altexps Integer or character vector of alternative experiments of \code{x} to aggregate.
#' The aggregated assay from each alternative experiment is determined by \code{assay.type}.
#' Only relevant if \code{x} is a \link[SingleCellExperiment]{SingleCellExperiment}.
#' @param copy.altexps Logical scalar indicating whether to copy the \code{colData} and \code{metadata} of the output SingleCellExperiment into each of its alternative experiments.
#' @param coldata \link[S4Vectors]{DataFrame} of column data, containing one row for each cell.
#' @param index Integer vector containing the index of the factor combination to which each cell in \code{coldata} was assigned.
#' @param number Integer specifying the total number of unique factor combinations.
#' @param only.atomic Logical scalar specifying whether to skip non-atomic, non-factor columns.
#' @param placeholder Placeholder value to store in the output column when a factor combination does not have a single unique value. 
#'
#' @details
#' \code{factors} is expected to be an ordinary list of categorical variables.
#' However, it may also be a \link[S4Vectors]{List} or \link[S4Vectors]{DataFrame}, which will be coerced to a corresponding list.
#' 
#' Alternatively, \code{factors} may be an atomic vector or factor representing a single variable.
#' In such cases, \code{output.prefix} is used as the name.
#' 
#' @return
#' For \code{aggregateAcrossCells.se}, a SummarizedExperiment is returned where each column corresponds to a factor combination.
#' Each row corresponds to a gene in \code{x} and the \code{\link[SummarizedExperiment]{rowData}} is taken from \code{x}.
#' The assays contain the sum of counts (\code{"sums"}) and the number of detected cells (\code{"detected"}) in each combination for each gene.
#' The \code{colData} contains:
#' \itemize{
#' \item The factor combinations, with column names prefixed by \code{output.prefix}.
#' \item The cell count for each combination, named by \code{counts.name}.
#' \item Additional \code{colData} from \code{x} if \code{include.coldata=TRUE}.
#' This is aggregated with \code{aggregateColData} on the combination indices.
#' }
#' The metadata contains a vector of length equal to the number of cells in \code{x},
#' where each value is an index of the factor combination (i.e., column of the output object) to which that cell was assigned.
#'
#' If \code{altexps} is specified, a SingleCellExperiment is returned instead.
#' The same aggregation for the main experiment is applied to each alternative experiment.
#' If \code{copy.altexps=TRUE}, the \code{colData} for each alternative experiment will contain a copy of the factor combinations and cell counts,
#' and the \code{metadata} will contain a copy of the index vector.
#'
#' For \code{aggregateColData}, a \link[S4Vectors]{DFrame} is returned with number of rows equal to \code{number}.
#' Each atomic or factor column in \code{coldata} is represented by a column in the output DFrame.
#' In each column, the \code{j}-th entry is equal to the unique value of all rows in \code{coldata} with \code{index} equal to \code{j},
#' or \code{placeholder} if there is not exactly one unique value.
#' If \code{only.atomic=FALSE}, any non-atomic/non-factor columns of \code{coldata} are represented in the output DFrame by a vector of \code{placeholder} values.
#' If \code{only.atomic=TRUE}, any non-atomic/non-factor columns of \code{coldata} are skipped.
#' 
#' @author Aaron Lun
#' @examples
#' sce <- getTestRnaData.se("start")
#' aggr <- aggregateAcrossCells.se(sce, sce$level1class)
#' head(assay(aggr))
#' colData(aggr)
#'
#' @export
#' @importFrom methods is
#' @importClassesFrom S4Vectors List
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData colData<- metadata<-
#' @importFrom SingleCellExperiment SingleCellExperiment altExpNames altExp<- altExp mainExpName
aggregateAcrossCells.se <- function(
    x,
    factors,
    num.threads = 1,
    more.aggr.args = list(),
    assay.type = "counts",
    output.prefix = "factors",
    counts.name = "counts",
    index.name = "index",
    include.coldata = TRUE,
    more.coldata.args = list(),
    altexps = NULL,
    copy.altexps = FALSE
) {
    if (is.list(factors)) {
        # this is fine.
    } else if (is(factors, "List")) {
        factors <- as.list(factors)
    } else {
        factors <- list(factors)
        names(factors) <- output.prefix
        output.prefix <- NULL
    }

    out <- do.call(
        scrapper::aggregateAcrossCells,
        c(
            list(assay(x, assay.type)),
            .collapse_args(
                list(factors=factors, num.threads=num.threads),
                more.aggr.args
            )
        )
    )

    CON <- SummarizedExperiment
    if (length(altexps)) {
        CON <- SingleCellExperiment
    }
    se <- CON(out[c("sums", "detected")], rowData=rowData(x))

    common.cd <- DataFrame(out$combinations, check.names=FALSE)
    colnames(common.cd) <- paste0(output.prefix, colnames(common.cd))
    if (!is.null(counts.name)) {
        common.cd[[counts.name]] <- out$counts
    }
    colData(se) <- common.cd
    if (include.coldata) {
        aggr.cd <- do.call(aggregateColData, c(list(colData(x), out$index, number=nrow(out$combinations)), more.coldata.args))
        colData(se) <- cbind(colData(se), aggr.cd)
    }

    if (!is.null(index.name)) {
        metadata(se)[[index.name]] <- out$index
    }

    if (length(altexps)) {
        mainExpName(se) <- mainExpName(x)
        ae.names <- altExpNames(x)
        if (is.numeric(altexps)) {
            altexps <- ae.names[altexps]
        }

        for (ae in altexps) {
            ae.se <- Recall(
                altExp(se,ae),
                out$index,
                num.threads=num.threads,
                more.aggr.args=more.aggr.args,
                assay.type=assay.type,
                altexps=NULL,
                output.prefix=NULL,
                counts.name=NULL,
                include.coldata=include.coldata
            )

            ae.cd <- colData(ae.se)[,-1,drop=FALSE] # remove uninteresting factor combination
            if (copy.altexps) {
                ae.cd <- cbind(common.cd, ae.cd)
            }
            colData(ae.se) <- ae.cd

            if (copy.altexps) {
                metadata(ae.se) <- metadata(se)
            }
            altExp(se, ae) <- ae.se
        }
    }

    se
}

#' @export
#' @rdname aggregateAcrossCells.se
#' @importFrom S4Vectors make_zero_col_DFrame
aggregateColData <- function(coldata, index, number, only.atomic = TRUE, placeholder = NA) {
    collected <- make_zero_col_DFrame(nrow=number)
    index <- factor(index, seq_len(number))

    for (cn in colnames(coldata)) {
        curcol <- coldata[[cn]]
        if (!is.atomic(curcol) && !is.factor(curcol)) {
            if (!only.atomic) {
                collected[[cn]] <- rep(placeholder, number)
            }
            next
        }

        grouped <- split(curcol, index)
        alloc <- rep(curcol[1], number)
        for (i in seq_along(grouped)) {
            u <- unique(grouped[[i]])
            if (length(u) != 1L) {
                alloc[i] <- placeholder
            } else {
                alloc[i] <- u
            }
        }

        collected[[cn]] <- alloc
    }

    collected
}
