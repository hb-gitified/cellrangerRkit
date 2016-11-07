#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

#' GeneBCMatrix class
#'

#' An S4 class to represent a matrix containing cell-barcoded gene counts
#' @slot barcode_filtered Whether the matrix contains only cell-containing barcodes (vs all sequenced barcodes)
#' @slot pipestance_path The path to the pipestance that this sample was loaded from
#' @slot summary Metrics summary
#' @slot subsampled The matrix was read-subsampled
#' @slot molecule_info Molecule info; used for read-subsampling
#' @import data.table
#' @import Biobase
GeneBCMatrix <- setClass("GeneBCMatrix",
                         contains="ExpressionSet",
                         slots=list(
                           barcode_filtered="logical",
                           pipestance_path="character",
                           summary="list",
                           subsampled="logical",
                           molecule_info="data.table"
                         ),
                         prototype = prototype( new("VersionedBiobase",
                                                    versions = c( classVersion("ExpressionSet"), GeneBCMatrix = "1.0.0" ) ))
)

#' Create a new GeneBCMatrix
#' @export
#' @param mat A sparse matrix whose rows are genes and columns are cell-barcodes
#' @param fd A dataframe with gene IDs as row names and additional metadata such as gene symbols as columns.
#' @param pd A dataframe with cell barcodes as row names and any additional metadata columns.
#' @param template A GeneBCMatrix to copy slots from
#' @return A new GeneBCMatrix object
#' @examples
#' \dontrun{
#' matrix = exprs(gbm1)
#' featureData = fData(gbm1)
#' phenoData = pData(gbm1)
#'
#' # Note including a template copies all other slots of gbm to new object
#' new_gbm = newGeneBCMatrix(matrix, fd=featureData, pd=phenoData, template=gbm1)
#'}
newGeneBCMatrix <- function(mat, fd, pd, template=NULL) {
  res <- new("GeneBCMatrix", assayData=assayDataNew( "environment", exprs=mat ), featureData=new("AnnotatedDataFrame", data=fd), phenoData=new("AnnotatedDataFrame", data=pd), molecule_info=data.table())

  if (!is.null(template)) {
    for (slot_name in setdiff(slotNames(res), c('assayData', 'featureData', 'phenoData'))) {
      slot(res, slot_name) <- slot(template, slot_name)
    }
  }

  res
}

#' dim for GeneBCMatrix objects
#' @export
#' @param x a GeneBCMatrix object
#' @return Result of dim function on exprs() on GeneBCMatrix
#' @examples
#' gene_count = dim(gbm1)[1]
#' cell_barcode_count = dim(gbm1)[2]
setMethod("dim", "GeneBCMatrix", function(x) dim(exprs(x)))

#' nrow for GeneBCMatrix objects
#' @export
#' @param x a GeneBCMatrix object
#' @return Result of nrow function on exprs() on GeneBCMatrix
#' @examples
#' gene_count = nrow(gbm1)
setMethod("nrow", "GeneBCMatrix", function(x) nrow(exprs(x)))

#' ncol for GeneBCMatrix objects
#' @export
#' @param x a GeneBCMatrix object
#' @return Result of ncol function on exprs() on GeneBCMatrix
#' @examples
#' cell_barcode_count = ncol(gbm1)
setMethod("ncol", "GeneBCMatrix", function(x) ncol(exprs(x)))

#' dimnames for GeneBCMatrix objects
#' @export
#' @param x a GeneBCMatrix object
#' @return Result of dimnames function on exprs() on GeneBCMatrix
#' @examples
#' gene_ids = dimnames(gbm1)[1]
#' cell_barcodes = dimnames(gbm1)[2]
setMethod("dimnames", "GeneBCMatrix", function(x) dimnames(exprs(x)))

#' rowSums for GeneBCMatrix objects
#' @export
#' @param x a GeneBCMatrix object
#' @param na.rm na.rm argument base function
#' @param dims dims argument passed to base function
#' @return Result of rowSums function on exprs() on GeneBCMatrix
#' @examples
#' umi_counts_per_gene = rowSums(gbm1)
setMethod("rowSums", "GeneBCMatrix", function(x, na.rm, dims=1) rowSums(exprs(x), na.rm, dims))

#' colSums for GeneBCMatrix objects
#' @export
#' @param x a GeneBCMatrix object
#' @param na.rm na.rm argument base function
#' @param dims dims argument passed to base function
#' @return Result of colSums function on exprs() on GeneBCMatrix
#' @examples
#' umi_counts_per_cell = colSums(gbm1)
setMethod("colSums", "GeneBCMatrix", function(x, na.rm, dims=1) colSums(exprs(x), na.rm, dims))

#' rowMeans for GeneBCMatrix objects
#' @export
#' @param x a GeneBCMatrix object
#' @param na.rm na.rm argument base function
#' @param dims dims argument passed to base function
#' @return Result of rowMeans function on exprs() on GeneBCMatrix
#' @examples
#' average_umis_per_gene = rowMeans(gbm1)
setMethod("rowMeans", "GeneBCMatrix", function(x, na.rm, dims=1) rowMeans(exprs(x), na.rm, dims))

#' colMeans for GeneBCMatrix objects
#' @export
#' @param x a GeneBCMatrix object
#' @param na.rm na.rm argument base function
#' @param dims dims argument passed to base function
#' @return Result of colMeans function on exprs() on GeneBCMatrix
#' @examples
#' average_umis_per_cell = colMeans(gbm1)
setMethod("colMeans", "GeneBCMatrix", function(x, na.rm, dims=1) colMeans(exprs(x), na.rm, dims))

######################################
# Non-exported functions for delegation
######################################
#' Delegate a binary operator to the matrix slot
#' @param e1 LHS
#' @param e2 RHS
#' @return result from operator
binary_op <- function(e1, e2) newGeneBCMatrix(mat=callGeneric(exprs(e1), exprs(e2)), fd=fData(e1), pd=pData(e1), e1)

#' Delegate a binary operator to the matrix slot
#' @param e1 LHS
#' @param e2 RHS
#' @return result from operator
binary_op_left <- function(e1, e2) newGeneBCMatrix(mat=callGeneric(exprs(e1), exprs(e2)), fd=fData(e1), pd=pData(e1), e1)

#' Delegate a binary operator to the matrix slot
#' @param e1 LHS
#' @param e2 RHS
#' @return result from operator
binary_op_right <- function(e1, e2) newGeneBCMatrix(mat=callGeneric(e1, exprs(e2)), fd=fData(e2), pd=pData(e2), e2)

#' Delegate a unary operator to the matrix slot
#' @param x GeneBCMatrix object
#' @return result from operator
unary_op <- function(x) newGeneBCMatrix(mat=callGeneric(exprs(x)), fd=fData(x), pData(x), x)

#' Delegate a unary arithmetic operator to the matrix slot
#' @param e1 GeneBCMatrix object
#' @param e2 dummy
#' @return result from operator
unary_arith_op <- function(e1, e2) newGeneBCMatrix(mat=callGeneric(exprs(e1)), fd=fData(e1), pd=pData(e1), e1)

######################################
# Non-exported functions for delegation
######################################

#' Unary negation of GeneBCMatrix objects
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 dummy
#' @return result from operator
setMethod("-", signature(e1 = "GeneBCMatrix", e2 = "missing"), unary_arith_op)

#' Unary identity of GeneBCMatrix objects
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 dummy
#' @return result from operator
setMethod("+", signature(e1 = "GeneBCMatrix", e2 = "missing"), unary_arith_op)

#' Arithmetic operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Arith", signature(e1 = "GeneBCMatrix", e2 = "GeneBCMatrix"), binary_op)

#' Arithmetic operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 logical
#' @return result from operator
setMethod("Arith", signature(e1 = "GeneBCMatrix", e2 = "logical"), binary_op_left)

#' Arithmetic operator for GeneBCMatrix
#' @export
#' @param e1 logical
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Arith", signature(e1 = "logical", e2 = "GeneBCMatrix"), binary_op_right)

#' Arithmetic operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 numeric
#' @return result from operator
setMethod("Arith", signature(e1 = "GeneBCMatrix", e2 = "numeric"), binary_op_left)

#' Arithmetic operator for GeneBCMatrix
#' @export
#' @param e1 numeric
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Arith", signature(e1 = "numeric", e2 = "GeneBCMatrix"), binary_op_right)


#' Compare operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Compare", signature(e1 = "GeneBCMatrix", e2 = "GeneBCMatrix"), binary_op)


#' Compare operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Compare", signature(e1 = "GeneBCMatrix", e2 = "numeric"), binary_op_left)

#' Compare operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Compare", signature(e1 = "numeric", e2 = "GeneBCMatrix"), binary_op_right)

#' Compare operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Compare", signature(e1 = "logical", e2 = "GeneBCMatrix"), binary_op_left)


#' Compare operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Compare", signature(e1 = "GeneBCMatrix", e2 = "logical"), binary_op_left)

#' Logical operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Logic", signature(e1 = "GeneBCMatrix", e2 = "GeneBCMatrix"), binary_op)

#' Logical operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Logic", signature(e1 = "GeneBCMatrix", e2 = "numeric"), binary_op_left)

#' Logical operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Logic", signature(e1 = "numeric", e2 = "GeneBCMatrix"), binary_op_right)

#' Logical operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Logic", signature(e1 = "GeneBCMatrix", e2 = "logical"), binary_op_left)

#' Logical operator for GeneBCMatrix
#' @export
#' @param e1 GeneBCMatrix
#' @param e2 GeneBCMatrix
#' @return result from operator
setMethod("Logic", signature(e1 = "logical", e2 = "GeneBCMatrix"), binary_op_right)

#' Math operator for GeneBCMatrix
#' @export
#' @param x GeneBCMatrix
#' @return result from operator
setMethod("Math", signature(x = "GeneBCMatrix"), unary_op)

#' Math2 operator for GeneBCMatrix
#' @export
#' @param x GeneBCMatrix
#' @return result from operator
setMethod("Math2", signature(x = "GeneBCMatrix"), unary_op)

#' Summary operator for GeneBCMatrix
#' @export
#' @param x GeneBCMatrix
#' @return result from operator
setMethod("Summary", signature(x = "GeneBCMatrix"), function(x) callGeneric(exprs(x)))
