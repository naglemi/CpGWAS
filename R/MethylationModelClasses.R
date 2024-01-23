#' MethylationBase Class
#'
#' This class represents an elastic net model for a single base position,
#' explaining methylation as a function of a SNP matrix.
#' It also contains information about methylation positions, associated models,
#' and statistical metrics.
#'
#' @slot methylationPosition numeric Base position for methylation site.
#' @slot windowSize numeric Window size for SNPs surrounding methylation site.
#' @slot n_SNPs numeric Number of SNPs in the window surrounding methylation site.
#' @slot glmnetModel list List of `glmnet` model objects.
#' @slot snpWeights list List of SNP weight vectors (obtained from `glmnet`).
#' @slot cor numeric Correlation coefficient between predicted, observed values.
#' @slot mse numeric Mean squared error of prediction for methylation.
#' @slot alpha numeric Alpha parameter from `glmnet` model (tuned with `glmnet_tune_alpha`).
#' @slot lambda numeric Lambda parameter from `glmnet` model (tuned with `cv.glmnet`).
#' @slot runtime numeric Computation time for the model (including tuning, cross-validation).
#' @export
#'
setClass(
  "MethylationBase",
  slots = c(
    methylationPosition = "numeric",
    windowSize = "numeric",
    n_SNPs = "numeric",
    glmnetModel = "ANY",
    snpWeights = "ANY",
    cor = "numeric",
    mse = "numeric",
    alpha = "numeric",
    lambda = "numeric",
    runtime = "numeric"
  )
)

setValidity("MethylationBase", function(object) {
  if (!inherits(object@glmnetModel, c("elnet", "glmnet"))) {
    return("glmnetModel must contain a glmnet model object")
  }
  TRUE
})


setValidity("MethylationBase", function(object) {
  if (!inherits(object@snpWeights, c("dgCMatrix"))) {
    return("snpWeights must contain a dgCMatrix object")
  }
  TRUE
})


#' MethylationScaff Class
#'
#' Represents a scaffold structure in methylation analysis, containing multiple
#'  `MethylationBase` models.
#' Each `MethylationScaff` is identified by a unique scaffold identifier.
#'
#' @slot scaffoldIdentifier character Unique identifier for the scaffold.
#' @slot models list List of MethylationBase objects.
#' @export
setClass(
  "MethylationScaff",
  slots = c(
    scaffoldIdentifier = "character",
    models = "list"
  )
)

#' Constructor for MethylationScaff Class
#'
#' Initializes a `MethylationScaff` object with a given identifier and a list of
#' `MethylationBase` models.
#' This is intended to store many `MethylationBase` models for the same scaffold,
#' typically a whole chromosome.
#'
#' @param scaffoldIdentifier character The identifier for the scaffold.
#' @param models list List of MethylationBase objects.
#' @return MethylationScaff object.
#'
#' @examples
#' \dontrun{
#' modelsList <- list(new("MethylationBase", ...), new("MethylationBase", ...))
#' scaff <- new("MethylationScaff", scaffoldIdentifier = "scaff1", models = modelsList)
#' }
#'

setMethod(
  "initialize",
  "MethylationScaff",
  function(.Object, scaffoldIdentifier, models) {
    if (!is.character(scaffoldIdentifier)) {
      stop("scaffoldIdentifier must be a character string.")
    }
    .Object@scaffoldIdentifier <- scaffoldIdentifier
    .Object@models <- lapply(models, function(model) {
      if (!inherits(model, "MethylationBase")) {
        stop("models must be of class MethylationBase.")
      }
      model
    })
    .Object
  }
)

#' Save MethylationScaff Object
#'
#' Saves a `MethylationScaff` object to an RDS file.
#' The file is named after the scaffold identifier and saved in the specified output directory.
#'
#' @param object MethylationScaff The MethylationScaff object to save.
#' @param outputDir character The directory to save the file, default is "output/".
#' @export
setGeneric("saveMethylationScaff", function(object, outputDir="output/") {
  standardGeneric("saveMethylationScaff")
})

#' Method to Save MethylationScaff Object
#'
#' @param object MethylationScaff The MethylationScaff object to save.
#' @param outputDir character The directory to save the file, default is "output/".
#' @rdname saveMethylationScaff
#' @export
setMethod("saveMethylationScaff", "MethylationScaff", function(object, outputDir="output/") {
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }

  fileName <- paste0(object@scaffoldIdentifier, ".rds")
  filePath <- file.path(outputDir, fileName)

  saveRDS(object, file = filePath)

  message(paste("MethylationScaff object saved to", filePath))
})

#' MethylationGenome Class
#'
#' Represents the genome-level structure in methylation analysis.
#' Contains multiple `MethylationScaff` objects, each associated with different scaffolds.
#'
#' @slot scaffoldIdentifier list List of scaffold identifiers.
#' @slot MethylationScaffs list List of MethylationScaff objects.
setClass(
  "MethylationGenome",
  slots = c(
    scaffoldIdentifier = "list",
    MethylationScaffs = "list"
  )
)

#' Constructor for MethylationGenome Class
#'
#' Initializes a `MethylationGenome` object with a list of scaffold identifiers and
#' corresponding `MethylationScaff` objects.
#'
#' @param scaffoldIdentifiers list List of scaffold identifiers.
#' @param MethylationScaffs list List of `MethylationScaff` objects.
#' @return MethylationGenome object.
#' @examples
#' \dontrun{
#' genome <- new("MethylationGenome", scaffoldIdentifiers = list("scaff1", "scaff2"),
#' MethylationScaffs = list(scaff1, scaff2))
#' }
setMethod(
  "initialize",
  "MethylationGenome",
  function(.Object, scaffoldIdentifiers, MethylationScaffs) {
    if (!is.list(scaffoldIdentifiers)) {
      stop("scaffoldIdentifiers must be a list.")
    }
    .Object@scaffoldIdentifiers <- scaffoldIdentifiers
    .Object@MethylationScaffs <- lapply(MethylationScaffs, function(scaff) {
      if (!inherits(scaff, "MethylationScaff")) {
        stop("MethylationScaffs must be of class MethylationScaff.")
      }
      scaff
    })
    .Object
  }
)

#' Convert MethylationModels to Data Frame
#'
#' Creates a data frame representation of `MethylationScaff`.
#' Each row corresponds to a `MethylationScaff`, with columns for each non-list attribute.
#'
#' @param object MethylationModels The `MethylationScaff` object to convert.
#' @return data.frame Data frame representation of `MethylationScaff`,
#' with column for each attribute except models and SNP weights.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df <- convertToDataFrame(MethylationScaff)
#' }
convertToDataFrame <- function(object) {
  if (!inherits(object, "MethylationScaff")) {
    stop("The object must be of class 'MethylationScaff'.")
  }

  modelsList <- lapply(object@models, function(model) {
    data.frame(
      scaffoldIdentifier = object@scaffoldIdentifier,  # Scaffold identifier is the same for all models
      methylationPosition = model@methylationPosition,
      windowSize = model@windowSize,
      nSNPs = model@n_SNPs,
      cor = model@cor,
      mse = model@mse,
      alpha = model@alpha,
      lambda = model@lambda,
      runtime = model@runtime
    )
  })

  do.call("rbind", modelsList)
}

