#' Convert output from vendor software to a format suitable for analysis with
#' HDXstats
#'
#' @param dataset The HDXMS dataset from an external software tool. Please
#'   make sure that for each peptide the same number of replicates are
#'   reported. If not the global threshold, hybrid significance test and the
#'   Ttest produce unreliable results
#' @param type Which software tool has been used to extract information from
#'   the spectra. The default value is DynamX cluster.
#'
#' @return A data frame with columns: peptide, start, end, label_time, state and
#'   deuterium_level
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @export
calculateDeuteriumLevels <- function(dataset, type = "DynamX cluster") {

    Sequence <- NULL
    if (type == "DynamX cluster") {
      sequences <- unique(dataset$Sequence)
      nbSeq <- length(sequences)
      dataset$State <- as.factor(dataset$State)
      stateData <- levels(dataset$State)
      positions <- unique(paste(dataset$Start, dataset$End, sep="_"))
      positions <- matrix(unlist(strsplit(positions,"_")),ncol=2, byrow=TRUE)
      df <- list()

      for (k in 1:nbSeq) {
        retained <- subset(dataset, Sequence == sequences[k])

        tmp <- list()

        for (l in 1:length(stateData)) {
          state <- stateData[l]
          idxState <- which(retained$State == state)

          if (length(idxState) != 0) {
            retainedCenter <- retained[idxState,]

            idx <- which(retainedCenter$Exposure == 0.0)
            if(length(idx)==0){
              stop("Labeling time points do not contain 0.0")
            }else{
              centerMass <- mean(retainedCenter$Center[idx] * retainedCenter$z[idx] - retainedCenter$z[idx] *              1.0079)

              deuteration <- (retainedCenter$Center * retainedCenter$z - retainedCenter$z * 1.0079) - centerMass
              exposure_time <- retainedCenter$Exposure
              tmp[[l]] <- data.frame(
                "peptide" = sequences[k],
                "start" = as.numeric(positions[k,1]),
                "end" = as.numeric(positions[k,2]),
                "label_time" = exposure_time,
                "state" = state,
                "deuterium_level" = deuteration
                )
            }
          }
        }

        df[[k]] <- do.call(rbind, tmp)
      }

      output <- do.call(rbind, df)
    }

    output
}

#' Convert the processed HDXMS dataset in long format to the wide format which is required for the moderated test-statistics.
#'
#' @param datasetLong The HDXMS dataset in long format with following columns: peptide, label_time, state and deuterium_level
#' @param nbRep The number of replicates. Make sure that each peptide has exactly the same number of replicates as specified here.
#'
#' @return The HDXMS dataset in wide format
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @export
convertLongToWide <- function(datasetLong, nbRep){

  datasetLong <- datasetLong[,c("peptide","label_time","state","deuterium_level")]
  nbTime <- length(unique(datasetLong$label_time))
  datasetLong$replicates <- rep(1:nbRep,dim(datasetLong)[1]/nbRep)
  datasetLong$id <- paste(paste0("state_",datasetLong$state), paste0("label_time_", datasetLong$label_time), paste0("Replicate", datasetLong$replicates), sep="_")
  df <- reshape(datasetLong, idvar="peptide", timevar="id", drop=c("label_time","state","replicates"), direction="wide")
  colnames(df) <- gsub("^.*?\\.","",colnames(df))
  colnames(df) <- gsub("_Replicate.","",colnames(df))
  df
}

#' Create the contrast matrix for the moderated test statistics
#'
#' @param states character vector specifying all protein states of interest
#' @param label_time character vector specifying all time points of interest
#' @param designMatrix matrix that corresponds to the design of the HDX-MS experiment
#'
#' @return Matrix which columns correspond to contrasts
#' @examples
#' #Please check the vignette of the HDXstats package.
#' @importFrom limma makeContrasts
#' @importFrom utils combn
#' @export
createContrastMatrix <- function(states, label_time, designMatrix){

  statePairwiseComparisons <- t(combn(states, 2))
  contrastV <- unlist(lapply(label_time, function(iTime) {
    apply(statePairwiseComparisons, 1, function(iComparison) {
      paste0(iComparison[1], ".", iTime, "-", iComparison[2], ".", iTime)
    })
  }))

  output <- makeContrasts(contrasts=contrastV, levels=designMatrix)
}
