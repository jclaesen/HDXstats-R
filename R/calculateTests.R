#' Calculate the moderated T-test for differential HDX-MS data
#'
#' @param design the design matrix of the differential HDX-MS experiment, with rows corresponding to the observations and columns to coefficients to be estimated
#' @param contrastMatrix contrast matrix of the experiment contrasts between a set of parameters as a numeric matrix. The matrix specifies which pairwise comparisons between the coefficients are to be extracted from the fit
#' @param mod.test which moderated test to use: the moderated t-test ('t.test') or the thresholded moderated t-test ('thresholded')
#' @param threshold required minimum difference in deuteration-levels to be tested by the thresholded moderated t-test
#' @param adjust.method the multiple testing correction method
#' @param dataset the HDXMS dataset in a long format with following columns: peptide, label_time, state and deuterium_level
#'
#' @details
#'
#' The design matrix can be constructed with the model.matrix function of the stats-package, while the function makeContrasts() of the limma-package can be used to create the contrast matrix.
#'
#'
#'
#' @return A list with the results from the moderated T-tests. Each element in the list corresponds to a specific comparison as defined by the contrast matrix.
#' @examples
#'
#' #Please check the vignette of the HDXstats package.
#'
#' @importFrom limma lmFit topTable contrasts.fit eBayes treat
#' @export
moderated.Ttest <- function(design = NULL, contrastMatrix = NULL, mod.test = NULL, threshold = NULL, adjust.method = "BH", dataset = NULL) {

    if (is.null(design))
      stop("No design matrix specified")

    if (is.null(contrastMatrix))
      stop("No contrast matrix specified")

    if (is.null(dataset))
      stop("No dataset provided")

    if (!is.null(mod.test)) {
      mod.test <- tolower(mod.test)
      fittedModel <- lmFit(dataset, design)
      contrastModel <-
        contrasts.fit(fittedModel, contrastMatrix)
      comparisons <- colnames(contrastMatrix)

      res <- list()

      if (mod.test == "t.test") {
        fit <- eBayes(contrastModel)

        for (k in comparisons) {
          res[[k]] <-
            topTable(
              fit,
              coef = k,
              number = Inf,
              sort.by = "P",
              adjust.method = adjust.method
            )[,-6]
        }

      } else if (mod.test == "thresholded") {
        if (!is.null(threshold)) {
          fit2 <- treat(contrastModel, lfc = threshold)
          for (k in comparisons) {
            res[[k]] <-
              topTable(
                fit2,
                coef = k,
                number = Inf,
                sort.by = "P",
                adjust.method = adjust.method
              )[,-6]
          }

        } else{
          stop("No threshold specified")
        }

      } else{
        stop(paste(
          mod.test,
          "is not supported. Please specify a correct moderated test statistic.",
          sep = " "
        ))
      }

    } else{
      stop("No moderated test statistic specified.")
    }
  res <- lapply(res, function(x) {
    colnames(x)[colnames(x)=="logFC"] <- "Deuteration_Difference"
    colnames(x)[colnames(x)=="AveExpr"] <- "Average_Deuteration"
    return(x)
  })
    return(res)
}

#' Calculate the moderated F-test for differential HDX-MS data
#'
#' @param design the design matrix of the differential HDX-MS experiment, with rows corresponding to the observations and columns to coefficients to be estimated
#' @param contrastMatrix contrast matrix of the experiment contrasts between a set of parameters as a numeric matrix. The matrix specifies which pairwise comparisons between the coefficients are to be extracted from the fit
#' @param comparisons a vector or data.frame or matrix that specifies which contrasts have to be tested with the moderated F-statistic
#' @param adjust.method the multiple testing correction method
#' @param dataset the HDXMS dataset in a long format with following columns: peptide, label_time, state and deuterium_level
#'
#' @details
#'
#' The design matrix can be constructed with the model.matrix function of the stats-package, while the function makeContrasts() of the limma-package can be used to create the contrast matrix.
#'
#' The dataset
#'
#' @return A list with the results from the moderated F-tests. Each element in the list corresponds to a specific comparison as defined by the parameter comparisons.
#' @examples
#'
#' #Please check the vignette of the HDXstats package.
#'
#' @importFrom limma lmFit topTable contrasts.fit eBayes
#' @export
moderated.Ftest <- function(design = NULL, contrastMatrix = NULL, comparisons = NULL, adjust.method = "BH",     dataset = NULL) {

    if (is.null(design))
      stop("No design matrix specified")

    if (is.null(contrastMatrix))
      stop("No contrast matrix specified")

    if (is.null(dataset))
      stop("No dataset provided")

    if (is.null(comparisons))
      stop("No comparisons provided")

    fittedModel <- lmFit(dataset, design)
    contrastModel <- contrasts.fit(fittedModel, contrastMatrix)
    fit <- eBayes(contrastModel)

    if (is.matrix(comparisons) | is.data.frame(comparisons)) {
      res <- list()
      nbComparisons <- dim(comparisons)[1]

      for (k in 1:nbComparisons) {
        res[[k]] <-
          topTable(
            fit,
            coef = comparisons[k, ] ,
            number = Inf,
            sort.by = "F",
            adjust.method = adjust.method
          )
      }
    } else if (class(comparisons) == "character") {
      res <-
        topTable(
          fit,
          coef = comparisons,
          number = Inf,
          sort.by = "F",
          adjust.method = adjust.method
        )
    }
    colnames(res)[colnames(res)=="logFC"] <- "Deuteration_Difference"
    colnames(res)[colnames(res)=="AveExpr"] <- "Average_Deuteration"
    return(res)
  }

#' Calculate the difference in deuteration between two states for each peptide at each labelling time point
#'
#' @param dataset the HDXMS dataset in a long format with following columns: peptide, label_time, state and deuterium_level
#' @param groups indicator variable that defines the two states to be compared
#'
#'
#' @return A data frame with the calculated differences for each peptide at each time point.
#' @examples
#'
#' #Please check the vignette of the HDXstats package.
#'
#' @export
calcDiff <- function(dataset, groups = NULL) {

  label_time <- peptide <- NULL
  if (is.null(groups))
    stop("No groups to compare provided")
  if (length(groups) != 2)
    stop("Incorrect number of groups to compare provided")

  timepoints <- unique(dataset$label_time[dataset$label_time!=0.0])
  peptides <- unique(dataset$peptide)

  nbTime <- length(timepoints)
  nbGroup <- length(groups)
  nbPeptides <- length(peptides)
  N <- nbPeptides * nbTime

  output <- data.frame(peptide = character(N), label_time = integer(N), comparison = character(N), differences = numeric(N)
    )
  id <- 1

  for (p in peptides) {
    for (tp in timepoints) {
      tmp <- subset(dataset, peptide == p & label_time == tp)
      diffD <- mean(tmp$deuterium_level[tmp$state == groups[1]] - tmp$deuterium_level[tmp$state == groups[2]])
      output$peptide[id] <- p
      output$label_time[id] <- tp
      output$comparison[id] <- paste0(groups[1], "vs", groups[2])
      output$differences[id] <- diffD
      id <- id + 1
    }
  }

  output
}

#' Calculate the global threshold test statistic
#'
#' @param dataset tTe HDXMS dataset in a long format with following columns: peptide, label_time, state and deuterium_level.
#' @param nbRep The number of replicates. Make sure that each peptide has exactly the same number of replicates as specified here.
#' @param group Which protein states need to be compared with each other.
#' @param alpha Significance level. Default is 0.05
#'
#'
#' @return A list with the results of the global threshold test statistic. Each element of the list corresponds to a specific comparison
#' @examples
#'
#' #Please check the vignette of the HDXstats package.
#'
#' @importFrom stats qt var
#' @export
globalThreshold <- function(dataset, nbRep, group, alpha = 0.05) {

  peptide <- state <- NULL
  timepoints <- unique(dataset$label_time[dataset$label_time!=0])
  peptides <- unique(dataset$peptide)

  nbTime <- length(timepoints)
  nbPeptides <- length(peptides)
  nbGroup <- length(group)
  critValue <- abs(qt(alpha / 2, 2 * (nbRep - 1)))

  varTMP <- matrix(NA, nrow = nbGroup, ncol = nbTime)

  for (z in 1:nbGroup) {
    intm0 <- subset(dataset, state == group[z])

    for (k in 1:nbTime) {
      varTMP[z, k] <- var(intm0$deuterium_level[intm0$label_time == timepoints[k]])

    }
  }

  pooledVar <- sum(varTMP) / (nbPeptides * 2 * nbTime)
  SEM <- sqrt(pooledVar / nbRep + pooledVar / nbRep)
  threshold <- critValue * SEM
  intm <- calcDiff(dataset, groups = group)
  intm$threshold <- threshold
  intm$selected <- (abs(intm$differences) >= threshold)
  output <- intm[order(abs(intm$differences), decreasing = TRUE), ]
  output

}

#' Calculate the hybrid significance test statistic
#'
#' @param alpha significance level. Default is 0.05
#' @param dataset the HDXMS dataset in a long format with following columns: peptide, label_time, state and deuterium_level
#' @param nbRep The number of replicates. Make sure that each peptide has exactly the same number of replicates as specified here.
#' @param group Which protein states need to be compared with each other.
#'
#'
#' @return A list with the results of the hybrid significance test statistic. Each element of the list corresponds to a specific comparison
#' @examples
#'
#' #Please check the vignette of the HDXstats package.
#'
#' @importFrom stats t.test
#' @export
hybridSignificanceTest <- function(dataset, nbRep, group, alpha = 0.05) {

  label_time <- peptide <- state <- NULL
  timepoints <-  unique(dataset$label_time[dataset$label_time!=0])
  peptides <- unique(dataset$peptide)

  nbTime <- length(timepoints)
  nbGroup <- length(group)
  nbPeptides <- length(peptides)

  step1 <- globalThreshold(dataset, nbRep, group, alpha)


  step1$significant <- FALSE
  id <- which(step1$selected == TRUE)
  step1$pval <- rep(1, nbPeptides * nbTime)

  for (k in id) {
    selPeptide <- step1$peptide[k]
    selTime <- step1$label_time[k]

    tmp <- subset(dataset, peptide == selPeptide & label_time == selTime & (state == group[1] | state == group[2]))
     step1$pval[k] <- t.test(deuterium_level ~ state, alternative = "two.sided", conf.level = 1 - alpha, tmp)$p.value
    if (step1$pval[k] <= alpha)
      step1$significant[k] <- TRUE
  }

  output <- step1[order(step1$pval), ]

  output
}

#' Compares the difference in deuterium-levels between two states against a fixed user-defined threshold
#'
#' @param dataset The HDXMS dataset in a long format with following columns: peptide, label_time, state and deuterium_level.
#' @param nbRep The number of replicates. Make sure that each peptide has exactly the same number of replicates as specified here.
#' @param group Which protein states need to be compared with each other.
#' @param threshold The minimum difference in deuterium-levels that is considered to be significant.
#'
#'
#' @return A list with the results of the fixed threshold approach. Each element of the list corresponds to a specific comparison
#' @examples
#'
#' #Please check the vignette of the HDXstats package.
#'
#' @importFrom utils combn
#' @export
fixedThreshold <- function(dataset, nbRep, group, threshold) {

  peptide <- state <- NULL
  timepoints <- unique(dataset$label_time)
  peptides <- unique(dataset$peptide)

  nbTime <- length(timepoints)
  nbGroup <- length(group)
  nbPeptides <- length(peptides)

  intm <- calcDiff(dataset, groups = group)
  intm$threshold <- threshold
  intm$selected <- (abs(intm$differences) >= threshold)
  output <- intm[order(abs(intm$differences), decreasing = TRUE), ]

  output
}

#' Calculates the two sample t-test for differential HDX-MS data
#'
#' @param dataset The HDXMS dataset in a long format with following columns: peptide, label_time, state and deuterium_level
#' @param nbRep The number of replicates. Make sure that each peptide has exactly the same number of replicates as specified here.
#' @param group Which protein states need to be compared with each other.
#' @param alpha Significance level. Default is 0.05
#'
#'
#'
#' @return A list with the results of the two sample t-test. Each element of the list corresponds to a specific comparison
#' @examples
#'
#' #Please check the vignette of the HDXstats package.
#'
#' @importFrom stats t.test
#' @export
ttestHDX <- function(dataset, nbRep, group, alpha = 0.05) {

  label_time <- peptide <- state <- NULL
  timepoints <- unique(dataset$label_time[dataset$label_time!=0])
  peptides <- unique(dataset$peptide)

  nbTime <- length(timepoints)
  nbPeptides <- length(peptides)

  N <- nbPeptides * nbTime

  intm <- data.frame(
    peptide = character(N),
    label_time = integer(N),
    comparison = character(N),
    tstat = numeric(N),
    pval = numeric(N),
    significant = logical(N)
    )
  id <- 1

  for (p in peptides) {
    for (tp in timepoints) {
      tmp <- subset(dataset, peptide == p & label_time == tp & (state == group[1] | state == group[2]))
      res <- t.test(deuterium_level ~ state, alternative = "two.sided", conf.level = 1 - alpha, tmp)
      intm$peptide[id] <- p
      intm$label_time[id] <- tp
      intm$comparison[id] <- paste0(group[1], "vs", group[2])
      intm$tstat[id] <- res$statistic
      intm$pval[id] <- res$p.value

      if (res$p.value <= alpha)
        intm$significant[id] <- TRUE
      id <- id + 1
    }
  }
  output <- intm[order(intm$pval), ]

  return(output)
}
