% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateTests.R
\name{moderated.Ttest}
\alias{moderated.Ttest}
\title{Calculate the moderated T-test for differential HDX-MS data}
\usage{
moderated.Ttest(
  design = NULL,
  contrastMatrix = NULL,
  mod.test = NULL,
  threshold = NULL,
  adjust.method = "BH",
  dataset = NULL
)
}
\arguments{
\item{design}{the design matrix of the differential HDX-MS experiment, with rows corresponding to the observations and columns to coefficients to be estimated}

\item{contrastMatrix}{contrast matrix of the experiment contrasts between a set of parameters as a numeric matrix. The matrix specifies which pairwise comparisons between the coefficients are to be extracted from the fit}

\item{mod.test}{which moderated test to use: the moderated t-test ('t.test') or the thresholded moderated t-test ('thresholded')}

\item{threshold}{required minimum difference in deuteration-levels to be tested by the thresholded moderated t-test}

\item{adjust.method}{the multiple testing correction method}

\item{dataset}{the HDXMS dataset in a long format with following columns: peptide, label_time, state and deuterium_level}
}
\value{
A list with the results from the moderated T-tests. Each element in the list corresponds to a specific comparison as defined by the contrast matrix.
}
\description{
Calculate the moderated T-test for differential HDX-MS data
}
\details{
The design matrix can be constructed with the model.matrix function of the stats-package, while the function makeContrasts() of the limma-package can be used to create the contrast matrix.
}
\examples{

#Please check the vignette of the HDXstats package.

}
