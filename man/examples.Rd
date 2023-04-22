\name{exams}
\alias{exams}
\title{
  Examples for Subgroup Test in Robust Linear Models
}
\description{
  Examples for robust test of the linear regression models.
}

\usage{exams(method = "wast", B = 1000, K = 1000)
}

\arguments{
  \item{method}{There are there methods, including the proposed 'wast', 'sst', and 'slrt'.}

  \item{B}{An integer, the number of bootstrap samples.}

  \item{K}{An integer, the number of threshold values for 'sst' and 'slrt'.}
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}
}


\references{
Andrews, D. W. K. and Ploberger, W. (1994). Optimal tests when a nuisance parameter is
present only under the alternative. Econometrica, 62(6):1383-1414.

Fan, A., Rui, S., and Lu, W. (2017). Change-plane analysis for subgroup detection and sample
size calculation. Journal of the American Statistical Association, 112(518):769-778.

Huang, Y., Cho, J., and Fong, Y. (2021). Threshold-based subgroup testing in logistic regression
models in two phase sampling designs. Journal of the Royal Statistical Society: Series C. 291-311.

Liu, X. (2022). Change-plane testing in the robust estimating equations. Manuscript.
}

\examples{
  pvals <- exams(method = "wast")
  pvals

  pvals <- exams(method = "wastcv")
  pvals
}