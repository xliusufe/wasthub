\name{pvalhuber}
\alias{pvalhuber}
\title{
  P-value for Subgroup Test in Robust Linear Models
}
\description{
  Provide p-value for subgroup test in robust linear models, including four methods 'wast', 'wastcv', 'sst', and 'slrt'.
}

\usage{
  pvalhuber(data, method = "wast", isWB = FALSE, B = 1000, K = 1000, typewgt = 1,
        isBeta = FALSE, shape1 = 1, shape2 = 1, N0 = 5000, MU = NULL, ZK = NULL, isMed = TRUE)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{X} (baseline variable), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}

  \item{method}{There are there methods, including the proposed 'wast', 'wastcv', 'wastapprox', 'sst', and 'slrt'.}

  \item{isWB}{A bool value. The wild bootstrap method is used if \code{isWB = TRUE}. Default is \code{isWB = FALSE}.}

  \item{B}{An integer, the number of bootstrap samples.}

  \item{K}{An integer, the number of threshold values for 'sst' and 'slrt' (Default is \code{K = 1000}), or the length of preset \eqn{\tau} for 'wastcv' (Default is \code{K = 100}).}

  \item{typewgt}{An integer. The standard multivariate normal density is chosen as the grouping parameter \eqn{\theta}'s' weight if \code{typewgt = 1}; and uniform distribution on \eqn{S^{q-1}} is chosen otherwise. Default is \code{typewgt = 1}.}

  \item{isBeta}{A bool value. The weight \eqn{w(\gamma)} is chosen to be Beta distribution if \code{isBeta = TRUE}, which can be used if the grouping difference variable is bounded in \eqn{[0,1]}. Default is \code{FALSE}.}

  \item{shape1}{The first parameter of Best distribution if \code{isBeta = TRUE}.}

  \item{shape2}{The second parameter of Best distribution if \code{isBeta = TRUE}.}

  \item{N0}{An integer, the number of samples to approximate \eqn{\omega_{ij}} for 'wastapprox'. Default is \code{N0 = 5000}. }

  \item{MU}{A vector with same length as \eqn{\boldsymbol{U}}, which is the mean of weight to approximate \eqn{\omega_{ij}} for 'wastapprox'. Default is \code{MU = NULL}, in which \code{MU = runif(p3) - 0.5}. }

  \item{ZK}{A vector with length \code{N0}, which is normal sample to approximate \eqn{\omega_{ij}} for 'wastapprox'. Default is \code{ZK = NULL}, in which \code{ZK	= rnorm(N0)}. }

  \item{isMed}{A bool value. Residual's median \code{tau0*median( abs(residual - median(residual)) )/qnorm(0.75)} is used to get initial \eqn{\tau} if \code{isAda = TRUE}. Response's median \code{tau0*median( abs(y - median(y)) )/qnorm(0.75)} is used to get initial \eqn{\tau} if \code{isAda = FALSE}. Default is \code{isAda = TRUE}.}
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}
}

\details{
  Generalized linear models
  \deqn{
  y_i = \boldsymbol{X}_i^T\boldsymbol{\alpha} + \boldsymbol{Z}_i^T\boldsymbol{\beta}{\boldsymbol 1}(\boldsymbol{U}_i^T\boldsymbol{\gamma}\geq 0) + \epsilon_i.
  }

The hypothesis test problem is
\deqn{
  H_0: \boldsymbol{\beta}=\boldsymbol{0}\quad versus\quad H_1: \boldsymbol{\beta}\neq \boldsymbol{0}.
}
}

\references{
Andrews, D. W. K. and Ploberger, W. (1994). Optimal tests when a nuisance parameter is
present only under the alternative. Econometrica, 62(6):1383-1414.

Fan, A., Rui, S., and Lu, W. (2017). Change-plane analysis for subgroup detection and sample
size calculation. Journal of the American Statistical Association, 112(518):769-778.

Huang, Y., Cho, J., and Fong, Y. (2021). Threshold-based subgroup testing in logistic regression
models in two phase sampling designs. Journal of the Royal Statistical Society: Series C. 291-311.

Liu, X. (2023). Change-plane testing in the generalized estimating equations. Manuscript.

Liu. X. (2023). Robust change-plane testing and learning based on Huber loss. Manuscript.
}

\examples{
  data(simulatedData_gaussian)
  pvals <- pvalhuber(data = data_gaussian, method = "wast")
  pvals

  data(simulatedData_quantile)
  pvals <- pvalhuber(data = data_quantile, method = "wast")
  pvals

}