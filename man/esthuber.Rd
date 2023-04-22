\name{esthub}
\alias{esthub}
\title{
  Estimation in Robust Linear Models with subgroups
}
\description{
  Provide estimators of coefficients in robust linear models with subgroups.
}

\usage{
  esthub(x, y, method = "adaptive", maxIter = 100,
                tol = 0.00001, M = 100, tauk = NULL)
}

\arguments{
  \item{x}{A matrix in \eqn{R^{n\times p}}. The design matrix.}

  \item{y}{A  vector in \eqn{R^{n}}. The response.}

  \item{method}{There are there methods, including 'adaptive', 'CV' and preset \code{qq = 1.345*median( abs(y - median(y)) )/qnorm(0.75)}. Default is \code{method = 'adaptive'}.}

  \item{maxIter}{An integer, the maximum number of iterations. Default is \code{maxIter = 100}.}

  \item{tol}{Convergence threshhold. Default is \code{tol = 0.00001}.}

  \item{M}{An integer, the length of \code{tauk}. Default is \code{M = 100}.}

  \item{tauk}{A numeric vector, which is the preset \eqn{\tau} if \code{method = 'CV'}. Default is \code{tauk = NULL}, which is \code{tauk 	= seq(0.01,10,length.out = M)}.}
}

\value{
  \item{alpha}{Estimator of the baseline parameter \eqn{\boldsymbol{\alpha}}.}
  \item{tau}{Optimal \eqn{\tau}.}
  \item{sigma2}{Estimator of the error's variance.}
}

\details{
  Robust linear models
  \deqn{
  y_i = \boldsymbol{X}_i^T\boldsymbol{\alpha} + \epsilon_i.
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
  fit <- esthub(x = data_gaussian$X, y = data_gaussian$Y, method = "adaptive")
  fit$alpha

  data(simulatedData_quantile)
  fit <- esthub(x = data_quantile$X, y = data_quantile$Y, method = "adaptive")
  fit$alpha
}