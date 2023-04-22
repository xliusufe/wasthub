\name{esthubmcpBoot}
\alias{esthubmcpBoot}
\title{
  Estimation in Robust Linear Models with subgroups
}
\description{
  Provide estimators of coefficients in robust linear models with subgroups.
}

\usage{
  esthubmcpBoot(data, ng = 2, method = "adaptive", smooth = 'sigmoid',
        weights = 'exponential', h = NULL, maxIter = 100, tol = 0.00001, B = 1000)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{X} (baseline variable), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}

  \item{ng}{An integer, which is the number of change-planes. Default is \code{ng = 2}.}

  \item{method}{There are there methods, including 'adaptive', 'CV' and preset \code{qq = 1.345*median( abs(y - median(y)) )/qnorm(0.75)}. Default is \code{method = 'adaptive'}.}

  \item{smooth}{The smooth function. Either "sigmoid" (the default), "pnorm", or "mixnorm", see details below.}

  \item{weights}{The weights. Either "exponential" (the default), "norm", or "bernoulli", see details below.}

  \item{h}{A numeric number, which is the bandwidth in the smooth function. Default is \code{h = NULL}, which is \code{h = log(n)/sqrt(n)}.}

  \item{maxIter}{An integer, the maximum number of iterations. Default is \code{maxIter = 100}.}

  \item{tol}{Convergence threshhold. Default is \code{tol = 0.00001}.}

  \item{B}{An integer, the number of bootstrap sample sets. Default is \code{B = 1000}.}
}

\value{
  \item{alpha}{Estimator of the baseline parameter \eqn{\boldsymbol{\alpha}}.}
  \item{beta}{Estimator of the grouping difference parameter \eqn{\boldsymbol{\beta}}.}
  \item{theta}{Estimator of the grouping parameter \eqn{\boldsymbol{\theta}}.}
  \item{delta}{A vector with length \eqn{n}. Estimator of the indicator function \eqn{I(\boldsymbol{U}^T\boldsymbol{\theta}\geq 0)}.}
  \item{tau}{Optimal \eqn{\tau}.}
  \item{ha}{Estimator of the threshholds \eqn{\{a_1, \cdots, a_S\}}, where \eqn{S} equals to \code{ng}.}
  \item{sigma2}{Estimator of the error's variance.}
  \item{std}{A vector with length \eqn{p+q+r-1}. The standard deviation (\code{sd}) of parameter \eqn{(\boldsymbol{\alpha}^T, \boldsymbol{\beta}^T,\boldsymbol{\gamma}_{-1}^T)^T},
  where \eqn{\boldsymbol{\gamma}_{-1}=(\gamma_2,\cdots,\gamma_r)^T}.}
  \item{alphaB}{A matrix in \eqn{R^{p_1\times B}}, each colum of which is the estimator of the baseline parameter \eqn{\boldsymbol{\alpha}}.}
  \item{betaB}{A matrix in \eqn{R^{p_2\times B}}, each colum of which is the estimator of the grouping difference parameter \eqn{\boldsymbol{\beta}}.}
  \item{thetaB}{A matrix in \eqn{R^{p_3\times B}}, each colum of which is the estimator of the grouping parameter \eqn{\boldsymbol{\theta}}.}
  \item{haB}{A matrix in \eqn{R^{S\times B}}, each colum of which is the estimator of the threshholds \eqn{\{a_1, \cdots, a_S\}}, where \eqn{S} equals to \code{ng}.}
}

\details{
  Robust linear models
  \deqn{
  y_i = \boldsymbol{X}_i^T\boldsymbol{\alpha} + \boldsymbol{Z}_i^T\boldsymbol{\beta}{\boldsymbol 1}(\boldsymbol{U}_i^T\boldsymbol{\theta}\geq 0) + \epsilon_i.
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
  fit <- esthubmcpBoot(data = data_gaussian)
  fit$alpha

  data(simulatedData_quantile)
  fit <- esthubmcpBoot(data = data_quantile)
  fit$beta
}