\name{estolsmcp}
\alias{estolsmcp}
\title{
  Estimation in the Linear Models with subgroups
}
\description{
  Provide estimators of coefficients in the linear models with subgroups.
}

\usage{
  estolsmcp(data, ng = 2, smooth = 'sigmoid', isBoot = FALSE, isWB = FALSE,
                            h = NULL, maxIter = 100, tol = 0.00001, B = 1000)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{X} (baseline variable), \eqn{Z} (grouping difference variable), and \eqn{U} (grouping variable).}

  \item{ng}{An integer, which is the number of change-planes. Default is \code{ng = 2}.}

  \item{smooth}{The smooth function. Either "sigmoid" (the default), "pnorm", or "mixnorm", see details below.}

  \item{isBoot}{A bool value.  A bootstrap method is used if \code{isBoot = TRUE}. Default is \code{isBoot = FALSE}.}

  \item{isWB}{A bool value. The wild bootstrap method is used if \code{isWB = TRUE}. Default is \code{isWB = FALSE}.}

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
  \item{ha}{Estimator of the threshholds \eqn{\{a_1, \cdots, a_S\}}, where \eqn{S} equals to \code{ng}.}
  \item{sigma2}{Estimator of the error's variance.}
}

\details{
  The linear models
  \deqn{
  y_i = \boldsymbol{X}_i^T\boldsymbol{\alpha} + \boldsymbol{Z}_i^T\sum_{s=1}^{S}\boldsymbol{\beta}_s{\boldsymbol 1}(U_i+\boldsymbol{U}_{2i}^T\boldsymbol{\theta}_{-1}\geq a_s) + \epsilon_i.
  }
  with the identifiable restraint that \eqn{a_1< a_2<\cdots< a_S}.
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
  fit <- estolsmcp(data = data_gaussian)
  fit$alpha

  data(simulatedData_quantile)
  fit <- estolsmcp(data = data_quantile)
  fit$alpha
}