\name{simulatedData}
\alias{simulatedData}
\docType{data}
\title{
Simulated data from robust linear regression models
}
\description{
Simulated data from the robust linear regression, including \code{simulatedData_gaussian} and \code{simulatedData_quantile}.
}
\usage{
data(simulatedData_gaussian)
}
\details{
We simulated data robust linear regression models
  \deqn{
  y_i = \boldsymbol{X}_i^T\boldsymbol{\alpha} + \boldsymbol{Z}_i^T\boldsymbol{\beta}{\boldsymbol 1}(\boldsymbol{U}_i^T\boldsymbol{\gamma}\geq 0) + \epsilon_i.
  }

\itemize{
\item Y: the response, an \eqn{n}-vector

\item X: the baseline variable with dimension \eqn{n\times p}

\item Z: the grouping difference variable with dimension \eqn{n\times q}

\item U: the grouping variable with dimension \eqn{n\times r}
}

}

\references{
Liu, X. (2022). Change-plane testing in the generalized estimating equations. Manuscript.
}
\examples{
data(simulatedData_gaussian)

y <- data_gaussian$Y[1:5]
x <- dim(data_gaussian$X)
z <- dim(data_gaussian$Z)
u <- dim(data_gaussian$U)

data(simulatedData_quantile)
y <- data_quantile$Y[1:5]
x <- dim(data_quantile$X)
z <- dim(data_quantile$Z)
u <- dim(data_quantile$U)

}
\keyword{datasets}