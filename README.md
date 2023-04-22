# wasthub
R package "wasthub" for learning classifier and calculating p-value of the test statistic for subgroup detecting in the robust linear regression models. In the paper Liu (2023), we propose a novel test statistic by taking the weighted average of the squared score test statistic (WASTHUB) over the nuisance parametric space. The proposed test statistics not only improve power, but also save dramatically computational time. Many common and useful models are considered, including models with change point or change plane. We propose a novel U-like test statistic to detect multiple change planes in the framework of EE.

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("xliusufe/wasthub")

# Usage

   - [x] [wasthub-manual.pdf](https://github.com/xliusufe/wasthub/blob/master/inst/wasthub-manual.pdf) ---------- Details of the usage of the package.
# Example
    library(wasthub)

    ## Testing

    data(simulatedData_gaussian)
    pvals <- pvalhuber(data = data_gaussian, method = "wast")
    pvals

    data(simulatedData_quantile)
    pvals <- pvalhuber(data = data_quantile, method = "wast")
    pvals


    ## Estimation

    data(simulatedData_gaussian)
    fit <- esthubcp(data = data_gaussian, method = "adaptive")
    fit$alpha

    data(simulatedData_quantile)
    fit <- esthubcp(data = data_quantile, method = "adaptive")
    fit$alpha
    fit$beta

# References
Andrews, D. W. K. and Ploberger, W. (1994). Optimal tests when a nuisance parameter is
present only under the alternative. Econometrica, 62(6):1383-1414.

Davies, R. B. (1977). Hypothesis testing when a nuisance parameter is present only under the
alternative. Biometrika, 64(2):247-254.

Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is present only under the
alternative. Biometrika, 74(1):33-43.

Fan, A., Rui, S., and Lu, W. (2017). Change-plane analysis for subgroup detection and sample
size calculation. Journal of the American Statistical Association, 112(518):769-778.

Huang, Y., Cho, J., and Fong, Y. (2021). Threshold-based subgroup testing in logistic regression
models in two phase sampling designs. Journal of the Royal Statistical Society: Series C. 70(2):291-311.

Liu, X. (2023). Change-plane testing in the generalized estimating equations. Manuscript.

Liu. X. (2023). Robust change-plane testing and learning based on Huber loss. Manuscript.

# Development
This R package is developed by Xu Liu (liu.xu@sufe.edu.cn).
