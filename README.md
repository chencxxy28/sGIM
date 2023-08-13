# sGIM
Integrating External Summary Information in the Presence of Prior Probability Shift


**Background**: Recent years have witnessed a rise in the popularity of information integration without sharing of raw data. By leveraging and incorporating summary information from external sources, internal studies can realize enhanced estimation efficiency and prediction accuracy. However, a noteworthy challenge in utilizing summary-level information is accommodating the inherent heterogeneity across diverse data sources. In this study, we delve into the issue of selection bias present in two cohorts, wherein the bias function depends on the outcome. We introduce a novel semi-parametric constrained optimization-based approach to integrate information within this framework, which has not been extensively explored in existing literature. Our proposed method tackles selection bias by considering the outcome-dependent bias function and effectively addresses the estimation uncertainty associated with summary information from the external source. Our approach facilitates valid inference even in the absence of a known variance-covariance estimate from the external source.


**sGIM** is a semi-parametric constrained optimization-based approach to integrate summary information from the external source while accounting for outcome-dependent selection bias.

# Installation

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("chencxxy28/sGIM")
```

# Vignettes

Please visit [Tutorial](https://chencxxy28.github.io/sGIM/articles/NAME-OF-VIGNETTE.html)



