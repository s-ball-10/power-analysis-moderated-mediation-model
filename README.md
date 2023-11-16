# power-analysis-multilevelmodel
This repository provides the code for doing power analyis for moderated mediation models using monte carlo simulations.

# Background
In order for an experiment to be meaningful, an appropriate sample size needs to be determined. Too small a sample size can result in unreliable effect sizes, inflated false significant effects, and low reproducibility due to low power (Rousselet, 2022). Although analytical methods to compute the adequate sample size for general linear models exist (e.g. G*Power by Faul et al., 2007), these approaches are not applicable to more complex models such as moderated mediation models. This is why the literature suggests performing Monte Carlo simulations in such situations (Donnelly et al., 2021). 

Using the package lavaan for structural equation modeling (Rosseel, 2012), the code provides a function that estimates the power at different sample sizes with the standard 0.05 alpha error probability based on 1000 simulations per tested sample size. The most critical assumption in the simulation are the expected effect sizes. You can obtain these expected effect sizes from other related studies, for example. The necessary sample size is then determined by the smallest effect size of interest.

## Cited papers

Donnelly, S., Jorgensen, T. D., & Rudolph, C. W. (2021). Power analysis for con- ditional indirect effects: A tutorial for conducting monte carlo simulations with categorical exogenous variables.

Faul, F., Erdfelder, E., Lang, A.-G., & Buchner, A. (2007). G* power 3: A flexi- ble statistical power analysis program for the social, behavioral, and biomedical sciences. Behavior research methods, 39(2), 175–191.

Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. Journal of Statistical Software, 48(2), 1–36.

Rousselet, G. (2022). Problems with small sample sizes. https://garstats.wordp ress.com/2017/02/04/small-sample-sizes/. Accessed March 15, 2022.

