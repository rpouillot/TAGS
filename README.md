# TAGS: evaluation of Tests in the Absence of a Gold Standard 

This repository include the R code for the <a href="http://www.thesearemyapps.com/Tags/">TAGS</a> R shiny application.

The TAGS software can be used to estimate the sensitivity of two or more diagnostic tests in the absence of a gold standard, provided two or more populations with differing prevalences can be cross-classified based on diagnostic test results. 

The algorithm in the TAGS software follows the frequentist paradigm and utilizes Newton-Raphson and EM algorithms to generate maximum likelihood estimates. 

The TAGS software can accommodate not only data of the "2 independent tests, 2 populations"-type, but also higher order combinations of numbers of tests and numbers of populations. Recognizing that in some instances, true prevalence may be known for some populations, TAGS is capable of utilizing "reference population data" (where one or more populations is of known disease status). Parameter estimation using TAGS becomes possible once the number of degrees of freedom given by the data is greater than the number of parameters to be estimated. A goodness-of-fit test and residual correlations, both of which are provided by TAGS output, provide a means of evaluating model adequacy. 

The algorithm includes 2 strong assumptions: (i) diagnostic tests are assumed to be conditionally independent, and (ii) test diagnostic values are considered constant when applied to different populations.

Reference: <a href="https://doi.org/10.1016/S0167-5877(01)00272-0">
R. Pouillot, G. Gerbier, I.A. Gardner. ''TAGS'', a program for the evaluation of test accuracy in the absence of a gold standard, Preventive Veterinary Medicine, 53 (1-2), pp 67-81.</a>
