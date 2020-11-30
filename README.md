# powerdist
Power distribution R Package (Version 0.4.0)
The previous version are named during the development process and will not be published


This is a package for calculate attained power and construted power distributions for unequal cluster-size, cross-sectional stepped-wedge and parallel cluster randomized trials, with or without stratification. Allowed outcome types are: continuous (Gaussian), binary and count.

<p>&nbsp;</p>

# Install the package
### Step one: install _"**devtools**"_ package


```
install.packages("devtools")
```
<p>&nbsp;</p>

### Step two: install _"**powerdist**"_ from Github

```
library(devtools)
install_github("https://github.com/douyangyd/powerdist")
```
<p>&nbsp;</p>

### Step three: Load the package
Remember to load the package before using the package

```
library(powerdist)
```

<p>&nbsp;</p>

# Shinyapp

A shinyapp has been created for non-R user (https://douyang.shinyapps.io/powerdist/)


<p>&nbsp;</p>

# Author information
This package is developed and maintained by Yongdong Ouyang, Liang Xu and Hubert Wong.
The Shinyapp is created by Yongdong Ouyang


<p>&nbsp;</p>

# Contact
We are welcome any suggestions and comments. If you noticed any bugs or errors, please let us know by creating a pull request or contacting Yongdong Ouyang (douyang@cheos.ubc.ca)



<p>&nbsp;</p>

# Reference
Ouyang Y, Karim ME, Gustafson P, et al. Explaining the variation in the attained power of a stepped-wedge trial with unequal cluster sizes. BMC Medical Research Methodology 2020; 20: 166.


Longford NT. Logistic regression with random coefficients. Computational Statistics & Data Analysis 1994; 17: 1–15.

Hussey MA, Hughes JP. Design and analysis of stepped wedge cluster randomized trials. Contemp Clin Trials 2007; 28: 182–191.
