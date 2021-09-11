[![Build Status](https://app.travis-ci.com/Ali-Mahzarnia/MFSGrp.svg?branch=main)](https://app.travis-ci.com/Ali-Mahzarnia/MFSGrp)
# MFSGrp  
This R package runs the Group Elastic Net (including lasso, ridge, elastic net, and ordinary least square) regression with scalar response values and observed functional covariates. In addition, it penalizes the curvature of the output by implementing a penalty on the second derivative of the estimated coefficient curves. One of the two algorithms of this package is ADMM (mostly developed in C++). ADMM is designed for parallel computations and is only recommended on systems equipped with many strong cores. This algorithm runs parallel on Linux, but it runs serial on Windows. The second algorithm uses the [fGMD](https://github.com/Ali-Mahzarnia/fGMD) package that is built exclusively for this package. The [fGMD](https://github.com/Ali-Mahzarnia/fGMD) package is a heavily modified version of the [gglasso](https://github.com/cran/gglasso) package. The features added to the original gglasso package are: the mixing parameter (alpha) and its net search cross-validation, the curvature penalization for functional regression and its net search cross-validation, the optimized Fortran core function to accelerate the curvature penalization updates, and the progress reports with time estimations. For this package to work, first install [fGMD](https://github.com/Ali-Mahzarnia/fGMD) as instructed below. The [fGMD](https://github.com/Ali-Mahzarnia/fGMD) package does not work independently from this package, and it does not interfere with the functions of the [gglasso](https://github.com/cran/gglasso) package due to slight name differences.
 
 
# Installation
## 1-Dependencies:
In order to have a successful installation, make sure you have all of the required dependencies installed on R. You can install these dependencies with the R commands:  
[Rcpp (>= 1.0.6)](https://cran.r-project.org/web/packages/Rcpp/index.html) ```install.packages("Rcpp")```  
[RcppArmadillo (>= 0.10.2.2.0)](https://cran.r-project.org/web/packages/RcppArmadillo/index.html) ```install.packages("RcppArmadillo")```     
[fda (>= 5.1.9)](https://cran.r-project.org/web/packages/fda/index.html) ```install.packages("fda")```    
[Matrix (>=1.3-2)](https://cran.r-project.org/web/packages/Matrix/index.html) ```install.packages("Matrix")```  
[pbmcapply (>= 1.5.0)](https://cran.r-project.org/web/packages/pbmcapply/index.html) ```install.packages("pbmcapply")```  


## 2-Install [fGMD](https://github.com/Ali-Mahzarnia/fGMD):  
You can install `fGMD` from [GitHub](https://github.com/Ali-Mahzarnia/fGMD) with the R command:
```R
install.packages("https://github.com/Ali-Mahzarnia/fGMD/archive/master.tar.gz", repos = NULL, type="source")
```
## 3-Install [MFSGrp](https://github.com/Ali-Mahzarnia/MFSGrp):
You can install `MFSGrp` from [GitHub](https://github.com/Ali-Mahzarnia/MFSGrp) with the R command:
```R  
install.packages("https://github.com/Ali-Mahzarnia/MFSGrp/archive/refs/heads/main.tar.gz",  repos = NULL, type="source")
```


# Alternative Instalation Methods
## 1- Installing the Development Version: 
This method most likely installs the required dependencies automatically. You can install the development version of `MFSGrp` and `fGMD` via [pacman](https://cran.r-project.org/web/packages/pacman/index.html) with the R commands:
``` R  
install.packages("pacman")
pacman::p_install_gh("Ali-Mahzarnia/fGMD")
pacman::p_install_gh("Ali-Mahzarnia/MFSGrp")
```

## 2-Installing from the Source files:
If the installation fails with the other methods, install the packages from the source files directly with the R commands:
``` R  
# fGMD  
install.packages("https://github.com/Ali-Mahzarnia/fGMD/raw/master/fGMD_1.0.tar.gz",  repos = NULL, type="source")
# MFSGrp
install.packages("https://github.com/Ali-Mahzarnia/MFSGrp/raw/main/MFSGrp_1.0.tar.gz",  repos = NULL, type="source")
```


# Manual and Example:
After installations, you can pull up the manual that includes a simulation example by the R command ```??MFSGrp```. Click on `MFSGrp::MFSGrp` under the help pages for the manual. If the manual cannot be pulled up, first try ```.rs.restartR()```, then try ```??MFSGrp```.   


# Main reference
Ali Mahzarnia, Jun Song. "Multivariate functional group sparse regression: functional predictor selection," Submitted for publication in 2021. [arXiv link](https://arxiv.org/abs/2107.02146)

