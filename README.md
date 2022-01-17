# MFSGrp  
This R package runs the Group Elastic Net (including lasso, ridge, elastic net, and ordinary least square) regression with scalar response values and observed functional covariates. In addition, it penalizes the curvature of the output by implementing a penalty on the second derivative of the estimated coefficient curves. One of the two algorithms of this package is ADMM (mostly developed in C++). ADMM is designed for parallel computations and is only recommended on systems equipped with many strong cores. This algorithm runs parallel on Linux, but it runs serial on Windows. The second algorithm uses the [fGMD](https://github.com/Ali-Mahzarnia/fGMD) package that is built exclusively for this package. The [fGMD](https://github.com/Ali-Mahzarnia/fGMD) package is a heavily modified version of the [gglasso](https://github.com/cran/gglasso) package. The features added to the original gglasso package are: the mixing parameter (alpha) and its net search cross-validation, the curvature penalization for functional regression and its net search cross-validation, the optimized Fortran core function to accelerate the curvature penalization updates, and the progress reports with time estimations. For this package to work, first install [fGMD](https://github.com/Ali-Mahzarnia/fGMD) as instructed below. The [fGMD](https://github.com/Ali-Mahzarnia/fGMD) package does not work independently from this package, and it does not interfere with the functions of the [gglasso](https://github.com/cran/gglasso) package due to slight name differences.
 
 
# Installation
It is highly recommended that the latest version of [R](https://www.r-project.org/), [Rstudio](https://www.rstudio.com/products/rstudio/download/), and [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html) are installed before installing the following dependencies on R, and the instructions on their pages are followed so they are activated. Sometimes, the best way to handle errors, such as that of gfortran, due to outdated versions of these three programs, is to uninstall all of these three, then install their latest versions in the above order. 
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

# Examples
The package manual that can be pulled up by  ```??MFSGrp``` on the RStudio console has a thorough explanation and set of examples.


```p=35 # number of functional predictors for each observation
n=200 # sample size
nt=500 # number of recorded time points for each functional covariate
# nt will be reduced to 100 after the inner products are computed below

X= array(NaN, c(p, n, nt)); # Brownian motion
for(j in 1:p){
  for (i in 1:n){
    X[j,i,]=cumsum(rnorm(nt,0,1)) }
}

# true nonzero coefs: beta_5, beta_8, beta_11, and the rest are zeros
# beta_5(t)=sin(3*pi*t), beta_8(t)=sin(5*pi*t/2) and beta_11(t)=t^2
beta5 = function(t){return(sin(3*pi*t/2))}
beta8 = function(t){return(sin(5*pi*t/2))}
beta11=function(t){return(t^2)}
b5=matrix(0, ncol=1, nrow=nt)
b8=matrix(0, ncol=1, nrow=nt)
b11=matrix(0, ncol=1, nrow=nt)

# evaluate population betas on (0,1) at five hundred time points
for(i in 0:nt){
  j=i
  b5[i]=beta5(j/nt)
  b8[i]=beta8(j/nt)
  b11[i]=beta11(j/nt)
}

# evaluate the inner products of Xs and beta 5 and 8 and 11 via Reiman sum
Xb5=matrix(0, ncol=n, nrow=1)
Xb8=matrix(0, ncol=n, nrow=1)
Xb11=matrix(0, ncol=n, nrow=1)

for(j in 1:n){
  Xb5[j]=(X[5,j,] %*%b5)/nt
  Xb8[j]=(X[8,j,]%*%b8)/nt
  Xb11[j]=(X[11,j,]%*%b11)/nt
}
# construct Y
Y=matrix(0, ncol=n, nrow=1)
# standard deviation of the noise term
sd=0.05
# noise term 
eps=matrix(0, ncol=n, nrow=1)
for(n in 1:n){
  eps[, n]=rnorm(1,0,sd)
}
Y=Xb5+Xb8+Xb11+eps
# the algorithm takes care of the intercept in the prediction
Y=Y+3; #intercept


# make the design matrix (pick every 5 elements), here nt becomes 100
X.obs = X[,,(1:100)*nt/100, drop=FALSE]

# observed times scaled to (0,1)
tt=(1:100)/100

# test and train sets (half, half)
trainIndex=sample(1:n, size = round(0.5*n), replace=FALSE)
Ytrain=Y[, trainIndex, drop = FALSE ]
Ytest=Y[, -trainIndex, drop = FALSE ]
Xtrain=X.obs[,trainIndex,, drop=FALSE]
Xtest=X.obs[,-trainIndex,, drop=FALSE]

# The model:
# total 35 functional predictors
# beta_5(t), beta_8(t), beta_11(t) are nonzero, and the others are zero

# plot X^1, ..., X^9   (out of total p=35)
par(mfrow=c(3,3))
for(j in 1:9){
  plot(tt,Xtrain[j,1,],type='l', ylim=c(-30,30), main=paste0("j=",j))
  for(k in 2:10)lines(tt, Xtrain[j,k,])
}
par(mfrow=c(1,1))```

