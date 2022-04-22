# MFSGrp-class2
# This is a different version of [MFSGrp](https://github.com/Ali-Mahzarnia/MFSGrp) with a different unbalanced example:
In this version, the unbalanced case is different. The number of observed time points are different across p functional covariates while they are the same for each functional covariate across observations. In the example of [MFSGrp](https://github.com/Ali-Mahzarnia/MFSGrp), the observed data points are uneven for each observation while they are the same for the p fucntional covariates within that observation.    

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
## 3-Install [MFSGrp-class2](https://github.com/Ali-Mahzarnia/MFSGrp-class2):
You can install `MFSGrp` from [GitHub](https://github.com/Ali-Mahzarnia/MFSGrp-class2) with the R command:
```R  
install.packages("https://github.com/Ali-Mahzarnia/MFSGrp-class2/archive/refs/heads/main.tar.gz",  repos = NULL, type="source")
```


# Alternative Instalation Methods
## 1- Installing the Development Version: 
This method most likely installs the required dependencies automatically. You can install the development version of `MFSGrp` and `fGMD` via [pacman](https://cran.r-project.org/web/packages/pacman/index.html) with the R commands:
``` R  
install.packages("pacman")
pacman::p_install_gh("Ali-Mahzarnia/fGMD")
pacman::p_install_gh("Ali-Mahzarnia/MFSGrp-class2")
```

## 2-Installing from the Source files:
If the installation fails with the other methods, install the packages from the source files directly with the R commands:
``` R  
# fGMD  
install.packages("https://github.com/Ali-Mahzarnia/fGMD/raw/master/fGMD_1.0.tar.gz",  repos = NULL, type="source")
# MFSGrp
install.packages("https://github.com/Ali-Mahzarnia/MFSGrp-class2/raw/main/MFSGrp_1.0.tar.gz",  repos = NULL, type="source")
```


# Manual and Example:
After installations, you can pull up the manual that includes a simulation example by the R command ```??MFSGrp```. Click on `MFSGrp::MFSGrp` under the help pages for the manual. If the manual cannot be pulled up, first try ```.rs.restartR()```, then try ```??MFSGrp```.   


# Examples
The package manual that can be pulled up by  ```??MFSGrp``` on the RStudio console has a thorough explanation and set of examples.


```R 
p=35 # number of functional predictors for each observation
n=200 # sample size
nt=500 # number of recorded time points for each functional covariate
# nt will be reduced to 100 after the inner products are computed below

X= array(NaN, c(p, n, nt)); # Brownian motion
for(j in 1:p){
  for (i in 1:n){
    X[j,i,]=cumsum(rnorm(nt,0,1)) }
}
```
<img src="https://render.githubusercontent.com/render/math?math=X^i_j(t)"> for <img src="https://render.githubusercontent.com/render/math?math=i=1,\cdots,n">, <img src="https://render.githubusercontent.com/render/math?math=j=1,\cdots,p">, and <img src="https://render.githubusercontent.com/render/math?math=t \in \{ t_1, \cdots, t_{nt} \}">

```R 
# true nonzero coefs: beta_5, beta_8, beta_11, and the rest are zeros
# beta_5(t)=sin(3*pi*t), beta_8(t)=sin(5*pi*t/2) and beta_11(t)=t^2
beta5 = function(t){return(sin(3*pi*t/2))}
```
<img src="https://render.githubusercontent.com/render/math?math=\beta^5(t)=\sin(\frac{5 \pi t}{2})"> for <img src="https://render.githubusercontent.com/render/math?math=t \in \[0,1 ]">.

```R 
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
``` 
<img src="https://render.githubusercontent.com/render/math?math=< X^5_i, \beta^5 >= \sum_{k=1}^{nt} \frac{X^5_i(t_k) \beta^5(t_k)}{nt} ">

```R 
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
```
<img src="https://render.githubusercontent.com/render/math?math=Y_i=\sum_{j=1}^p <X^j_i,\beta^j>%2B3%2B\epsilon">


```R 
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
```
![Figures of the first 9 covariates](https://github.com/Ali-Mahzarnia/222/blob/main/readme%20pics/a.png)

```R
par(mfrow=c(1,1))

# load the library 
library(MFSGrp) 

# run
m=17 #basisino
part=rep(m,p) # partition

# green: true beta  (only beta5, beta8, beta11 are the nonzero functions)
# black: estimated betas

# lasso 
```
<img src="https://render.githubusercontent.com/render/math?math=\arg\min_{\beta} \frac{1}{2} || Y-<X,\beta>||^2 %2B \frac{\lambda_{der}}{2} || \beta''||^2 %2B  \lambda \sum_{j} || \beta^j||">

```R 
# in order to see all figures after the run, use the "previous plot" arrow on Rstudio
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
           Ypred=Ytest, Penalty = "glasso" , bspline=TRUE, sixplotnum="max" , 
           lamdermax=1e-3, lamdermin = 1e-5)
```           
Chosen lambdader is 0.0005994843 and Maximum Estimated Time: 12.2067  seconds
```R
sqrt(results$MSEpredict)  # test Root MSE
```
0.5452372
```R
sum(results$coef==0)/m    # number of zero functional coefficients among 35
```
30

```R 
results$lambda # the regularized lambda
```
0.3617951

![Figures of the nonzero funcitonal predictors ](https://github.com/Ali-Mahzarnia/222/blob/main/readme%20pics/a2.png) ![](https://github.com/Ali-Mahzarnia/222/blob/main/readme%20pics/a3.png) ![](https://github.com/Ali-Mahzarnia/222/blob/main/readme%20pics/b.png)

# Unbalanced time points:
The package manual has multiple examples about this topic.



```R
# We'd generate functional data that are observed at a different number of time points.
# For example, the first covariate is observed at  91 time points (or 91*5) and the second at 112 (or 112*5) 
p=19 # number of functional predictors for each observation
n=300 # sample size
nt=5*sample(80:120, p, replace=T)  # vector of number of recorded *5 and later we pick 1 of each 5
#time points for each functional covariate that are different 
#and here are chosen randmoly
X=vector("list",p) # X is a list of p, each item is a mtrix of n*nt(j): X[[1]],...X[[p]]
#(p, n, nt)
for(j in 1:p){# first make matrices of n*nt(j) for each X[j] : j=1,...,p
  X[[j]]=matrix(0, n, nt[j]);
} 

for(j in 1:p){# Brownian motion
  for (i in 1:n){
    X[[j]][i,]=as.numeric(cumsum(rnorm(nt[j],0,1))) }
}
```

<img src="https://render.githubusercontent.com/render/math?math=X^i_j(t)"> for <img src="https://render.githubusercontent.com/render/math?math=i=1,\cdots,n">, <img src="https://render.githubusercontent.com/render/math?math=j=1,\cdots,p">, and <img src="https://render.githubusercontent.com/render/math?math=t \in \{ t_1, \cdots, t_{nt(j)} \}">. Note that <img src="https://render.githubusercontent.com/render/math?math=nt(j)\}"> depends on <img src="https://render.githubusercontent.com/render/math?math=j=1, \cdots, p">, becasue the observed time points are different from one covariate to another.


```R
# true nonzero coefs: beta_1, beta_2, beta_3, and the rest are zeros
# beta_1(t)=sin(3*pi*t), beta_2(t)=sin(5*pi*t/2) and beta_3(t)=t^2
beta1 = function(t){return(sin(3*pi*t/2))}
beta2 = function(t){return(sin(5*pi*t/2))}
beta3=function(t){return(t^2)}
b1=matrix(0, ncol=1, nrow=nt[1])
b2=matrix(0, ncol=1, nrow=nt[2])
b3=matrix(0, ncol=1, nrow=nt[3])

# evaluate population betas on (0,1)  grids
for (d in 1:3) {
  temp=get(paste0("b", d, sep=""))
  for(k in 1:nt[d]){
    temp[k]=get(paste("beta", d, sep=""))(k/nt[d]);}
  assign((paste0("b", d, sep="")),temp )
}
# evaluate the inner products of Xs and beta 5 and 8 and 11 via Reiman sum
Xb1=matrix(0, ncol=n, nrow=1)
Xb2=matrix(0, ncol=n, nrow=1)
Xb3=matrix(0, ncol=n, nrow=1)

for(j in 1:n){
  Xb1[j]=(X[[1]][j,] %*%b1)/nt[1]
  Xb2[j]=(X[[2]][j,]%*%b2)/nt[2]
  Xb3[j]=(X[[3]][j,]%*%b3)/nt[3]
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
Y=Xb1+Xb2+Xb3+eps
# the algorithm takes care of the intercept in the prediction
Y=Y+3; #intercept


# make the design matrix (pick every 5 elements), here nt[j] becomes nt[j]/5
X.obs=vector("list",p);
for (j in 1:p) { 
  X.obs[[j]]=X[[j]][,seq(5, nt[j], 5), drop=FALSE]; 
  nt[j]=dim(X.obs[[j]])[2]
}
#X.obs = X[,,(1:100)*nt/100, drop=FALSE]
#X.obs=X
# observed times scaled to (0,1)
tt=vector("list",p);
for (j in 1:p) {
  tt[[j]]= (1:nt[j])/nt[j]
}


# test and train sets (half, half)
trainIndex=sample(1:n, size = round(0.5*n), replace=FALSE)
Ytrain=Y[, trainIndex, drop = FALSE ]
Ytest=Y[, -trainIndex, drop = FALSE ]
Xtrain=vector("list",p); Xtest=vector("list",p);
for (j in 1:p) {
  Xtrain[[j]]=X.obs[[j]][trainIndex,, drop=FALSE]
  Xtest[[j]]=X.obs[[j]][-trainIndex,, drop=FALSE]
  
}

# The model:
# total 3 functional predictors
# beta_1(t), beta_2(t), beta_3(t) are nonzero, and the others are zero

# plot X^1, ..., X^9   (out of total p=35)
par(mfrow=c(3,3))
for(j in 1:9){
  plot(tt[[j]],Xtrain[[j]][1,],type='l', ylim=c(-30,30), main=paste0("j=",j))
  for(k in 2:10)lines(tt[[j]], Xtrain[[j]][k,])
}
par(mfrow=c(1,1))
```

![Figures of the first 9 covariates](https://github.com/Ali-Mahzarnia/222/blob/main/readme%20pics/c.png)

```R
# load the library 
library(MFSGrp) 

# run
m=17 #basisino
part=rep(m,p) # partition



# green: true beta  (only beta5, beta8, beta11 are the nonzero functions)
# black: estimated betas

# lasso 
# in order to see all figures after the run, use the "previous plot" arrow on Rstudio
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
               Ypred=Ytest, Penalty = "glasso" , bspline=TRUE, sixplotnum="max" , 
               lamdermax=1e-3, lamdermin = 1e-5, unbalanced=TRUE)
```
Chosen lambdader is 0.0005994843 and Maximum Estimated Time: 5.701051  seconds            
```R 
sqrt(results$MSEpredict)  # test Root MSE
```
0.2485125
```R 
sum(results$coef==0)/m    # number of zero functional coefficients out of 19
```
16

```R 
results$lambda # the regularized lambda
```
 0.4207415
 
 ![Figures of the nonzero funcitonal predictors ](https://github.com/Ali-Mahzarnia/222/blob/main/readme%20pics/c2.png) ![](https://github.com/Ali-Mahzarnia/222/blob/main/readme%20pics/c3.png) ![](https://github.com/Ali-Mahzarnia/222/blob/main/readme%20pics/c4.png)


# Main reference
Mahzarnia A, Song J (2022) Multivariate functional group sparse regression: Functional predictor selection. PLoS ONE 17(4): e0265940. [DOI link](https://doi.org/10.1371/journal.pone.0265940)
