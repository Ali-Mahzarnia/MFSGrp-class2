\name{MFSGrp-package}
\alias{MFSGrp}
\docType{package}
\title{
\packageTitle{MFSGrp}
}
\description{
\packageDescription{MFSGrp}
}
\details{
}
\usage{
MFSGrp(
  Ytrain,
  Xtrain, 
  basisno=5,
  tt, 
  lambda=NULL, 
  alpha=NULL,
  part, 
  Xpred=NULL,
  Ypred=NULL, 
  Silence=FALSE, 
  bspline=FALSE, 
  Penalty=NULL, 
  lambdafactor=0.005,
  nfolds=5, 
  predloss="L2",
  eps = 1e-08,
  maxit = 3e+08,
  nlambda=100,
  forcezero=FALSE, 
  forcezeropar=0.001,
  sixplotnum=1, 
  lambdaderivative=NULL,
  nfolder=5,
  nalpha=9,
  nlamder=10,
  lamdermin=1e-9,
  lamdermax=1e-3,
  alphamin=0 ,
  alphamax=1,
  a=3.7, 
  ADMM=FALSE,
  numcores=NULL,
  rho=1,
  unbalanced=FALSE
  )
}
\arguments{
\item{Ytrain}{a vector of scalar response value of size \eqn{n} in the train set.}

\item{Xtrain}{a three-dimensional matrix (array) of \eqn{p} functional covariates observed at \eqn{nt} time points and sampled \eqn{n} times in the train set. The size of the array is \eqn{p \times n \times  nt}{p*n*nt}. It is recommended that the time series signals are of the same type. For instance, an fMRI data set should be all BOLD contrast activities of the human brain that are recorded \eqn{nt} times during a time period \eqn{T} at \eqn{p} Regions of Interests of the Brain Atlas for \eqn{n} human subjects. In an economic data set example, the covariates should be time series signals that are all measured in percentages or all measured in per capita unit. }

\item{basisno}{the number of bases (of Fourier or B-spline) used in the process of converting the (time series) signals, i.e: observed functional covariates to their functional objects. If \code{bspline=FALSE}, this value must be an odd integer.}

\item{tt}{a vector of  the size \eqn{nt} that includes grid points within 0 to 1. This vector represents the time of observations during the whole period of time \eqn{T}, scaled to the interval \eqn{(0,1)}.}

\item{lambda}{regularization parameter of the (group) lasso penalty term (given alpha=0). By leaving this argument unspecified, the program computes the optimized value in a \code{nfolds}-fold cross-validation on the train set. The cross-validation takes place for values within \code{lambdamax} (the computed maximum lambda by the KKT condition) and \code{lambdafactor}*\code{lambdamax}, for \code{nlambda} possible values, which are equally spaced in the interval based on the \eqn{log()} of the interval boundaries. If the curve of RMSE vs. log(\code{lambda}) is not U-shaped consider a smaller value for \code{lambdafactor} than the default. If \code{ADMM=FALSE}, the value of \code{nfolds} must be greater than \eqn{2}.}

\item{alpha}{the Elastic Net mixing parameter. For the (group) lasso, \code{alpha} becomes \eqn{0} automatically; for the ridge penalty it  becomes 1. If it is 1 the penalty becomes ridge automatically. The only way to incorporate a user-defined \code{alpha} value between 0 and \eqn{1} is to set \code{penalty="gelast"}. When \code{penalty="gelast"}, and leaving this argument unspecified, the program computes the optimized value in a \code{nfolder}-fold net search for a grid of ten \code{lambda}s between the maximum and the minimum \code{lambda} on the train set. The net search cross-validation of \code{alpha} takes place for values within \eqn{0} \eqn{<} \code{alphamin} \eqn{<} \code{alphamax} \eqn{<} \eqn{1} and for \code{nalpha} possible values, which are equally spaced in the interval. If \code{ADMM=FALSE}, the value \code{nfolder} must be greater than \eqn{2}.}


\item{part}{a vector of size \eqn{p} with entries of \code{basisno}. This vector determines the partitions associated to each functional covariate.}

\item{Xpred, Ypred}{the test set versions of \code{Xtrain} and \code{Ytrain}.}

\item{Silence}{by default is \code{FALSE}. Given that it is \code{FALSE}, if in the simulation the population coefficient curves (objects b1,b2,b3,b4,b5 etc.) exist, they will be displayed by green curves along with their estimations. Otherwise, the estimated curves and the number of curves from \eqn{1} to \eqn{p} will be displayed in groups of six figures. There will be \code{sixplotnum} sets of six figures. If \code{sixplotnum = "max"}, there will be as many figures as needed to display all of the estimated nonzero curves. At least six curves will be displayed if \code{Silence=FALSE}. If in the simulation b1,b2,b3,b4,b5, etc. exist, remove them before running the regression on a real data set. See the example below.}

\item{bspline}{by default is \code{FALSE}. If \code{FALSE} the functional conversions and estimations will be processed with the Fourier bases. Otherwise, they will be processed with the B-spline bases.}

 \item{Penalty}{takes five values:
\itemize{ \item{  \code{"glasso"} for the (group) lasso.} 
 \item{\code{"ridge"} for the Tikhonov (aka ridge regression).}
 \item{\code{"OLS"} for the ordinary least square. With this penalty, when also penalizing the curvature \code{ADMM=T} might be faster on both Windows and Linux due to the nature of the \code{C++}. On the other hand, when the second derivative is not penalized, \code{ADMM=F} is accurate and fast.}
 \item{\code{"gelast"} for the (group) elastic net. See \code{alpha} for more information. }
 \item{\code{"gscad"} for the (group) smoothly clipped absolute deviation. Only works when \code{ADMM=TRUE}. It seems to be working and has slightly different results than that of the (group) lasso penalty with different values of \code{a}; however, there is no guarantee for the convergence in theory. This is  because this penalty is non-convex, and the ADMM algorithm is proven to converge for convex penalties. The parameter of the (group) scad is \code{a} which is by default \eqn{3.7}. See the reference below for more information.}}}
 
\item{predloss}{the prediction loss used only in the process of net search and cross-validations. Inherited from the \code{gglasso} package. By default, it is set as the mean square error \code{"L2"}. It can only be changed to other values, such as the mean absolute error \code{"L1"} if \code{ADMM=FALSE}.}

\item{eps}{inherited from the \code{gglasso} package. The threshold of the inner loops. By default, it is \code{1e-08}. If \code{ADMM=TRUE}, the absolute tolerance \code{ABSTOL=eps*1e+4}, and the relative tolerance \code{RELTOL=eps*1e+5}. See the reference below for more information about the stopping criteria of ADMM. }

\item{maxit}{inherited from the \code{gglasso} package. The maximum iterations of the inner loops. By default, it is \code{3e+08}. If \code{ADMM=TRUE}, this value (and its default) is divided by \code{3e+6}.}

\item{forcezero}{by default is \code{FALSE}. If \code{TRUE}, the estimated coefficient curves, the norm of which are less than  \code{forcezeropar}, will become zero curves.}

\item{lambdaderivative}{the regularization parameter of the second derivative penalty of the estimated curves. The value \eqn{0} does not penalize the second derivative. By leaving this argument unspecified, the program computes the optimized value in a \code{nfolder}-fold net search for a grid of ten \code{lambda}s between the maximum and the minimum \code{lambda} on the train set. The net search cross-validation of \code{lambdaderivative} takes place for values within \code{lamdermin} \eqn{<} \code{lamdermax} and for \code{nlamder} possible values, which are equally spaced in the interval based on the \eqn{log()} of the interval boundaries.
If both \code{alpha} and \code{lambdaderivative} need cross-validations, first the \code{lambdaderivative} is tuned with \code{alpha=0}, then the \code{alpha} is tuned with the regularized \code{lambdaderivative}. Both of these processes take place before the regularizations of \code{lambda}. Note that the value \code{nfolder} must be greater than \eqn{2} if \code{ADMM=FALSE}.
}

\item{ADMM}{by default is \code{FALSE}. If \code{TRUE} the algorithm is the Alternating Direction Method of Multiplier. This algorithm runs parallel on Linux, but it runs serial on Windows. It is recommended only on a system equipped with many strong cores. The number of cores to be used is \code{numcores}.  If \code{TRUE}, by leaving \code{numcores} unspecified, the program assigns and uses the maximum number of possible cores on the system. If the input \code{numcores} is greater than the maximum number of available cores, the maximum number of possible cores on the system will be assigned as \code{numcores}. The augmented parameter (the step size) is \code{rho}, which is by default \eqn{1}. See the reference below for more information. }

\item{unbalanced}{by default is \code{FALSE}. If \code{TRUE}, the number of observed time points are different from one functional covariate to another. In this case, \code{Xtrain}, \code{Xtest} and \code{tt} must be lists of \eqn{p}. Each element of \code{Xtrain} and  \code{Xtest} must be a matrix of  \eqn{n \times  nt[j]}{n*nt[j]} for \eqn{j=1,...,p}, and each element of \code{tt} must be a vector of  the size \eqn{nt[j]} that includes grid points within 0 to 1. The \eqn{nt} is a vector of sizes of observed time points for the \eqn{p} functional covariates. See the example below.}
}
\value{
\item{coef}{vector of size \code{basisno}\eqn{*p}. The coefficients of the \eqn{p} estimated functional curves based on \code{basisno} bases.}

\item{predict}{predicted values associated with \code{Ytest}.}

\item{MSEpredict}{mean squared error of the predictions.}

\item{lambda}{the regularized \code{lambda}.}
}
\examples{

p=35 # number of functional predictors for each observation
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
  Xb5[j]=(X[5,j,] \%*\%b5)/nt
  Xb8[j]=(X[8,j,]\%*\%b8)/nt
  Xb11[j]=(X[11,j,]\%*\%b11)/nt
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
par(mfrow=c(1,1))

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
           lamdermax=1e-3, lamdermin = 1e-5)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients
results$lambda # the regularized lambda


# lasso with forcezero 
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
          Ypred=Ytest, Penalty = "glasso" , bspline=TRUE, sixplotnum="max" , 
           lamdermax=1e-3, lamdermin = 1e-5, forcezero=TRUE)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients


# lasso with forcezero and without penalizing the second derivative
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
          Ypred=Ytest, Penalty = "glasso" , bspline=TRUE, sixplotnum="max" , 
           lambdaderivative=0, forcezero=TRUE)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients


# elastic net with different settings
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,Ypred=Ytest, 
          Penalty = "gelast" , bspline=FALSE, nlamder=15, nalpha=11, nfolder=4 )
sum(results$coef==0)/m # number of zero functional coefficients
sqrt(results$MSEpredict) # test Root MSE


# oracle
# as explained in the manual, because of C++, this regression is faster with ADMM, even on Windows 

active=c(5,8,11) # the active set: set of indices of the three nonzero coefficients
results=MFSGrp(Y=Ytrain,X=Xtrain[active,,],basisno=m,tt, part=part,Xpred
              =Xtest[active,,],Ypred=Ytest, Penalty ="OLS", bspline = TRUE, 
              sixplotnum = "max" ,lambdaderivative = NULL, lamdermax = 1e-3,
              lamdermin=1e-5, nfolder=10, ADMM=TRUE)
sqrt(results$MSEpredict)  # test Root MSE


# regression with no shrinkage or sparse penalty, only with second derivative penalty
results=MFSGrp(Ytrain,X=Xtrain,basisno=m,tt, part=part,Xpred=Xtest,Ypred=Ytest, 
          Penalty ="OLS", bspline = TRUE, sixplotnum = 3, ADMM=TRUE)
sqrt(results$MSEpredict) # test Root MSE
sum(results$coef==0)/m  # number of zero functional coefficients (naturally zero)


#  regression with no penalty
results=MFSGrp(Ytrain,X=Xtrain,basisno=m,tt, part=part,Xpred=Xtest,Ypred=Ytest, 
          Penalty ="OLS", bspline = TRUE, sixplotnum = 2, lambdaderivative=0)
sqrt(results$MSEpredict) # test Root MSE
sum(results$coef==0)/m # number of zero functional coefficients (naturally zero)


# a penelized example with ADMM
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
           Ypred=Ytest, Penalty = "glasso" , bspline=TRUE, sixplotnum="max" , 
           lamdermax=1e-3, lamdermin = 1e-5, ADMM=TRUE)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients.
results$lambda # the regularized lambda

# a penelized example with ADMM and a different rho
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
           Ypred=Ytest, Penalty = "glasso" , bspline=TRUE, sixplotnum="max" , 
           lamdermax=1e-3, lamdermin = 1e-5, ADMM=TRUE, rho=0.5)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients
results$lambda # the regularized lambda


# remove the population curves before data application
remove(b5,b8,b11)


# Unbalanced time points:
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



# true nonzero coefs: beta_5, beta_8, beta_11, and the rest are zeros
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
  Xb1[j]=(X[[1]][j,] \%*\%b1)/nt[1]
  Xb2[j]=(X[[2]][j,]\%*\%b2)/nt[2]
  Xb3[j]=(X[[3]][j,]\%*\%b3)/nt[3]
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

# load the library 
library(MFSGrp) 

# run
m=17 #basisino
part=rep(m,p) # partition



# green: true beta  (only beta1, beta2, beta3 are the nonzero functions)
# black: estimated betas

# lasso 
# in order to see all figures after the run, use the "previous plot" arrow on Rstudio
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
               Ypred=Ytest, Penalty = "glasso" , bspline=TRUE, sixplotnum="max" , 
               lamdermax=1e-3, lamdermin = 1e-5, unbalanced=TRUE)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients
results$lambda # the regularized lambda


# lasso with forcezero 
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
               Ypred=Ytest, Penalty = "glasso" , bspline=TRUE, sixplotnum="max" , 
               lamdermax=1e-3, lamdermin = 1e-5, 
               forcezero=TRUE, unbalanced=TRUE)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients

# elastic net with different settings
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,Ypred=Ytest, 
               Penalty = "gelast" , bspline=FALSE, nlamder=15, nalpha=11, 
               nfolder=4,unbalanced=TRUE )
sum(results$coef==0)/m # number of zero functional coefficients
sqrt(results$MSEpredict) # test Root MSE


# oracle
# as explained in the manual, because of C++, this regression is faster with ADMM, even on Windows 

active=c(1,2,3) # the active set: set of indices of the three nonzero coefficients
results=MFSGrp(Y=Ytrain,X=Xtrain[active],basisno=m,tt, part=part,Xpred
               =Xtest[active],Ypred=Ytest, Penalty ="OLS", bspline = TRUE, 
               sixplotnum = "max" ,lambdaderivative = NULL, lamdermax = 1e-3,
               lamdermin=1e-5, nfolder=10, ADMM=TRUE,unbalanced=TRUE)
sqrt(results$MSEpredict)  # test Root MSE


# regression with no shrinkage or sparse penalty, only with second derivative penalty
results=MFSGrp(Ytrain,X=Xtrain,basisno=m,tt, part=part,Xpred=Xtest,Ypred=Ytest, 
               Penalty ="OLS", bspline = TRUE, 
               sixplotnum = 3, ADMM=TRUE,unbalanced=TRUE)
sqrt(results$MSEpredict) # test Root MSE
sum(results$coef==0)/m  # number of zero functional coefficients (naturally zero)


#  regression with no penalty
results=MFSGrp(Ytrain,X=Xtrain,basisno=m,tt, part=part,Xpred=Xtest,Ypred=Ytest, 
               Penalty ="OLS", bspline = TRUE, 
               sixplotnum = 2, lambdaderivative=0, unbalanced=TRUE)
sqrt(results$MSEpredict) # test Root MSE
sum(results$coef==0)/m # number of zero functional coefficients (naturally zero)


# a penelized example with ADMM
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
               Ypred=Ytest, Penalty = "glasso" , bspline=TRUE, sixplotnum="max" , 
               lamdermax=1e-3, lamdermin = 1e-5, ADMM=TRUE, unbalanced=TRUE)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients.
results$lambda # the regularized lambda

# a penelized example with ADMM and a different rho
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
               Ypred=Ytest, Penalty = "glasso" , bspline=TRUE, sixplotnum="max" , 
               lamdermax=1e-3, lamdermin = 1e-5, ADMM=TRUE, rho=0.5, unbalanced=TRUE)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients
results$lambda # the regularized lambda


# remove the population curves before data application
remove(b1,b2,b3)


}
\references{
Ali Mahzarnia, Jun Song. "Multivariate functional group sparse regression: functional predictor selection," Submitted for publication in  2021.\cr URL: \url{https://arxiv.org/abs/2107.02146}\cr

Yang, Y. and Zou, H. (2015), "A Fast Unified Algorithm for
Computing Group-Lasso Penalized Learning Problems," \emph{Statistics and
Computing}. 25(6), 1129-1141.

Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized linear models via coordinate descent," \emph{Journal of
Statistical Software}, 33, 1.

S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein (2010) "Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers"  \emph{Foundations and Trends R in Machine Learning}, 3,1.

Fan, J and Li, R. (2001) "Variable Selection via Nonconcave Penalized Likelihood and its Oracle Properties" \emph{Journal of the American Statistical Association}, 96, 456. 


}
\seealso{
\code{fGMD} package. Install the \code{fGMD} package in order to have \code{MFSGrp} package work.
\cr URL: \url{https://github.com/Ali-Mahzarnia/fGMD}\cr
}
\author{
Ali Mahzarnia and Jun Song\cr Maintainer: Ali Mahzarnia <amahzarn@uncc.edu> <ali.mahzarnia@duke.edu>
}


