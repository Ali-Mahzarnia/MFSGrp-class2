
#setwd(paste0(getwd(),'/R'))
#f=list.files()
#lapply(f,source);


get.fd = function(xraw,tt,basisname=NULL,nbasis=15,ncv=10,basis=NULL){
  nt=dim(xraw)[1];n=dim(xraw)[2];p=dim(xraw)[3]
  if(is.na(p)) p = 1
  tmin=min(tt);tmax=max(tt)
  if(is.null(basis)){
    if(is.null(basisname)) stop('error: either basis or basisname required.')
    if(basisname=='bspline'){
      basis=fda::create.bspline.basis(c(tmin,tmax),nbasis=nbasis)
      kt =fda::bsplinepen(basis,0)
    }else if(basisname=='fourier'){
      if(nbasis %% 2 ==0) nbasis=nbasis +1   # odd number of basis functions
      basis=fda::create.fourier.basis(c(tmin,tmax), nbasis=nbasis)
      kt=diag(nbasis)
    }
  }
  nbasis=basis$nbasis
  
  if(p==1){
    lambda_all<- exp(seq(log(5*10^(-15)), log(5*10^(0)), len=ncv))
    nlam<-length(lambda_all)
    gcv_all<-matrix(0,n, nlam)
    ms_all_m<-array(0, c(nbasis, n,nlam))
    for(i in 1:nlam){
      lambda=lambda_all[i]
      fd_par<-fda::fdPar(fdobj=basis,Lfdobj=2,lambda=lambda)
      tmp<-fda::smooth.basis(tt,xraw,fdParobj=fd_par)
      gcv_all[,i]<-tmp$gcv
      ms_all_m[,,i]<-tmp$fd$coef
    }
    lam.ind = apply(gcv_all,1,which.min)
    lambdas<-lambda_all[lam.ind]
    xcoo = matrix(0,nbasis,n)
    for(i in 1:n){
      xcoo[,i]= ms_all_m[,i,lam.ind[i]]
    }
  }else if(p>1){
    lambda_all<- exp(seq(log(5*10^(-15)), log(5*10^(0)), len=ncv))
    nlam<-length(lambda_all)
    gcv_all<-array(0,c(n,p, nlam))
    ms_all_m<-array(0, c(nbasis, n,p,nlam))
    for(i in 1:nlam){
      lambda=lambda_all[i]
      fd_par<-fda::fdPar(fdobj=basis,Lfdobj=2,lambda=lambda)
      tmp<-fda::smooth.basis(tt,xraw,fdParobj=fd_par)
      gcv_all[,,i]<-tmp$gcv
      ms_all_m[,,,i]<-tmp$fd$coef
    }
    lam.ind = apply(gcv_all,c(1,2),which.min)
    
    lambdas<-matrix(0,n,p)
    xcoo = array(0,c(nbasis,n,p))
    for(i in 1:n){
      lambdas[i,] = lambda_all[lam.ind[i,]]
      xcoo[,i,]= ms_all_m[,i,,lam.ind[i]]
    }
  }
  
  return(list(coef=xcoo,lambdas=lambdas, basis=basis, lambdamin=max(lambdas[lam.ind])  ))
}



#' @export
# to test through after generating data
# basisno=5 ; X=X.obs ; lambda=NULL; Y=Y 
MFSGrp =function(Ytrain,Xtrain, basisno=5 ,tt, lambda=NULL, alpha=NULL ,
                part, Xpred=NULL,Ypred=NULL, Silence=FALSE, bspline=FALSE, Penalty=NULL, 
                lambdafactor=0.005, nfolds=5, 
                predloss="L2", eps = 1e-08, maxit = 3e+08, nlambda=100, forcezero=FALSE, 
                forcezeropar=0.001, sixplotnum=1, lambdaderivative=NULL,
                nfolder=5, nalpha=9, nlamder=10, lamdermin=1e-9, lamdermax=1e-3,alphamin=0 ,alphamax=1,
                a=3.7, ADMM=FALSE,numcores=NULL, rho=1 , unbalanced=FALSE){
  
  
  loss="ls";
  Y=Ytrain; X=Xtrain;
#######if( ) stop("orthaganlization must be one (orth=T) if ")
    if(bspline==T) {fpca=T} else {fpca=F}
  #Y=Ytrain;X=Xtrain;basisno=m; mm =m ;tt=tt; part=part; lambda=0.6;Xpred=Xtest;Ypred=Ytest; Silence=TRUE; fpca=T; bspline=T;  Penalty="glasso"; rho=1; lambdaderivative =NULL;alpha=0.5
  if (unbalanced==FALSE){
  p=dim(X)[1] #how many curves
  n=dim(X)[2] #how many random sample of each curve
  nt=dim(X)[3]} #number of grid points observed for each of n observed curve
  
  if (unbalanced==TRUE) {# X is a list
    p=length(lengths(X)); 
    n=dim(X[[1]])[1];
    nt=rep(0,p);
    for (j in 1:p) nt[j]=dim(X[[j]])[2]
  }
  
  K=length(part)# number of blocks
  m=basisno # number of baisis
  # next line allow user to use bspline instead of foriur bases if they chose bspline=true
  # deafualt = foriur
  if (bspline==TRUE) {B.basis=fda::create.bspline.basis(rangeval=c(0,1), nbasis=m)} else {B.basis=fda::create.fourier.basis(rangeval=c(0,1), nbasis=m)}
  # gram matrix by default is Foriur and so is I, unles bspline=TRUE
  if (bspline==TRUE){Gram=fda::bsplinepen(B.basis, 0)} else{Gram=fda::fourierpen(B.basis, 0)}
  # transforming  to coordinate representation wrt bases:
  oldGram=Gram
  #stop("here")
  #### if  fpca
  #lambdader=0
  if (fpca==TRUE){
    T= array(NaN, c(m,m,p));
    Xcoeforig=matrix(NaN, n,m*p);
    Xcoef=matrix(NaN, n,m*p);
    if (unbalanced==FALSE){
    for (j in 1:p){
      XX=get.fd(t(X[j,,]),tt, basis = B.basis)
      Xcoeforig[,((j-1)*m+1): (j*m)]=t(XX$coef)
      Xpca=fda::pca.fd( fda::fd(XX$coef, XX$basis)   , nharm=m)
      T[,,j]=  t(Xpca$harmonics$coefs)%*%Gram
      Xcoef[,((j-1)*m+1): (j*m)]= t(T[,,j]%*%t(Xcoeforig[,((j-1)*m+1): (j*m)]))
      #print(XX$lambdamin)
      #lambdader=max(lambdader,XX$lambdamin)
    }}
    if (unbalanced==TRUE){ # X and tt are lists now
      for (j in 1:p){
        XX=get.fd(t(X[[j]]),tt[[j]], basis = B.basis)
        Xcoeforig[,((j-1)*m+1): (j*m)]=t(XX$coef)
        Xpca=fda::pca.fd( fda::fd(XX$coef, XX$basis)   , nharm=m)
        T[,,j]=  t(Xpca$harmonics$coefs)%*%Gram
        Xcoef[,((j-1)*m+1): (j*m)]= t(T[,,j]%*%t(Xcoeforig[,((j-1)*m+1): (j*m)]))
        #print(XX$lambdamin)
        #lambdader=max(lambdader,XX$lambdamin)
      }}
    
    
    Gram=diag(m)
  } else {
    Xcoef=matrix(NaN, n,m*p);
    
    if (unbalanced==FALSE){ 
    for (j in 1:p){
      XX=get.fd(aperm(X[j,,]),tt, basis = B.basis)
      Xcoef[,((j-1)*m+1): (j*m)]=t(XX$coef)
      #lambdader=max(lambdader,XX$lambdamin)
    }}
    if (unbalanced==TRUE){ 
      for (j in 1:p){
        XX=get.fd(aperm(X[[j]]),tt[[j]], basis = B.basis)
        Xcoef[,((j-1)*m+1): (j*m)]=t(XX$coef)
        #lambdader=max(lambdader,XX$lambdamin)
      }}
    
    
  }



    lambdader=lambdaderivative; 

   


  ######
  # mean and centerizing X through subtracting coordinate represenation's mean
  meanX=colMeans(Xcoef)
  Xcoefc=t(t(Xcoef)-colMeans(Xcoef))
  
  
  #for (j in 1:(m*p)){
   # Xcoef[, j ]=Xcoef[, j ]/sd(Xcoef[, j ]) }

  
  # part is vector . Each element is number of covariate in each group
  # sum of all of them = number of all covariates
  cum_part = cumsum(part)
  # Centerize Y through subtracting mean
  Ymean=mean(Y)
  Y=Y-Ymean

  # in beta update this is used   #diagonal of blocks of Graham
  GG=Matrix::bdiag(replicate(p, Gram , simplify = FALSE))%*%diag(m*p)
  GG=as.matrix(GG)
  ## GG of second derivative penalty
  if (bspline==TRUE){Gramder=fda::bsplinepen(B.basis, 2)} else{Gramder=fda::fourierpen(B.basis, 2)}
  
  if (fpca==FALSE){
    GGder=Matrix::bdiag(replicate(p, Gramder , simplify = FALSE))%*%diag(m*p);
  }
  else {
    TT=list(  t(solvecpp(T[,,1]))%*%Gramder%*%solvecpp(T[,,1])   )
    for (j in 2:p){
      listemp=list(t(solvecpp(T[,,j]))%*%Gramder%*%solvecpp(T[,,j]))
      TT=append(TT , listemp)
    }
    GGder=Matrix::bdiag(TT)%*%diag(m*p)
  }
  GGder=as.matrix(GGder)
  
  
  
  
  if(ADMM==FALSE){
  #lamderfactor=1e5
#print(lamderfactor )

  #############################
  if (Penalty=="OLS"){
    ######################################
    
    
    
    if(is.null(lambdader))

{
      #nfolder=3
      #nlamder=10
      #if (bspline==TRUE) { lambdaders=exp(seq(log(1e-6), log(1e-3),len=nlamder))}
      #else {lambdaders=exp(seq(log(1e-8), log(1e-3),len=nlamder))}
      lambdaders=exp(seq(log(lamdermin), log(lamdermax),len=nlamder))
      #lambdaders=append(0, lambdaders)
      #nlamder=nlamder+1
      MSEall =matrix(0, nfolder, nlamder)
      
      folds = rep_len(1:nfolder, n)
      folds = sample(folds, n)
      inv.GG = solvecpp(GG)
      
      for(k in 1:nfolder){
        
        
        
        fold.ind = which(folds == k)
        
        x.train = Xcoefc[-fold.ind, ]
        x.test = Xcoefc[fold.ind, ]
        y.train = Y[-fold.ind]
        y.test = Y[fold.ind]
    
        #print(dim(x.train))
        #print(length(y.train))
        
      
        
        YXcoec=t(y.train)%*%x.train
        ridgecoef=t(x.train)%*%x.train%*%GG
        
        
        for(ilam in 1:nlamder){
          cat(k,"th fold of ", nfolder ," for the ", ilam, "th lambdader of", nlamder  , "\r" )
          ridgecoef=ridgecoef+lambdaders[ilam]*(inv.GG%*%GGder)
          rdgcoefinv=solvecpp(ridgecoef)
          olscoef=rdgcoefinv%*%t(YXcoec)
          pred = as.vector(x.test%*%GG%*%olscoef+Ymean)
          MSEall[k,ilam] = mean((pred-y.test)^2)  
        }
      }
      MSEs=apply(MSEall,2,mean)
      #print(MSEs)
      #print(lambdaders)
      #cat("\n")
      par(mfrow=c(1,1))
      #lambdader[1]=exp(-50)
      plot(y=MSEs, x=log(lambdaders) )
      lambdader=lambdaders[which.min(apply(MSEall,2,mean))]
      
      cat("\r  Chosen lambdader is", lambdader, "                             \r \n" )
      

      
      }
  
    
    #########################
    # last term of ols regression, admm beta update, and ridge is always X^t Y :
    YXcoec=Y%*%Xcoefc
    ridgecoef=t(Xcoefc)%*%Xcoefc%*%GG
    ridgecoef=ridgecoef+lambdader*(solvecpp(GG)%*%GGder)
    rdgcoefinv=solvecpp(ridgecoef)
    olscoef=rdgcoefinv%*%t(YXcoec)
    ols=olscoef
    beta=ols
  }

  else if  (Penalty== "glasso" |Penalty== "gelast" | Penalty== "ridge" ) {
    # find a lambda that almost work to achive sparsity
    
    if (Penalty=="glasso" ) {alpha=0;}
    if(Penalty=="ridge"){alpha=1;}
    
  

   #print(GGder[1:5,1:5]) 
     #lambdader=lambdader*p

    beta = rep(0,m*p);
    group <- rep(1:p,each=m)
    cv <- fGMD::cv.fGMD(Xcoefc, Y, group=group, loss=loss,
                     pred.loss=predloss ,nfolds=nfolds, intercept = F,
                     lambda.factor=lambdafactor , alpha=alpha, lambda=NULL, 
                     maxit=maxit, eps=eps, nlambda=nlambda, GGder=GGder, lambdader=lambdader,
                     nfolder=nfolder, nalpha=nalpha, nlamder=nlamder, 
                     lamdermin=lamdermin, lamdermax=lamdermax,alphamin=alphamin , alphamax=alphamax) #,  # lambda-lambda
    par(mfrow=c(1,1))
    plot(cv)
    lambda=cv$lambda.min
    
    beta = coef(cv$fGMD.fit, s = cv$lambda.min)[-1]
    
    
    if (forcezero==TRUE){
    start_ind = 1;
    for (i in 1:K){
      sel = start_ind:cum_part[i];
      betasel=beta[sel]
      if (  sqrt(t(betasel)%*%Gram%*%betasel) < forcezeropar) beta[sel]=0
      start_ind = cum_part[i] + 1;} }
    
    
    ols=beta
    olscoef=beta
    
  } else {beta=NULL; ols=NULL;olscoef=NULL}
  
  #OLS estimated coefs
  ################################

  
  
  
  }##if ADMM is flase end
else{
  
  
  if (.Platform$OS.type == "windows") {
    numcores = 1
  } else {
    coremax=parallel::detectCores()
    if(!is.null(numcores)) { if(numcores>coremax) numcores=coremax;}
    if(is.null(numcores))  { numcores=coremax; }
    
  }
  
  
  
  YXcoec=Y%*%Xcoefc
  maxit = maxit/3e+6
  
  if (Penalty=="glasso" |  Penalty=="gscad") {alpha=0;}
  if(Penalty=="ridge" | Penalty=="OLS"){alpha=1;}
  if(Penalty=="OLS"){lambda=0;}
  
  euc=TRUE; if ( bspline==TRUE & fpca==F ) {euc=FALSE};
  group <- rep(1:p,each=m)
  
  
  
  if (is.null(lambda)) {
    lambdas=0
    start_ind = 1
    ridge=t(YXcoec)
    for (i in 1:K){
      sel = start_ind:cum_part[i]
      lambdas[i] = normcpp(ridge[sel],Gram);
      start_ind = cum_part[i] + 1;}
    maximlam = max(lambdas)/5;
    lam=exp(seq(log(maximlam),log(lambdafactor*maximlam), length.out=nlamder))/length(Y)
  }
  else{lam=lambda}
  
  
  #print(lambda.factor)
  
  #print(lam)
  
  #print(lam)
  
  #nfolder=3
  if ( is.null(alpha) | is.null(lambdader) )
  {    
    
    if( !is.null(alpha) | !is.null(lambdader))   {par(mfrow=c(1,1))} else {par(mfrow=c(1,2))}
    #############################
    #print(lamdermax)
    #print(lamdermin)
    #print(nlamder)
    #lambdader=NULL
    if(is.null(lambdader))
    {
      if (is.null(alpha)) {alp=0;} else {alp=alpha}
      #lam=(1-alp)*lam
      lambdaders=exp(seq(log(lamdermin), log(lamdermax),len=nlamder))
      #print(lambdaders)
      lambdadersmse=rep(NA, nlamder)
      systimes=rep(NA, nlamder)
      
      for (j in seq(nlamder)) {
        start_time <- Sys.time()
        #pbmcapply::pbmclapply
        #parallel::mclapply
        #cat("\r The ", j, "th of the ", nlamder  ," lamdaders " )
        mse=pbmcapply::pbmclapply(lam,function(lambdas){foldcpp(Y=Y,X=Xcoef, basisno=basisno ,tt=0, lambda=lambdas, 
                                                   alpha=alp , part=part, rho=rho , 
                                                   Penalty=Penalty, GG=GG, lambdader =lambdaders[j], 
                                                   GGder=GGder,K=K, Gram=Gram,oldGram=oldGram,n=n,p=p,m=m,
                                                   kf=nfolder, euc=euc, Path=FALSE,
                                                   maxit=maxit, eps=eps, a=a, id=j, idmax=nlamder) }, mc.cores=numcores)
        mse=unlist(mse, recursive = F, use.names = T)
        lambdadersmse[j]=min(mse)
        #print(lambdadersmse)
        
        systimes[j]=Sys.time()-start_time
        
      }
      #print(lambdaders)
      #print(systimes)
      # lamders[1]=exp(-50)
      plot(y=lambdadersmse, x=log(lambdaders)) 
      chose= which(lambdadersmse==min(lambdadersmse) )[1]
      #print(chose[1])
      lambdader=min(lambdaders[chose])
      ET=nfolds*nlambda*systimes[chose]/(nfolder*nlamder)
      
      #lambdader=5e-4
      #cat(systimes[chose[1]])
      #cat(systimes[chose], "\n")
      cat("\r Chosen lambdader is", lambdader, "and Maximum Estimated Time:",  ET  ," seconds                       \r \n" )
      #stop("here")
      
    }
    #####
    
    if(is.null(alpha))
    {
      
      #nalpha=9
      alphasmse=rep(NA, nalpha)
      nalphaseq=nalpha+2
      alphas=seq(alphamin, alphamax, len=nalphaseq )[-c(1,nalphaseq)]
      #alphas=seq(0.1,0.9, len=nalpha )
      systimes=rep(NA, nalpha)
      
      #print(alphas)
      for (j in seq(nalpha)) 
      {          start_time <- Sys.time()
      #pbmcapply::pbmclapply
      #parallel::mclapply
      #cat("\r The ", j, "th of the ", nalpha  ," alphas " )
      
      mse=pbmcapply::pbmclapply(lam,function(lambdas){foldcpp(Y=Y,X=Xcoef, basisno=basisno ,tt=0, lambda=lambdas, 
                                                 alpha=alphas[j] , part=part, rho=rho , 
                                                 Penalty=Penalty, GG=GG, lambdader =lambdader, 
                                                 GGder=GGder,K=K, Gram=Gram,oldGram=oldGram,n=n,p=p,m=m,
                                                 kf=nfolder, euc=euc, Path=FALSE,
                                                 maxit=maxit, eps=eps, a=a,
                                                 id=j, idmax=nalpha, alphanet=TRUE)}, mc.cores=numcores)
      
      mse=unlist(mse, recursive = F, use.names = T)
      alphasmse[j]=min(mse)
      #print(lambdadersmse)
      
      systimes[j]=Sys.time()-start_time
      
      }
      
      #print(lambdadersmse)
      #print(lambdaders)
      #print(systimes)
      
      plot(y=alphasmse, x=alphas )
      chose= which(alphasmse==min(alphasmse) )[1]
      #print(chose[1])
      alpha=min(alphas[chose])
      #lambdader=5e-4
      #cat(systimes[chose[1]])
      ET=nfolds*nlambda*systimes[chose]/(nfolder*nalpha)
      #cat(systimes[chose], "\n")
      cat("\r Chosen alpha is ", alpha, " and Maximum Estimated Time:  ",  ET  ," seconds                       \r \n" )
      #stop("here")
      
    }
    ############################### 
    
    
    
    
    
  }  
  
  par(mfrow=c(1,1))
  
  
  if( Penalty!="OLS"){
    
    if(is.null(lambda)) {
      lambdas=exp(seq(log(lambdafactor*maximlam), log(maximlam), length.out=nlambda))/length(Y)} 
    else{lambdas=lambda}
    #print(lambdas)
    mse=pbmcapply::pbmclapply(lambdas,function(lambd){foldcpp(Y=Y,X=Xcoef, basisno=basisno ,tt=0, lambda=lambd, 
                                                   alpha=alpha , part=part, rho=rho , 
                                                   Penalty=Penalty, GG=GG, lambdader =lambdader, 
                                                   GGder=GGder,K=K, Gram=Gram,oldGram=oldGram,n=n,p=p,m=m,
                                                   kf=nfolds, euc=euc, Path=TRUE,
                                                   maxit=maxit, eps=eps, a=a)}, mc.cores=numcores)
    
    
    mses=unlist(mse, recursive = F, use.names = T)
    #mses=unlist(mses[names(mses) %in% "MSE"])
    plot(y=mses, x=log(lambdas) )
    
    chose= which(mses==min(mses) )
    lambda=max(lambdas[chose])} else{lambda=0}
  #print(lambda)
  
  
  final=foldcpp(Y=Y,X=Xcoef, basisno=basisno ,tt=0, lambda=lambda, 
                alpha=alpha , part=part, rho=rho , 
                Penalty=Penalty, GG=GG, lambdader =lambdader, 
                GGder=GGder,K=K, Gram=Gram,oldGram=oldGram,n=n,p=p,m=m,
                kf=1, euc=euc, Path=TRUE,
                maxit=maxit, eps=eps, a=a)  
  
  final=unlist(final, recursive = F, use.names = T)
  
  beta=final
  ols=beta
  olscoef=beta
  

  
  
}
  
  
  #OLS estimated coefs
  ################################
  ### if fpca
  if (fpca==TRUE) {
    betaold=beta
    for (j in 1:p)
    {
      beta[((j-1)*m+1): (j*m)]=solvecpp(T[,,j])%*%beta[((j-1)*m+1): (j*m)]
    }
  }
  
  ols=beta
  olscoef=beta
  
  

  # This sileice =FALSE is compatble with  first three curves in simulation study
  if(Silence==FALSE){

    if (sixplotnum=="max"){sixplotnum= max(  ceiling( sum(beta!=0)/(6*m) ) , 1) ;}
    par(mfrow=c(2,3))
    if (p==3) {par(mfrow=c(1,3))}
    ################################
    i=0
    for (j in 0:(p-1)){
      ind=(j*m+1):((j+1)*m)
      beta=ols[ind]
      
      if(sum(abs(beta))!=0 ){ 
        i=i+1; if(p>3 & ceiling(i/6)>sixplotnum) {break;}
        fbeta=fda::fd(beta, B.basis);plot(fbeta, xlab = paste0(" Time. This is the plot of the Coef. ",j+1));   
        a=paste0("b",j+1)
      if( exists(a) ) {
        if ( unbalanced==FALSE){lines(y= eval(parse(text=a)),x=(1:(5*nt))*1/(5*nt), col="green")}
        if ( unbalanced==TRUE){lines(y= eval(parse(text=a)),x=(1:(5*nt[j+1]))*1/(5*nt[j+1]), col="green")}
        
        }
      }
      
      if( p>3 & j>(p-4) & i<6){i=i+1;if(ceiling(i/6)>sixplotnum) {break;} ; fbeta=fda::fd(beta, B.basis);plot(fbeta, xlab = paste0(" Time. This is the plot of the Coef. ",j+1));}

        }
  }
  
  
  # given a data set predict and find mean square errors, only activates if Xpred and Ypred are not null
  Ypred=as.vector(Ypred)
  predict=0
  MSEpredict=0
  if(!is.null(Xpred)) {
    Gram=oldGram
    #gram diag block
    #if (dim(Xpred)[2]!=length(Ypred)) {print("error dim test")}
    if (unbalanced==FALSE){
    ptest=dim(Xpred)[1]
    ntest=dim(Xpred)[2]
    nttest=dim(Xpred)[3]
    GG=Matrix::bdiag(replicate(ptest, Gram , simplify = FALSE))%*%diag(m*p)
    GG=as.matrix(GG)
    Xcoeftest=matrix(NaN, ntest,m*ptest);
    for (j in 1:p){
      XX=get.fd(t(Xpred[j,,]),tt, basis = B.basis)
      Xcoeftest[,((j-1)*m+1): (j*m)]=t(XX$coef)
    }
    }
    
    if (unbalanced==TRUE){
      ptest=length(lengths(Xpred))
      ntest=dim(Xpred[[1]])[1]
      nttest=rep(0,p); for (j in 1:p) nttest[j]=dim(Xpred[[j]])[2]
      GG=Matrix::bdiag(replicate(ptest, Gram , simplify = FALSE))%*%diag(m*p)
      GG=as.matrix(GG)
      Xcoeftest=matrix(NaN, ntest,m*ptest);
      for (j in 1:p){
        XX=get.fd(t(Xpred[[j]]),tt[[j]], basis = B.basis)
        Xcoeftest[,((j-1)*m+1): (j*m)]=t(XX$coef)
      }
    }
    
    
    # in prediction do the oppiste of centerizing that was doine on original date set
    Xcoeftestc=t(t(Xcoeftest)-meanX)
    
    #for (j in 1:(m*p)){
    #  Xcoeftestc[, j ]=Xcoeftestc[, j ]/sd(Xcoeftestc[, j ]) }
    
    
    predict=as.vector(Xcoeftestc%*%GG%*%olscoef+Ymean)
    # error vector
    E=Ypred-predict
    # Sum square error
    SSE=t(E)%*%E
    # means square error
    MSEpredict=SSE/ntest
  }

  
  # have a list to ask for what needed
  return(list("coef"=olscoef,"predict"=predict, "MSEpredict"=MSEpredict, "lambda"=lambda))
  
}







