// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//


// [[Rcpp::export]]
arma::mat solvecpp(const arma::mat& A){
    return(pinv(A, 1.0e-25));}




// [[Rcpp::export]]
double maxeigencpp( const arma::mat& X) {
    return(max(eig_sym( X )));
}


// [[Rcpp::export]]
double normcpp(const arma::vec& x, const arma::mat& A, bool euc=true) {
    if (euc == true)  {return(norm(x,2));}
    else {return(as_scalar(sqrt(trans(x) *A* x) ));}
}





// [[Rcpp::export]]
arma::rowvec softcpp(const arma::vec& x, const arma::mat& A, double lambda, bool euc=true) {
  if (euc==true){ double b;
    if(normcpp(A*x,A,euc) <= lambda) b=0; 
    else {b=1-lambda/normcpp(x,A, euc); }   return(trans(b*x));}
  else{double b;
    b=1-lambda/normcpp(x,A, euc);
    b=std::max(b,0.0);
    if (normcpp(A*x,A,euc) <= lambda) {
      b=0;}
    else {b=b;} return(trans(b*x));}
}




// [[Rcpp::export]]
arma::rowvec  softscadcpp(const arma::vec& x, const arma::mat& Gram, double lambda, double rho, double a=3.7, bool euc=true) {
  double normvec=normcpp(x,Gram, euc);
  arma::rowvec ans;
  if   (normvec <= lambda* (1+1/rho) ) { ans=softcpp(x,Gram,lambda/rho, euc);}
  if   (  (lambda*(1+1/rho) < normvec) & (normvec <= a*lambda)  ) { double cons=(a-1)*rho-a*lambda/normvec; double cons2=cons/(rho*a-rho-1) ; ans=cons2*x.t();}
  if   (normvec > a*lambda) { ans=x.t();}
  return ans;}




// [[Rcpp::export]]
arma::vec foldcpp(arma::vec Y, const arma::mat& X, int basisno, const arma::rowvec& tt , double lambda, 
                  double alpha, const arma::mat& GG, double lambdader , 
                  const arma::rowvec& part, double rho ,
                  std::string Penalty, const arma::mat& GGder,
                  int K, const arma::mat& Gram, const arma::mat& oldGram, int n, int p, int m
                    , int kf=10, bool euc=true, bool Path=false, double eps=1e-08,
                      int maxit=1000, double a=3.7, int id=1, int idmax=1,
                      bool alphanet=false){

  arma::mat  Xmain=X;
  arma::vec  Ymain=Y;
  int step=floor(n/kf);
  double MSE=0;
  arma::vec Ypred;
  arma::mat Xcoef;
  arma::mat Xcoefc;
  arma::mat Xpred;
  arma::rowvec meanX;
  arma::rowvec cum_part;
  double Ymean;
  arma::rowvec YXcoec;
  arma::mat ridgecoef;
  arma::mat rdgcoefinv;
  arma::vec olscoef;
  arma::vec ols;
  arma::vec beta(m*p);
  arma::vec output;
  double ABSTOL=eps*1e+4;
  double RELTOL=eps*1e+5;
  
  
  if(kf>1){
    
    for(int l=1; l<=kf;l++){ 
      
      int start=(l-1)*step;
      int end=l*step-1;
      arma::uvec testind=arma::regspace<arma::uvec>(start , 1 ,end);
      Y=Ymain;
      Y.shed_rows(testind);
      Xcoef=Xmain;
      Xcoef.shed_rows(testind);
      Ypred=Ymain.rows(testind);
      Xpred=Xmain.rows(testind);
      
      
      meanX=mean(Xcoef,0);
      Xcoefc=Xcoef;
      Xcoefc.each_row()-= meanX; 
      
      cum_part = cumsum(part);
      Ymean=mean(Y);
      Y-=Ymean;
      
      YXcoec=Y.t()*Xcoefc;
      
      if (Penalty=="OLS"|| Penalty=="ridge"){rho=0; }
      ridgecoef=Xcoefc.t()*Xcoefc*GG+rho*arma::eye(m*p,m*p);
      ridgecoef+=lambdader*(inv(GG)*GGder);
      rdgcoefinv=inv(ridgecoef);
      
      
      
      
      if (Penalty=="OLS"){
        olscoef=rdgcoefinv*YXcoec.t();
        ols=olscoef;  
        beta=ols;
      }  else if(Penalty=="ridge") {
        ridgecoef+=lambda*arma::eye(m*p,m*p);
        rdgcoefinv=inv(ridgecoef);
        olscoef=rdgcoefinv*YXcoec.t();
        ols=olscoef;
        beta=ols;
      } else if  (Penalty== "glasso" ||Penalty== "gelast" || Penalty== "gscad") {
        
        
        
        
        arma::vec gamma=arma::zeros<arma::vec>(m*p);
        arma::vec u=arma::zeros<arma::vec>(m*p);
        beta.zeros();
        
        
        double softcoef=rho/(rho+2*alpha*lambda);
        double softpar=lambda*(1-alpha)/rho;
        
        
        
        for(int k=1; k<maxit+1;k++)
        {
          arma::vec q=YXcoec.t() + rho*(gamma - u);
          beta=rdgcoefinv*q;
          arma::vec gammaold=gamma; 
          int start_ind = 0;  
          
          for(int i=0; i<K;i++)
          { 
            arma::uvec sel =arma::regspace<arma::uvec>(start_ind, 1, (cum_part(i)-1) );
            if (Penalty=="glasso") {gamma.elem(sel) = softcpp( beta.elem(sel) + u.elem(sel), Gram,lambda/rho);}
            if (Penalty=="gelast") {gamma(sel) = softcoef*softcpp(beta(sel) + u(sel), Gram,softpar);}
            if (Penalty=="gscad") {gamma(sel) = softscadcpp(beta(sel) + u(sel), Gram,lambda, rho, a=a);}
            start_ind = cum_part(i) ;
            
          }
          
          
          u += (beta - gamma);
          arma::vec r=beta-gamma;
          arma::vec s=rho*GG*(gammaold-gamma);
          double epsprim=sqrt(p)*RELTOL+ ABSTOL* std::max(normcpp(beta,GG, euc),normcpp(gamma,GG, euc));
          double epsdual= sqrt(p)*RELTOL+ (ABSTOL/rho)*normcpp(u,GG, euc); 
          
          if (  (normcpp(r,GG, euc) <= epsprim) & (normcpp(s,GG, euc) <=  epsdual) ) {break;}
          
          
          
          
        }
        
        
        beta = gamma; // test
        ols=beta;
        olscoef=beta;
        
      }
      ols=beta;
      olscoef=beta;
      
      arma::vec predict;
      double MSEpredict;
      int ntest=testind.size();
      
      Xpred.each_row()-=meanX;
      predict=vectorise(Xpred*GG*olscoef, 0);
      predict+=arma::as_scalar(Ymean);
      arma::vec E=Ypred-predict;
      double SSE=norm(E,2);
      SSE=SSE*SSE;
      MSEpredict=SSE/ntest;
      MSE+=MSEpredict;

    
    
    if(Path==false){ 
      if(alphanet==false){Rcpp::Rcout <<'\r' << "     The    " <<id<< "th   of   the   " <<idmax << "   lamdaderivatives "  ;}
      else{Rcpp::Rcout <<'\r' << "     The    " <<id<< "th   of   the   " <<idmax << "   alphas " ;} }
    else { Rcpp::Rcout <<'\r' << "        The      final   regularization    for         lambda    " ;}
    
    }
    
    
    MSE/=kf;
    
    output=MSE;
  }
  
  if (kf==1){
    
    
    
    
    
    
    Y=Ymain;
    Xcoef=Xmain;
    
    
    meanX=mean(Xcoef,0);
    Xcoefc=Xcoef;
    Xcoefc.each_row()-= meanX; 
    
    cum_part = cumsum(part);
    Ymean=mean(Y);
    Y-=Ymean;
    
    YXcoec=Y.t()*Xcoefc;
    
    if (Penalty=="OLS"|| Penalty=="ridge"){rho=0; }
    ridgecoef=Xcoefc.t()*Xcoefc*GG+rho*arma::eye(m*p,m*p);
    ridgecoef+=lambdader*(inv(GG)*GGder);
    rdgcoefinv=inv(ridgecoef);
    
    
    
    
    if (Penalty=="OLS"){
      olscoef=rdgcoefinv*YXcoec.t();
      ols=olscoef;  
      beta=ols;
    }  else if(Penalty=="ridge") {
      ridgecoef+=lambda*arma::eye(m*p,m*p);
      rdgcoefinv=inv(ridgecoef);
      olscoef=rdgcoefinv*YXcoec.t();
      ols=olscoef;
      beta=ols;
    } else if  (Penalty== "glasso" ||Penalty== "gelast" || Penalty== "gscad") {
      
      
      
      arma::vec gamma=arma::zeros<arma::vec>(m*p);
      arma::vec u=arma::zeros<arma::vec>(m*p);
      beta.zeros();
      
      
      double softcoef=rho/(rho+2*alpha*lambda);
      double softpar=lambda*(1-alpha)/rho;
      
      
      
      for(int k=1; k<maxit+1;k++)
      {
        arma::vec q=YXcoec.t() + rho*(gamma - u);
        beta=rdgcoefinv*q;
        arma::vec gammaold=gamma; 
        int start_ind = 0;  
        
        for(int i=0; i<K;i++)
        { 
          arma::uvec sel =arma::regspace<arma::uvec>(start_ind, 1, (cum_part(i)-1) );
          if (Penalty=="glasso") {gamma.elem(sel) = softcpp( beta.elem(sel) + u.elem(sel), Gram,lambda/rho);}
          if (Penalty=="gelast") {gamma(sel) = softcoef*softcpp(beta(sel) + u(sel), Gram,softpar);}
          if (Penalty=="gscad") {gamma(sel) = softscadcpp(beta(sel) + u(sel), Gram,lambda, rho, 3);}
          start_ind = cum_part(i) ;
          
        }
        
        
        u += (beta - gamma);
        arma::vec r=beta-gamma;
        arma::vec s=rho*GG*(gammaold-gamma);
        double epsprim=sqrt(p)*RELTOL+ ABSTOL* std::max(normcpp(beta,GG, euc),normcpp(gamma,GG, euc));
        double epsdual= sqrt(p)*RELTOL+ (ABSTOL/rho)*normcpp(u,GG, euc); 
        
        if (  (normcpp(r,GG, euc) <= epsprim) & (normcpp(s,GG, euc) <=  epsdual) ) {break;}
        
        
        
        
      }
      
      
      beta = gamma; // test
      ols=beta;
      olscoef=beta;
      
    }
    ols=beta;
    olscoef=beta;
    
    
    

      
    
    
    
    
    
    output=beta;
    
    
  }
  
  return(output);
  
  
}


