
#'@title Cross validation for psiLearn learning with variable selection
#'@author MingyangLiu <liux3941@umn.edu>
#'@description Return the psi-learning models with best tuning parameters.
#'@usage  cv.psiITRVS(X,A,R,m=5,kernel='linear',kappa.ratio=0.01,kappa.max=1,nkappa=10,tau=0.1,lambda.ratio=0.1,lambda.max=2,nlam=10,tau2=0.2,maxit=100,tol=1e-5,res=FALSE)
#'@param X \eqn{n} by \eqn{p} input matrix.
#'@param A a vector of n entries coded 1 and -1 for the treatment assignments.
#'@param R a vector of outcome variable, larger is more desirable.
#'@param m m-folds cross validation
#'@param kernel the kernel used in the decision, here we only use the linear kernel.
#'@param tau tuning parameter for the ridge penalty in psi-Learning.
#'@param kappa.ratio the ratio between the max kappa and the min kappa which  controls the complexity of the decision function.
#'@param kappa.max max kappa in the tunning parameter seq.
#'@param nkappa num of kappa for grid search.
#'@param lambda.ratio the ratio between the max lambda and the min lambda which  controls the complexity of the decision function in the TLP penalty.
#'@param lambda.max max lambda in the tunning parameter seq
#'@param nlam num of kappa for grid search
#'@param tau2 tunning parameter to control the margin used in TLP penalty,default is 0.2
#'@param maxit number of max iteration used in \code{\link{psiITR}}
#'@param tol tolerance used in \code{\link{psiITR}}
#'@param res Whether to estimate the residual as the outcome for interaction effect, default is FALSE
#'@return It returns the estimated coefficients in the decision funcion after cross validation
#'  \item{w}{the coefficent for the decision function.}
#'  \item{bias}{the intercept in both the linear case and the kernel case.}
#'  \item{sigma}{if the kernel is rbf then the optimal sigma is returned }
#'@export
#'@import foreach
#'@importFrom parallel detectCores
#'@examples
#'           n=100;p=5
#'           X=replicate(p,runif(n, min = -1, max = 1))
#'           A=2*rbinom(n, 1, 0.5)-1
#'           T=cbind(rep(1,n,1),X)%*%c(1,2,1,0.5,rep(0,1,p-3))
#'           T0=(cbind(rep(1,n,1),X)%*%c(0.54,-1.8,-1.8,rep(0,1,p-2)))*A
#'           R=as.vector(rnorm(n,mean=0,sd=1)+T+T0)
#'           cv_psi_LinearVS<-cv.psiITRVS(X,A,R,m=5,kernel='linear',kappa.ratio=0.01,kappa.max=1,nkappa=10,tau=0.1,lambda.ratio=0.1,lambda.max=2,nlam=10,tau2=0.2,maxit=100,tol=1e-5,res=FALSE)


cv.psiITRVS<-function(X,A,R,m=5,kernel='linear',kappa.ratio=0.01,kappa.max=1,nkappa=10,tau=0.1,lambda.ratio=0.1,lambda.max=2,nlam=10,tau2=0.2,maxit=100,tol=1e-5,res=FALSE){

  eps0=1e-3
  n=length(A)
  #X has no constant column
  # X=cbind(X,rep(1,dim(X)[1]))
  p=dim(X)[2]
  pai=A*sum(A==1)/length(A)+(1-A)/2
  wt=R/pai
  kappas <- exp(seq(log(kappa.max),log(kappa.ratio*kappa.max),len=nkappa))
  lambdas<-exp(seq(log(lambda.max),log(lambda.ratio*lambda.max),len=nlam))
  if (res == TRUE){
    r=resEst(X,R,pai)
  }else {r=R}
  wtr=r/pai

  if(kernel=='rbf'){
    stop(gettextf("in this variable selection, only linear kernel is accepted "))
  }

    V=foreach(fold=cv,.combine='cbind')%dopar% {
      Xtrain=X[fold,]
      rtrain=r[fold]
      Atrain=A[fold]
      Xtest=X[-fold,]
      rtest=r[-fold]
      Atest=A[-fold]
      paitest=pai[-fold]
      wtrtest=rtest/paitest
      w0=psi_Init(Xtrain,Atrain,rtrain,kernel='linear')
      res=matrix(0,nrow = nkappa,ncol=nlam)
      for (j in 1:nkappa){
        for (s in 1:nlam){
          fitobj<-psiITR_VS(Xtrain,Atrain,w0=w0,rtrain,tau=tau,kappa=kappas[j],lambda=lambdas[s],maxit=maxit, tol=tol,tau2=tau2)
          Apre=sign(Xtest%*%fitobj$w+fitobj$bias)
          res[j,s]=sum(wtrtest*(Atest==Apre))/sum(1/paitest*(Atest==Apre))

        }}
    }
    mimi=apply(V,c(1,2),mean)
    best=which(mimi==max(mimi),arr.ind=TRUE)
    best_kappa=kappas[best[1]]
    best_lambda=lambdas[best[2]]
    fitobj_best<-psiITR_VS(X,A,R,tau=tau,kappa=best_kappa,lambda=best_lambda,maxit=maxit, tol=tol,tau2=tau2)
    fitobj_best$kappa=best_kappa
    fitobj_best$lambda=best_lambda

    return(fitobj_best)

}
