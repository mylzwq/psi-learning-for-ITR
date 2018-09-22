#'@title Cross validation for psiLearn learning
#'@author MingyangLiu <liux3941@umn.edu>
#'@description Return the psi-learning models with best tuning parameters.
#'@usage  cv.psiITR(X,A,R,m=5,kernel='linear',sigma=NULL,kappa.ratio=0.01,kappa.max=1.5,nkappa=10,tau=0.1,maxit=100,tol=1e-5,res=FALSE)
#'@param X \eqn{n} by \eqn{p} input matrix.
#'@param A a vector of n entries coded 1 and -1 for the treatment assignments.
#'@param R a vector of outcome variable, larger is more desirable.
#'@param m m-folds cross validation
#'@param kernel kernel function used in the decision function
#'@param sigma  bindwidth for 'rbf' kernel it can be provided by the user, if not, it can be estimated from \code{\link{Sig_est}}
#'@param kappa.ratio the ratio between the max kappa and the min kappa which  controls the complexity of the decision function.
#'@param kappa.max max kappa in the tunning parameter seq
#'@param nkappa num of kappa for grid search
#'@param tau tuning parameter for the loss function in psi-Learn
#'@param maxit number of max iteration used in \code{\link{psiITR}}
#'@param tol tolerance used in \code{\link{psiITR}}
#'@param res Whether to estimate the residual as the outcome for interaction effect, default is FALSE
#'@return It returns the estimated coefficients in the decision funcion after cross validation
#'  \item{w}{the coefficent for the decision function.}
#'  \item{bias}{the intercept in both the linear case and the kernel case.}
#'  \item{sigma}{if the kernel is rbf then the optimal sigma is returned }
#'@export
#'@import foreach
#'@import caret
#'@importFrom parallel detectCores
#'@examples
#'           n=100;p=5
#'           X=replicate(p,runif(n, min = -1, max = 1))
#'           A=2*rbinom(n, 1, 0.5)-1
#'           T=cbind(rep(1,n,1),X)%*%c(1,2,1,0.5,rep(0,1,p-3))
#'           T0=(cbind(rep(1,n,1),X)%*%c(0.54,-1.8,-1.8,rep(0,1,p-2)))*A
#'           R=as.vector(rnorm(n,mean=0,sd=1)+T+T0)
#'           cv_psi_Linear<-cv.psiITR(X,A,R,m=5,kernel='linear',kappa.ratio=0.01,kappa.max=1.5,nkappa=10,tau=0.1,maxit=100, tol=1e-5)
#'           cv_psi_Rbf<-cv.psiITR(X,A,R,m=5,kernel='rbf',kappa.ratio=0.01,kappa.max=1,nkappa=10,tau=0.1,maxit=100, tol=1e-5)



cv.psiITR<-function(X,A,R,m=5,kernel='linear',sigma=NULL,kappa.ratio=0.01,kappa.max=1.5,nkappa=10,tau=0.1,maxit=100, tol=1e-5,res=FALSE){

  eps0=1e-3
  n=length(A)
  #X has no constant column
  # X=cbind(X,rep(1,dim(X)[1]))
  p=dim(X)[2]
  pai=A*sum(A==1)/length(A)+(1-A)/2
  wt=R/pai
  kappas <- exp(seq(log(kappa.max),log(kappa.ratio*kappa.max),len=nkappa))
  if (res == TRUE){
    r=resEst(X,R,pai)
  }
  else r=R

  wtr=r/pai

  #rand=sample(m,n,replace=TRUE)

  num_cores <- parallel::detectCores() - 1
  cl<-makeCluster(num_cores)
  registerDoParallel(cl)

  if (kernel=='linear'){
      cv <- caret::createFolds(A, k=m, returnTrain = TRUE)
      V=foreach(fold=cv,.combine='cbind')%dopar% {
      Xtrain=X[fold,]
      rtrain=r[fold]
      Atrain=A[fold]
      Xtest=X[-fold,]
      rtest=r[-fold]
      Atest=A[-fold]
      paitest=pai[-fold]
      wtrtest=rtest/paitest
      w0=psi_Init(Xtrain,Atrain,rtrain,kernel=kernel)
      res=rep(0,nkappa)
      for (j in 1:nkappa){
        fitobj<-psiITR(Xtrain,Atrain,rtrain,w0,tau=tau,kappa=kappas[j])
        #kern_test=kernelFunc(Xtest,Xtrain,kern=kernel,param=kParam)
        Apre=sign(Xtest%*%fitobj$w+fitobj$bias)
        res[j]=sum(wtrtest*(Atest==Apre))/sum(1/paitest*(Atest==Apre))
      }
    }
    mini=colMeans(V)
    best=which.max(mini)
    kappa_best=kappas[best]
    fitobj_best<-psiITR(X,A,R,tau=tau,kappa=kappa_best,maxit=maxit, tol=tol,kernel='linear')
    fitobj_best$kappa=kappa_best
    return(fitobj_best)
  }


  if (kernel=='rbf'){
    if (is.null(sigma)){
      sigma=Sig_est(X,A)
    }
    cv <- caret::createFolds(A, k=m, returnTrain = TRUE)
    V=foreach(fold=cv,.combine='cbind')%dopar% {
      Xtrain=X[fold,]
      rtrain=r[fold]
      Atrain=A[fold]
      Xtest=X[-fold,]
      rtest=r[-fold]
      Atest=A[-fold]
      paitest=pai[-fold]
      wtrtest=rtest/paitest
      w0=psi_Init(Xtrain,Atrain,Rtrain,kernel='rbf')
      res=rep(0,nkappa)
      for (j in 1:nkappa){
        rbf=kernlab::rbfdot(sigma=sigma)
        fitobj<-psiITR(Xtrain,Atrain,Rtrain,w0,tau=tau,kappa=kappas[j],maxit=maxit, tol=tol,kernel='rbf',sigma=sigma)
        Ktest=kernlab::kernelMatrix(rbf,Xtest,Xtrain)
        Apre=sign(Ktest%*%fitobj$w+fitobj$bias)
        res[j]=sum(wtrtest*(Atest==Apre))/sum(1/paitest*(Atest==Apre))
      }
    }

    mini=colMeans(V)
    best=which.max(mini)
    kappa_best=kappas[best]
    fitobj_best<-psiITR(X,A,R,tau=tau,kappa=kappa_best,maxit=maxit,tol=tol,kernel='rbf',sigma=sigma)
    fitobj_best$kappa=kappa_best
    return(fitobj_best)
  }
}

