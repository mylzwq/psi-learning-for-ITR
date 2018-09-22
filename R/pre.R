#'@title predict treatment in psi-Learning
#'@description This function predicts the optimal treatments with model of psi-Learning,which is estimated from \code{\link{psiITR}} or \code{\link{psiITR_VS}}
#'@usage pre.psi(fitobj,Xnew,X=NULL,kernel='linear')
#'@param fitobj model of class from psi-Learning
#'@param Xnew new data for prediction, mainly \eqn{n} by \eqn{p} input matrix.
#'@param X in linear case this is NULL and in 'rbf' kernel case this is the training data
#'@param kernel kernel function for psi-Learning, can be \code{'linear'} or \code{'rbf'} (radial basis kernel), default is \code{'linear'}.When using \code{'rbf'} , the bandwidth parameter
#'            sigma is asked to provide
#'@return predicted treatment for each new observation
#'@export
#'@importFrom kernlab rbfdot kernelMatrix
#'@examples
#'           n=100;p=5;ntest=200
#'           X=replicate(p,runif(n, min = -1, max = 1))
#'           A=2*rbinom(n, 1, 0.5)-1
#'           T=cbind(rep(1,n,1),X)%*%c(1,2,1,0.5,rep(0,1,p-3))
#'           T0=(cbind(rep(1,n,1),X)%*%c(0.54,-1.8,-1.8,rep(0,1,p-2)))*A
#'           R=as.vector(rnorm(n,mean=0,sd=1)+T+T0)
#'           cv_psi_Linear<-cv.psiITR(X,A,R,m=5,kernel='linear',kappa.ratio=0.01,kappa.max=1.5,nkappa=10,tau=0.1,maxit=100, tol=1e-5)
#'           cv_psi_Rbf<-cv.psiITR(X,A,R,m=5,kernel='rbf',kappa.ratio=0.01,kappa.max=1,nkappa=10,tau=0.1,maxit=100, tol=1e-5)
#'           Xtest=replicate(p,runif(ntest, min = -1, max = 1))
#'           pre.psiLin=pre.psi(cv_psi_Linear,Xtest)
#'           sig=Sig_est(X,A)
#'           pre.psiRbf=pre.psi(cv_psi_Rbf,Xtest,X,kernel='rbf',sigma=sig)

pre.psi<-function(fitobj,Xnew,X=NULL,kernel='linear',sigma=NULL){
  if (kernel=='linear'){
    pre=sign(Xnew%*%fitobj$w+fitobj$bias)
  }else if (kernel=='rbf'){
    if (is.null(sigma)){
      stop(gettextf("in the rbf kernel case the bandwith parameter is asked to provide"))
    }
    rbf=kernlab::rbfdot(sigma=sigma)
    K=kernlab::kernelMatrix(rbf,Xnew,X)
    pre=sign(K%*%fitobj$w+fitobj$bias)
  }
  return(pre)
}


