#'@title psi-Learning in individulized treatment rule
#'@author MingyangLiu <liux3941@umn.edu>
#'@description Given the tunning parameters return the psiLearning model to estimate the optimal ITR
#'@usage psiITR(X,A,R,w0=NULL,tau=0.1,kappa=0.1,maxit=100,tol=1e-4,kernel='linear',sigma=NULL,res=FALSE)
#'@param X \eqn{n} by \eqn{p} input matrix.
#'@param A a vector of n entries coded 1 and -1 for the treatment assignments.
#'@param R a vector of outcome variable, larger is more desirable.
#'@param w0 Inital estimate for the coefficients from  \code{\link{psi_Init}} or can be provided by the user.
#'@param tau tuning parameter for the loss function in psi-Learn
#'@param kappa tunning parameter to control the complexity of the decision function.
#'@param maxit maximum iterations
#'@param kernel kernel function for psi-Learning, can be \code{'linear'} or \code{'rbf'} (radial basis kernel), default is \code{'linear'}.When \code{'rbf'} is specified, one can specify the \code{sigma} parameter for rbf kernel.
#'@param sigma when using the rbf kernel, the bandwidth parameter  for 'rbf' kernel, default is 0.5.
#'@param tol tolerance error bound
#'@param res Whether to estimate the residual as the outcome for interaction effect, default is FALSE
#'@seealso \code{\link{psi_Init}}
#'@return It returns the estimated coefficients in the decision funcion and the fitted value
#'  \item{w}{the coefficent for the decision function, if in the linear case it is p-dimension and if in the rbf kernel case, it is n-dimension.}
#'  \item{bias}{the intercept in both the linear case and the kernel case.}
#'  \item{fit}{a vector of estimated values for \eqn{\hat{f(x)}} in training data, in the linear case it is \eqn{fit=bias+X*w}  and in the kernel case \eqn{fit=bias+K(X,X)w}.}
#'@export
#'@importFrom glmnet cv.glmnet
#'@examples
#'           n=100;p=5
#'           X=replicate(p,runif(n, min = -1, max = 1))
#'           A=2*rbinom(n, 1, 0.5)-1
#'           T=cbind(rep(1,n,1),X)%*%c(1,2,1,0.5,rep(0,1,p-3))
#'           T0=(cbind(rep(1,n,1),X)%*%c(0.54,-1.8,-1.8,rep(0,1,p-2)))*A
#'           R=as.vector(rnorm(n,mean=0,sd=1)+T+T0)
#'           w0.Linear=psi_Init(X,A,R,kernel='linear')
#'           psi_Linear<-psiITR(X,A,R,w0.Linear,tau=0.1,kappa=0.5,maxit=100,tol=1e-4,kernel='linear')
#'           w0.rbf=psi_Init(X,A,R,kernel='rbf')
#'           sigma=Sig_est(X,A)
#'           psi_rbf<-psiITR(X,A,R,w0.rbf,tau=0.1,kappa=0.1,maxit=100,tol=1e-4,kernel='rbf',sigma=sigma)


psiITR<-function(X,A,R,w0=NULL,tau=0.1,kappa=0.1,maxit=100,tol=1e-4,kernel='linear',sigma=NULL,res=FALSE){
  n=dim(X)[1]; p=dim(X)[2]
  pai=A*sum(A==1)/length(A)+(1-A)/2

  # if (abs(sum(R))>1e-7){
  #   cvfit=cv.glmnet(X,R,nfolds=5)
  #   co=as.matrix(predict(cvfit,s="lambda.min",type="coeff"))
  #   R=R-cbind(rep(1,n),X)%*%co
  # }
  #
  if (res == TRUE){
    R=resEst(X,R,pai)
  }

  wt=drop(R/pai)
  tt=0; m=1; tol_s=tol;err=1
  Cost=Inf

  if(is.null(w0)){
    if (kernel=='rbf'){
      if (is.null(sigma)){
        sigma=Sig_est(X,A)
      }
      w0_init=psi_Init(X,A,R,sigma=sigma,kernel='rbf')
    }else{
      w0_init=psi_Init(X,A,R)
    }
    w0=c(w0_init$bias,w0_init$w)
  }else{w0=c(w0$bias,w0$w)}

w_old=w0
if (kernel=='linear'){
  while(m<maxit & err>tol){

    w_new_old=w_old[2:length(w_old)]
    bias_new_old=w_old[1]
    err_s=2*tol_s+0.1
    tt=0
    temp=drop(X%*%w_new_old)
    u_old=A*(temp+bias_new_old)
    tkern=A*X
    dudv=tkern; dudb=A
    lims=(wt>0)
    V1_old=1/tau*colMeans((-(u_old<=0)*lims-(1-u_old>=0)*(1-lims))*abs(wt)*dudv)
    V2_old=1/tau*mean((-(u_old<=0)*lims-(1-u_old>=0)*(1-lims))*abs(wt)*dudb)

    while(tt<maxit & err_s>tol_s){
      t_tt=0.1/(1+tt)
      u_new_old=drop(X%*%w_new_old+bias_new_old)

      dphi1du=-(1-A*u_new_old>=0)
      dphi2du=-(-A*u_new_old>=0)


      g_oldv=kappa*w_new_old+1/(n*tau)*colSums(((dphi1du*lims+dphi2du*(1-lims))*abs(wt))*dudv)-V1_old
      g_oldb=1/(n*tau)*sum(((dphi1du*lims+dphi2du*(1-lims))*abs(wt))*dudb)-V2_old

      w_new_new=w_new_old-t_tt*g_oldv
      bias_new_new=bias_new_old-t_tt*g_oldb
      err_s=sum((c(bias_new_new,w_new_new)-c(bias_new_old,w_new_old))^2)
      w_new_old=w_new_new
      bias_new_old=bias_new_new
      tt=tt+1;
    }
    w_new=c(bias_new_new,w_new_new);
    err=sum((w_old-w_new)^2);
    w_old=w_new;
    m=m+1;

    cost=tlp_S(X,A,R,wt,w_old,tau=tau,kappa=kappa,kernel='linear')
    if(abs(cost-Cost)<tol){
      break;
    } else if (cost<Cost){

      Cost = cost
      w_opt=w_new
    }

  }

 fit=X%*%w_opt[2:length(w_opt)]+w_opt[1]
 return(list(w=w_opt[2:length(w_opt)],bias=w_opt[1],fit=fit))
}else if (kernel=='rbf'){

  rbf=kernlab::rbfdot(sigma=sigma)
  K=kernlab::kernelMatrix(rbf,X)
  if(is.null(sigma)){
    sigma=Sig_est(X,A)
  }
  while(m<maxit & err>tol){
    w_new_old=w_old[2:length(w_old)]
    bias_new_old=w_old[1]
    err_s=2*tol_s+0.1
    tt=0
    temp=drop(K%*%w_new_old)
    u_old=A*(temp+bias_new_old)
    tkern=A*K
    dudv=tkern; dudb=A
    lims=(wt>0)
    V1_old=1/tau*colMeans((-(u_old<=0)*lims-(1-u_old>=0)*(1-lims))*abs(wt)*dudv)
    V2_old=1/tau*mean((-(u_old<=0)*lims-(1-u_old>=0)*(1-lims))*abs(wt)*dudb)

    while(tt<maxit & err_s>tol_s){
      t_tt=0.1/(1+tt)
      u_new_old=drop(K%*%w_new_old+bias_new_old)


      dphi1du=-(1-A*u_new_old>=0)
      dphi2du=-(-A*u_new_old>=0)


      g_oldv=kappa*(K%*%w_new_old)+1/(n*tau)*colSums(((dphi1du*lims+dphi2du*(1-lims))*abs(wt))*dudv)-V1_old
      g_oldb=1/(n*tau)*sum(((dphi1du*lims+dphi2du*(1-lims))*abs(wt))*dudb)-V2_old

      w_new_new=w_new_old-t_tt*g_oldv
      bias_new_new=bias_new_old-t_tt*g_oldb
      err_s=sum((c(bias_new_new,w_new_new)-c(bias_new_old,w_new_old))^2)
      w_new_old=w_new_new
      bias_new_old=bias_new_new
      tt=tt+1;
    }
    w_new=c(bias_new_new,w_new_new);
    err=sum((w_old-w_new)^2)
    w_old=w_new;
    m=m+1;


    cost=tlp_S(K,A,R,wt,w_old,tau=tau,kappa=kappa,kernel='rbf')
    if(norm(cost-Cost)<tol){
      break;
    } else if (cost<Cost){

      Cost = cost
      w_opt=w_new
    }

  }
  fit=K%*%w_opt[2:length(w_opt)]+w_opt[1]
  return(list(w=w_opt[2:length(w_opt)],bias=w_opt[1],fit=fit))
}
  }


