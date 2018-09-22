
#'@title psi-Learning in individulized treatment rule in the linear case with variable selection
#'@author MingyangLiu <liux3941@umn.edu>
#'@description Given the tunning parameters return the psiLearning model to estimate the optimal ITR with variable selection
#'@usage psiITR_VS(X,A,R,w0=NULL,tau=0.1,kappa=0.2,lambda=0.8,maxit=100,tol=1e-4,tau2=0.2,res=FALSE)
#'@param X \eqn{n} by \eqn{p} input matrix.
#'@param A a vector of n entries coded 1 and -1 for the treatment assignments.
#'@param R a vector of outcome variable, larger is more desirable.
#'@param w0 Inital estimate for the coefficients from  \code{\link{psi_Init}} or can be provided by the user.
#'@param tau tuning parameter for the loss function in psi-Learn
#'@param kappa tunning parameter to control the complexity of the decision function in the ridge penaly
#'@param lambda tunning parameter to control the complexity of the decision function in the TLP penalty
#'@param maxit maximum iterations allowed
#'@param tol tolerance error bound
#'@param tau2 tunning parameter to control the margin used in TLP penalty
#'@param res  Whether to estimate the residual as the outcome for interaction effect, default is FALSE
#'@seealso \code{\link{psi_Init}}
#'@return It returns the estimated coefficients in the decision funcion and the fitted value
#'  \item{w}{the coefficent for the decision function, if in the linear case it is p-dimension and if in the rbf kernel case, it is n-dimension.}
#'  \item{bias}{the intercept in both the linear case and the kernel case.}
#'  \item{fit}{a vector of estimated values for \eqn{\hat{f(x)}} in training data, in the linear case it is \eqn{fit=bias+X*w}  and in the kernel case \eqn{fit=bias+K(X,X)w}.}
#'@export
#'@examples
#'           n=100;p=5
#'           X=replicate(p,runif(n, min = -1, max = 1))
#'           A=2*rbinom(n, 1, 0.5)-1
#'           T=cbind(rep(1,n,1),X)%*%c(1,2,1,0.5,rep(0,1,p-3))
#'           T0=(cbind(rep(1,n,1),X)%*%c(0.54,-1.8,-1.8,rep(0,1,p-2)))*A
#'           R=as.vector(rnorm(n,mean=0,sd=1)+T+T0)
#'           w0.Linear=psi_Init(X,A,R,kernel='linear')
#'           psi_Linear<-psiITR_VS(X,A,R,w0.Linear,tau=0.1,kappa=0.3,lambda=1,maxit=100,tol=1e-4,tau2=0.2)



psiITR_VS<-function(X,A,R,w0=NULL,tau=0.1,kappa=0.2,lambda=0.8,maxit=100,tol=1e-4,tau2=0.2,res=FALSE){
  n=dim(X)[1]; p=dim(X)[2]
  pai=A*sum(A==1)/length(A)+(1-A)/2

  tt=0; m=1; tol_s=tol
  Cost=Inf

  if (res == TRUE){
    R=resEst(X,R,pai)
  }

  wt=drop(R/pai)
  if(is.null(w0)){
    w0_init=psi_Init(X,A,R)
    w0=c(w0_init$bias,w0_init$w)
  }else{w0=c(w0$bias,w0$w)}
w_old=w0
outIt=0
while (outIt<5){
  err=1
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
    lims_wold=(abs(w_new_old)<=tau2)
    V1_old=1/tau*colMeans((-(u_old<=0)*lims-(1-u_old>=0)*(1-lims))*abs(wt)*dudv)
    V2_old=1/tau*mean((-(u_old<=0)*lims-(1-u_old>=0)*(1-lims))*abs(wt)*dudb)
    V3_old=lambda/tau2*lims_wold*sign(w_new_old)


     while(tt<maxit & err_s>tol_s){
      t_tt=0.1/(1+tt)
      # t_tt = 1/(2*sqrt(n*(1+tt)));
      u_new_old=drop(X%*%w_new_old+bias_new_old)

      dphi1du=-(1-A*u_new_old>=0)
      dphi2du=-(-A*u_new_old>=0)


      g_oldv=kappa*w_new_old+1/(n*tau)*colSums(((dphi1du*lims+dphi2du*(1-lims))*abs(wt))*dudv)-V1_old+V3_old
      g_oldb=1/(n*tau)*sum(((dphi1du*lims+dphi2du*(1-lims))*abs(wt))*dudb)-V2_old

      w_new_new=w_new_old-t_tt*g_oldv
      bias_new_new=bias_new_old-t_tt*g_oldb
      err_s=sum((c(bias_new_new,w_new_new)-c(bias_new_old,w_new_old))^2)
      w_new_old=w_new_new
      bias_new_old=bias_new_new
      tt=tt+1;
    }
    w_new=c(bias_new_new,w_new_new)
    err=sum((w_old-w_new)^2)
    w_old=w_new
    m=m+1
  }

    cost=tlp_S(X,A,R,wt,w_new,tau=tau,kappa=kappa)
   # cost=tlp_S_VS(X,A,R,wt,w_old,tau=tau,kappa=kappa,lambda=lambda,kernel='linear')
    if(abs(cost-Cost)<1e-7){
      break;
    } else if (cost<Cost){

      Cost = cost
      w_opt=w_new
    }


}
  fit=X%*%w_opt[2:length(w_opt)]+w_opt[1]
  return(list(w=w_opt[2:length(w_opt)],bias=w_opt[1],fit=fit))
}

