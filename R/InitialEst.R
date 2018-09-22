#'@title Initial value estimate for the psiLearn.
#'@author MingyangLiu <liux3941@umn.edu>
#'@description Initial value estimate for both psi-linear,psi-kernel,and psi-linear-VS.
#'@usage psi_Init(X,A,R,tau=0.1,kappa=seq(0.01,1,length.out=10),sigma=NULL,kernel='linear',lambda=NULL,VS=FALSE)
#'@param X \eqn{n} by \eqn{p} input matrix.
#'@param A a vector of n entries coded 1 and -1 for the treatment assignments.
#'@param R a vector of outcome variable, larger is more desirable.
#'@param tau tuning parameter for the loss function in psi-Learn
#'@param kappa tunning parameter to control the complexity of the decision function, a seq of kappa is set as the default.
#'@param sigma bandwidth parameter when the kernel is 'rbf',if not provided,can be estimated from \code{\link{Sig_est}}
#'@param kernel kernel function for pai-Learn, can be 'linear' or 'rbf' (radial basis kernel), default is 'linear'.
#'              When 'rbf' is specified, one can specify the sigma parameter for radial basis kernel or estimate from \code{\link{Sig_est}}.
#'@param lambda when using the variable selection method for inital estimate, lambda is regarded as the tunning parameter in controlling the complexity in L1 penalty
#'@param VS wheter to use the variable selection method as the inital estimate,default is FALSE
#'@return It returns initial estimated coefficients for psiLearn with the best tuning parameters picked by cross validation.
#'  \item{w}{the coefficent for the decision function.}
#'  \item{bias}{the intercept in both the linear case and the kernel case.}
#'@seealso \code{\link{tlp_S}}
#'@import CVXR
#'@export



psi_Init<-function(X,A,R,tau=0.01,kappa=seq(0.01,1,length.out=10),sigma=NULL,kernel='linear',lambda=NULL,VS=FALSE){
  n=dim(X)[1]
  p=dim(X)[2]
  pai=A*sum(A==1)/length(A)+(1-A)/2
  fitval=MainEff(X,A,R)
  r=R-fitval
  wt=R/pai
  wtr=r/pai
  Value=-Inf
  if (is.null(sigma) & kernel=='rbf'){
    sigma=Sig_est(X,A)
  }
  if(is.null(lambda) & kernel=='linear'&VS==TRUE){
    lambda=0.1
  }


    for (i in 1:length(kappa)){
      kappai=kappa[i]
      c=1/(tau*kappai*n)
      if (VS==FALSE){
      result=DTRlearn::wsvm(X,A,wtr,kernel,sigma=sigma,C=c)
      #result=DTRlearn::Olearning_Single(X,A,R,kernel=kernel,sigma=sigma,m=4,clinear=c,e=1e-6)
      if(kernel=='linear'){
        w_old=result$beta
      } else{
        w_old=result$alpha1
      }

      bias_old=result$bias
      fit=sign(result$fit)

      }else{
       # suppressMessages(suppressWarnings(library(CVXR)))

        if(VS==TRUE & kernel!='linear'){
          stop(gettextf("in this variable selection, only linear kernel is accepted "))
        }
          w_old=Variable(p)
          bias_old=Variable()
          loss1<-function(X,A,w_old,bias_old){
          pos(1-A*(X%*%w_old+bias_old))
        }
        obj<-kappai*sum(w_old^2)+1/tau*mean(wtr*loss1(X,A,w_old,bias_old))+lambda*p_norm(w_old,1)
        prob<-Problem(Minimize(obj))
        result<-solve(prob)
        w_old<-result$getValue(w_old)
        bias_old<-result$getValue(bias_old)
        fit<-sign(X%*%w_old+bias_old)
    }
      value=sum(R*(A==fit))/sum(A==fit)
      if (Value<value){
        Value=value
        w=w_old
        bias=bias_old
     }
    }

  return(list(w=w,bias=bias))
}

#'@title  main effects estimate  in psiLearn
#'@author MingyangLiu <liux3941@umn.edu>
#'@description 	Suppose the main effects is linear, then return the fitted values in the weighted linear regression which
#'              is the main effects. Then the residual will be the interaction effects between the treatment and the "reward" plus the
#'              zero mean random effects.
#'@usage MainEff(X,A,R)
#'@param X n*p input matrix.
#'@param A a vector of n entries coded 1 and -1 for the treatment assignments.
#'@param R a vector of outcome variable, larger is more desirable.
#'@return It returns  the main effects estimated using weighted linear regression.
#'  \item{fitval}{the estimated main effects}
#'@export

MainEff<-function(X,A,R){
  wholedata=data.frame(X,R)
  p=A*sum(A==1,na.rm=TRUE)/length(A)+(1-A)/2
  model=lm(R~.,weights=1/(2*p),data=wholedata)
  fitval=model$fitted.values
  return(fitval)
}


#'@title  inverse bandwidth parameter estimate for the rbf kernel in psiLearn
#'@author MingyangLiu <liux3941@umn.edu>
#'@description 	The inverse kernel width used in Gaussian kernel
#'@usage Sig_est(X,A)
#'@param X n*p input matrix.
#'@param A a vector of n entries coded 1 and -1 for the treatment assignments.
#'@return It returns  estimated bandwidth parameter sigma for rbf kernel when it is not provided by the user
#'  \item{sigma}{the estimated bandwidth parameter.}
#'@export
Sig_est<-function(X,A){
  posID=(A>0)
  negID=1-posID
  Xpos=X[posID,]
  Xneg=X[negID,]
  m  <- nrow(Xpos); n <- nrow(Xneg)
  XY <- Xpos %*% t(Xneg)
  XX <- matrix( rep(apply(Xpos*Xpos, 1, sum), n), m, n, byrow=F )
  YY <- matrix( rep(apply(Xneg*Xneg, 1, sum), m), m, n, byrow=T )

  dis<-sqrt(pmax(XX + YY - 2*XY, 0))
  sigma<-2*(median(dis/2,na.rm=TRUE))^2
  sigma=1/sigma
  return(sigma)
}

#'@title Cost function in psi-Learn
#'@author MingyangLiu <liux3941@umn.edu>
#'@description  calculate the cost in psi-Learn with ridge penalty
#'@usage tlp_S(X,A,R,wt,w,tau=0.1,kappa=0.1,kernel='linear')
#'@param X \eqn{n} by \eqn{p} input matrix.
#'@param A a vector of n entries coded 1 and -1 for the treatment assignments.
#'@param R a vector of outcome variable, larger is more desirable.
#'@param wt a vector of weights for each observation.
#'@param w coefficients for the decision function, the first element is the bias
#'@param tau tuning parameter for the loss function in psi-Learn
#'@param kappa tunning parameter to control the complexity of the decision function.
#'@param kernel kernel function for pai-Learn, can be 'linear' or 'rbf' (radial basis kernel), default is 'linear'.
#'@return It returns the cost value in psi-Learn
#'  \item{cost}{the cost value in psi-Learn}
#'@export


tlp_S<-function(X,A,R,wt,w,tau=0.1,kappa=0.1,kernel='linear'){
  p=dim(X)[2]
  bee=w[1]
  wstar=w[2:length(w)]

  temp=drop(X%*%wstar)
  u=A*(temp+bee);
  phi1=pmax(1-u,0);
  phi2=pmax(-u,0);
  if (kernel=='linear'){
  cost=kappa/2*sum(wstar^2)+mean((phi1-phi2)*wt/tau);
  }else{

    cost=kappa/2*t(wstar)%*%X%*%wstar+mean((phi1-phi2)*wt/tau);
  }
  return(cost)
}

#'@title Cost function in psi-Learn with variable selection
#'@author MingyangLiu <liux3941@umn.edu>
#'@description  calculate the cost in psi-Learn with ridge penalty and TLP penalty
#'@usage tlp_S_VS(X,A,R,wt,w,tau=0.1,kappa=0.1,lambda=0.1,tau2=0.1,kernel='linear')
#'@param X \eqn{n} by \eqn{p} input matrix.
#'@param A a vector of n entries coded 1 and -1 for the treatment assignments.
#'@param R a vector of outcome variable, larger is more desirable.
#'@param wt a vector of weights for each observation.
#'@param w coefficients for the decision function, the first element is the bias
#'@param tau tuning parameter for the loss function in psi-Learn
#'@param kappa tunning parameter to control the complexity of the decision function in the ridge penaly
#'@param lambda tunning parameter to control the complexity of the decision function in the TLP penalty
#'@param tau2 tunning parameter to control in margin in the TLP penalty
#'@param kernel kernel function for pai-Learn, can be 'linear' or 'rbf' (radial basis kernel), default is 'linear'.
#'@return It returns the cost value in psi-Learn
#'  \item{cost}{the cost value in psi-Learn}
#'@export
#'@importFrom Rmpfr pmin


tlp_S_VS<-function(X,A,R,wt,w,tau=0.1,kappa=0.1,lambda=0.1,tau2=0.1,kernel='linear'){
  p=dim(X)[2]
  bee=w[1]
  wstar=w[2:length(w)]

  temp=drop(X%*%wstar)
  u=A*(temp+bee);
  phi1=pmax(1-u,0);
  phi2=pmax(-u,0);
  if (kernel!='linear'){
    stop(gettextf("Only linear kernel is used here "))
  }else{
    cost=kappa/2*sum(wstar^2)+mean((phi1-phi2)*wt/tau)+lambda/tau2*sum(Rmpfr::pmin(abs(wstar),tau2))
  }
  return(cost)
}
#'@title Residual estimate after subtracting the linear main effect
#'@author MingyangLiu <liux3941@umn.edu>
#'@description  calculate the residual  after subtracting the linear main effect
#'@usage resEst(X,R,p)
#'@param X \eqn{n} by \eqn{p} input matrix.
#'@param R a vector of outcome variable, larger is more desirable.
#'@param p a vector of probability of P(A|X) for each observation.
#'@return It returns the residual in psi-Learn
#'  \item{res}{the residual in psi-Learn}
#'@export
#'@importFrom glmnet cv.glmnet
resEst<-function(X,R,p){

  model=cv.glmnet(X,R,weights=p,nfolds=5)

  fit=predict(model,type="response",newx=X)
  res=R-fit

  return(res)

}

