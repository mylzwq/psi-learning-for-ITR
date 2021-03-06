% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitialEst.R
\name{tlp_S_VS}
\alias{tlp_S_VS}
\title{Cost function in psi-Learn with variable selection}
\usage{
tlp_S_VS(X,A,R,wt,w,tau=0.1,kappa=0.1,lambda=0.1,tau2=0.1,kernel='linear')
}
\arguments{
\item{X}{\eqn{n} by \eqn{p} input matrix.}

\item{A}{a vector of n entries coded 1 and -1 for the treatment assignments.}

\item{R}{a vector of outcome variable, larger is more desirable.}

\item{wt}{a vector of weights for each observation.}

\item{w}{coefficients for the decision function, the first element is the bias}

\item{tau}{tuning parameter for the loss function in psi-Learn}

\item{kappa}{tunning parameter to control the complexity of the decision function in the ridge penaly}

\item{lambda}{tunning parameter to control the complexity of the decision function in the TLP penalty}

\item{tau2}{tunning parameter to control in margin in the TLP penalty}

\item{kernel}{kernel function for pai-Learn, can be 'linear' or 'rbf' (radial basis kernel), default is 'linear'.}
}
\value{
It returns the cost value in psi-Learn
 \item{cost}{the cost value in psi-Learn}
}
\description{
calculate the cost in psi-Learn with ridge penalty and TLP penalty
}
\author{
MingyangLiu <liux3941@umn.edu>
}
