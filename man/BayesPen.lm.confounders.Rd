\name{BayesPen.lm.confounders}
\alias{BayesPen.lm.confounders}
\title{
Confounder Selection via Penalized Credible Regions for Linear Models
}
\description{
Fits a linear outcome and exposure model using code based on BLR (Campos and Rodriguez 2012) and performs confounder selection via penalized credible regions (Wilson and Reich 2014+).
}
\usage{
BayesPen.lm.confounders(y, x, u, prior, nIter, burnIn, thin, update, force,
  	    max.steps = NULL, max.refit, saveAt = "", include.me = FALSE, 
  	    z.score = FALSE)
}
\arguments{
  \item{y}{
A n-vector of responses.
}
  \item{x}{
For single exposures this is a n-vector of exposures.  For multiple exposures this is a design matrix of the exposures with n rows and one column for each exposure.
}
  \item{u}{
A design matrix with potential confounders and other covariates. This should include a column of ones if an intercept is desired.
}
  \item{prior}{
A list containing the priors for the regression coefficients and the error variance. Each object in the list is also a list. 
The list prior$varE contains the degrees of freedom (prior$varE$df) and scale (prior$varE$S) for the inverse-X^2 prior on the error variance. 
The list prior$varBR contains the degrees of freedom (prior$varBR$df) and scale (prior$varBR$S) for the inverse-X^2 prior on the variance of the regression coefficients. 
In both cases the parameterization used has prior mean variance S/(df-2). 
See the documentation for BLR for additional details. 
Only one prior can be specified and will be used for both the outcome and exposure models.
}
  \item{nIter}{
The number of MCMC iterations (integer).
}
  \item{burnIn}{
The number of MCMC iterations to be discarded as burnin (integer).}
  \item{thin}{
Thinning number for the MCMC chain (integer)
}
  \item{update}{ 
Integer specifying how often to print an update on the progress of the MCMC.
}
  \item{force}{
  An optional vector indexing which confounding variables should be forced into the model, for example an intercept may be forced into the model. 
  This vector indexes the columns of u. All exposures are forced into the model automatically.
}
  \item{max.steps}{
Maximum number of steps to be performed in the LARS algorithm (Hastie and Efron 2013).
}
  \item{max.refit}{ The maximum number of models to be refit.}
  \item{saveAt}{
  An optional string that is a pre-fix for the filenames for the files saved while the program runs. 
}
  \item{include.me}{
  Indicator for the multiple exposure case indicating if the exposure models should include the other exposures as covariates.
}
  \item{z.score}{
Indicator for using z-scores in the penalty instead of regression coefficients. 
If the error variance for the exposure and outcome models or the scale of the regression coefficients is very different then using the z-scores make the confounder and exposure components of the weights more comparable.} 
}
\value{
\item{joint.path}{A complete solution path for the joint credible regions approach. Each row is a model in the solution path with a 1 indicating a variable is included and a 0 indicating it is not included.}
\item{marginal.path}{Not used for confounder selection.}
\item{order.path}{The action returned from lars that shows when each covariate is added to the model.}
\item{order.marg}{Not used for confounder selection.}
 \item{joint }{Always TRUE for confounder selection.}
 \item{force}{Vector of variables forced into the model including the exposures.}
 \item{coefs }{A matrix of regression coefficients for each model in the solution path. The regression coefficients for parameters omitted from a model are set to 0.}
\item{SSE }{SSE of each refitted model.}
\item{dev }{Deviance of each refitted model.}
\item{df }{ Error degrees of freedom from each refitted model.}
 \item{lm}{ Full fitting of the outcome model.}
 \item{confounder.weights}{Confounder weights used in BayesPen to fit the model.}
}
\references{
Gustavo de los Campos and Paulino Perez Rodriguez, (2012). BLR: Bayesian Linear Regression. R package version 1.3. http://CRAN.R-project.org/package=BLR

Wilson, A. and Reich, B.J. (2014+). Confounder selection via penalized credible regions.
}
\author{
Ander Wilson, Howard D. Bondell, and Brian J. Reich
}

\examples{
set.seed(1234)
dat <- SimExample(500,model="WPD2")
fit.BRnew <- BayesPen.lm.confounders(y=dat$y,x=dat$X, u=dat$U)
BayesPen.plot(fit.BRnew)
}

\keyword{ regression }