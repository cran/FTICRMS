\name{baseline}
\alias{baseline}
\title{Calculate Baselines for FT-ICR MS Spectra}
\description{
Computes an estimated baseline curve for a spectrum by a method of Rocke and Xi.
}
\usage{
baseline(spect, sm.fac = 10^15, neg.pen = sqrt(pi/2), 
         max.iter = 30, frac.changed = 0.001, binsize = 64, 
         xvals = 1:length(spect))
}
\arguments{
  \item{spect}{vector containing the intensities of the spectrum}
  \item{sm.fac}{smoothing factor}
  \item{neg.pen}{negativity penalty}
  \item{max.iter}{maximum number of iterations before giving up}
  \item{frac.changed}{maximal fraction of indicators that can switch for convergence}
  \item{binsize}{bin size used to estimate noise variance}
  \item{xvals}{vector containing the mass locations of the spectrum; currently not used in the algorithm}
}
\details{
The algorithm works by maximizing the objective function 
\deqn{F(\{b_{i}\}) = \sum_{i=1}^{n}b_{i} - A_{1}\sum_{i=2}^{n-1}(b_{i-1}-2b_{i}+b_{i+1})^{2} - A_{2}
\sum_{i=1}^{n}[\max\{b_{i}-y_{i},0\}]^{2}}{F({b[i]}) = sum_{i=1}^{n}b[i] - A[1] sum_{i=2}^{n-1}(b[i-1]-2b[i]+b[i+1])^2 - A[2]
\sum_{i=1}^n[max{b[i]-y[i],0}]^2}

The value of \code{xvals} is currently irrelevant; future implementations may use it.

Experimentally, the smoothing factor \code{sm.fac} should be as large as possible without giving
obviously incorrect results.

The default value of \code{neg.pen} is calculated under the assumption of iid \eqn{N(0,1)} noise.

The default value of \code{frac.changed} is the smallest one for which the algorithm reliably terminates
before \code{max.iter} iterations are reached for any reasonable value of \code{max.iter}. 
}

\value{
  A list containing the following items:
  \item{baseline}{The computed baseline}
  \item{noise}{The estimated noise standard deviation}
  \item{iter}{The number of iterations for convergence}
  \item{changed}{Numeric vector of length \code{iter} containing the number of indicator variables that switched value on each iteration}
}
\references{
Barkauskas, D.A. \emph{et al}. (2008) \dQuote{Detecting glycan cancer biomarkers in serum 
samples using MALDI FT-ICR mass spectrometry data}.  Submitted to \emph{Bioinformatics}

Xi, Y. and Rocke, D.M. (2006) \dQuote{Baseline Correction for NMR Spectroscopic Metabolomics Data Analysis}, unpublished.
}
\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}
\note{The algorithm used was developed by Yuanxin Xi and David Rocke (unpublished).  This program is based
on a Matlab program by Yuanxin Xi.}
\seealso{\code{\link{run.baselines}}}
\examples{}
