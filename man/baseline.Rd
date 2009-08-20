\name{baseline}
\alias{baseline}
\title{Calculate Baselines for Spectroscopic Data}

\description{
Computes an estimated baseline curve for a spectrum using the \dQuote{BXR
algorithm,} a method of Xi and Rocke generalized by Barkauskas and Rocke.
}

\usage{
baseline(spect, init.bd, sm.par = 1e-11, sm.ord = 2, max.iter = 20, tol = 5e-8,
         sm.div = NA, sm.norm.by = c("baseline", "overestimate", "constant"),
         neg.div = NA, neg.norm.by = c("baseline", "overestimate", "constant"),
         rel.conv.crit = TRUE, zero.rm = TRUE, halve.search = FALSE)
}

\arguments{
    \item{spect}{vector containing the intensities of the spectrum}
    \item{init.bd}{initial value for baseline; default is flat baseline at median height}
    \item{sm.par}{smoothing parameter for baseline calculation}
    \item{sm.ord}{order of derivative to penalize in baseline analysis}
    \item{max.iter}{convergence criterion in baseline calculation}
    \item{tol}{convergence criterion; see below}
    \item{sm.div}{smoothness divisor in baseline calculation}
    \item{sm.norm.by}{method for smoothness penalty in baseline analysis}
    \item{neg.div}{negativity divisor in baseline calculation}
    \item{neg.norm.by}{method for negativity penalty in baseline analysis}
    \item{rel.conv.crit}{logical; whether convergence criterion should be relative to size of current baseline estimate}
    \item{zero.rm}{logical; whether to replace zeros with average of surrounding values}
    \item{halve.search}{logical; whether to use a halving-line search if step leads to smaller value of function}
}

\details{
If the spectrum is given by \eqn{y_{i}}{y[i]}, then the algorithm works by
maximizing the objective function \deqn{F(\{b_{i}\}) = \sum_{i=1}^{n}b_{i} -
\sum_{i=2}^{n-1}A_{1,i}(b_{i-1}-2b_{i}+b_{i+1})^{2} - \sum_{i=1}^{n}A_{2,i}
[\max\{b_{i}-y_{i},0\}]^{2}}{F({b[i]}) = sum_{i=1}^{n}b[i] -
sum_{i=2}^{n-1}A[1,i]*(b[i-1]-2b[i]+b[i+1])^2 -
\sum_{i=1}^n A[2,i]*[max{b[i]-y[i],0}]^2} using Newton's method (with embedded
halving line search if \code{halve.search == TRUE}) using starting value
\code{b[i] = init.bd[i]} for all \eqn{i}.  The middle term controls the
smoothness of the baseline and the last term applies a \dQuote{negativity
penalty} when the baseline is above the spectrum.

The smoothing factor \code{sm.par} corresponds to \eqn{A_{1}^{*}}{A[1]^{*}} in
Barkauskas (2009) and controls how large the estimated \emph{n}th derivative of
the baseline is allowed to be (for \code{sm.ord = n}).  From a practical
standpoint, values of \code{sm.ord} larger than two do not seem to adequately
smooth the baseline because the Hessian becomes computationally singular for any
reasonable value of \code{sm.par}.

The parameters \code{sm.div}, \code{sm.norm.by}, \code{neg.div}, and
\code{neg.norm.by} determine the methods used to normalize the smoothness and
negativity terms.  The general forms are
\eqn{A_{1,i} = n^{4}A_{1}^{*}/M_{i}/p}{A[1,i] = n^4 * A[1]^{*}/M[i]/p} and
\eqn{A_{2,i} = 1/M_{i}/p}{A[2,i] = 1/M[i]/p}.  Here, \code{n = length(spect)};
\eqn{p} is \code{sm.div} or \code{neg.div}, as appropriate; and
\eqn{M_{i}}{M[i]} is determined by \code{sm.norm.by} or \code{neg.norm.by}, as
appropriate.  Values of \code{"baseline"} make
\eqn{M_{i} = b_{i}'}{M[i] = b[i]'}, where \eqn{b_{i}'}{b[i]'} is the currently
estimated value of the baseline; values of \code{"overestimate"} make
\eqn{M_{i} = b_{i}'-y_{i}}{M[i] = b[i]'-y[i]}; and values of \code{"constant"}
make \eqn{M_{i} = \sigma}{M[i] = \sigma}, where \eqn{\sigma} is an estimate of
the noise standard deviation.

The values of \code{sm.norm.by} and \code{neg.norm.by} can be abbreviated and
both have default value \code{"baseline"}.  The default values of \code{NA} for
\code{sm.div} and \code{neg.div} are translated by default to
\code{sm.div = 0.5223145} and \code{neg.div = 0.4210109}, which are the
appropriate parameters for the FT-ICR mass spectrometry machine that generated
the spectra which were used to develop this package.  It is distinctly possible
that other machines will require different parameters, and almost certain that
other spectroscopic technologies will require different parameters; see
Barkauskas (2009a) for a description for how these parameters were obtained.

If \code{zero.rm == TRUE} and \eqn{y_{a},\ldots,y_{a+k}=0}{y[a],\dots,y[a+k] = 0},
then these values of the spectrum are set to be
\eqn{(y_{a-1}+y_{a+k+1})/2}{(y[a-1]+y[a+k+1])/2}.  (For typical MALDI FT-ICR
spectra, a spectrum value of zero indicates an erased harmonic and should not be
considered a real data point.)
}

\value{
A list containing the following items:
    \item{baseline}{The computed baseline}
    \item{iter}{The number of iterations for convergence}
    \item{changed}{Numeric vector of length \code{iter} containing the number of indicator variables that switched value on each iteration}
    \item{hs}{Numeric vector of length \code{iter} containing the number of halving line-searches done on each iteration}
}

\references{
Barkauskas, D.A. and D.M. Rocke.  (2009a) \dQuote{A general-purpose baseline 
estimation algorithm for spectroscopic data}.  to appear in \emph{Analytica 
Chimica Acta}.  doi:10.1016/j.aca.2009.10.043   

Barkauskas, D.A. \emph{et al}. (2009b) \dQuote{Analysis of MALDI FT-ICR mass 
spectrometry data: A time series approach}.  \emph{Analytica Chimica Acta}, 
\bold{648}:2, 207--214.

Barkauskas, D.A. \emph{et al}. (2009c) \dQuote{Detecting glycan cancer 
biomarkers in serum samples using MALDI FT-ICR mass spectrometry data}.  
\emph{Bioinformatics}, \bold{25}:2, 251--257.

Xi, Y. and Rocke, D.M. (2008) \dQuote{Baseline Correction for NMR Spectroscopic 
Metabolomics Data Analysis}.  \emph{BMC Bioinformatics}, \bold{9}:324.
}

\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}

\note{
The original algorithm was developed by Yuanxin Xi and David Rocke.  The code in
this package was first adapted from a Matlab program by Yuanxin Xi, then
modified to account for the new methodology in Barkauskas (2009a).

\code{halve.search = FALSE} is recommended unless both
\code{sm.norm.by == "constant"} and \code{neg.norm.by == "constant"}.
}

\seealso{
\code{\link{run.baselines}}
}

\examples{}
