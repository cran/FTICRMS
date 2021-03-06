\name{locate.peaks}
\alias{locate.peaks}
\title{Locate Peaks in a FT-ICR MS Spectrum}
\description{Locates peaks in FT-ICR MS spectra assuming that the peaks are roughly parabolic on the log scale.}

\usage{
locate.peaks(peak.base, num.pts = 5, R2.thresh = 0.98, 
             oneside.min = 1, peak.method = c("parabola", "locmaxes"), 
             thresh = -Inf)
}

\arguments{
  \item{peak.base}{numeric matrix with two columns containing the masses and the transformed spectrum intensities}
  \item{num.pts}{minimum number of points needed to have a peak}
  \item{R2.thresh}{minimum \eqn{R^2} needed to have a peak}
  \item{oneside.min}{minimum number of points needed on each side of the local maximum}
  \item{peak.method}{how to locate peaks}
  \item{thresh}{only local maxes that are larger than this will be checked to see if they are peaks}
}

\details{If \code{peak.method == "parabola"}, the algorithm works by locating local maxima in the spectrum,
then seeing if any \code{num.pts} consecutive points with at least \code{oneside.min} point(s) on each side 
of the local maximum have a coefficient of determination (\eqn{R^2}) of at least \code{R2.thresh} when fitted with a 
quadratic.  If, in addition, the coefficient of the squared term is negative, then this is declared a peak 
and the vertex of the corresponding parabola is located.  The coordinates of the vertex give the components 
\code{Center_hat} and \code{Max_hat} in the return value.  The \code{Width_hat} component is the negative 
reciprocal of the coefficient of the squared term.

If \code{peak.method == "locmax"}, then the algorithm merely returns the set of local maxima larger than
\code{thresh}, and the \code{Width_hat} component of the return value is \code{NA}.
}

\value{
A data frame with columns
    \item{Center_hat}{estimated mass of peak}
    \item{Max_hat}{estimated intensity of peak}
    \item{Width_hat}{estimated width of peak}
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
}

\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}

\note{
An extremely large value for \code{Width_hat} most likely indicates a bad fit.

\code{peak.method} can be abbreviated.  Using \code{peak.method = "locmax"} will 
vastly speed up the runtime, but may affect the quality of the analysis.

As noted in both papers in the References, a typical FT-ICR MS spectrum has far
more peaks than can be accounted for by actual compounds.  Thus, defining a good
value of \code{thresh} will vastly speed up the computation without materially
affecting the analysis.
}

\seealso{
\code{\link{run.peaks}}
}

\examples{}
