\name{run.peaks}
\alias{run.peaks}
\title{Locate Potential Peaks in FT-ICR MS Spectra}
\description{Takes baseline-corrected data and locates potential peaks in the spectra.}

\usage{
run.peaks(trans.method = c("shiftedlog", "glog", "none"),
          add.par = 0, subtract.base = FALSE, root.dir = ".", 
          base.dir, peak.dir, overwrite = FALSE,
          use.par.file = FALSE, par.file = "parameters.RData",
          num.pts = 5, R2.thresh = 0.98, oneside.min = 1,
          peak.method = c("parabola", "locmaxes"), 
          calc.all.peaks = FALSE, gengamma.quantiles = TRUE, 
          peak.thresh = 3.798194)
}

\arguments{
    \item{trans.method}{type of transformation to use on spectra before statistical analysis}
    \item{add.par}{additive parameter for \code{"shiftedlog"} or \code{"glog"} options for \code{trans.method}}
    \item{subtract.base}{logical; whether to subtract calculated baseline from spectrum}
    \item{root.dir}{directory for parameters file and raw data}
    \item{base.dir}{directory for baseline files; default is \code{paste(root.dir, "/Baselines", sep = "")}}
    \item{peak.dir}{directory for peak location files; default is \code{paste(root.dir, "/All_Peaks", sep = "")}}
    \item{overwrite}{logical; whether to replace existing files with new ones}
    \item{use.par.file}{logical; if \code{TRUE}, then parameters are read from \code{par.file} in directory \code{root.dir}}
    \item{par.file}{string containing name of parameters file}
    \item{num.pts}{number of consecutive points needed for peak fitting}
    \item{R2.thresh}{\eqn{R^2} value needed for peak fitting}
    \item{oneside.min}{minimum number of points on each side of local maximum for peak fitting}
    \item{peak.method}{method for locating peaks}
    \item{calc.all.peaks}{logical; whether to calculate all possible peaks or only sufficiently large ones}
    \item{gengamma.quantiles}{logical; whether to use generalized gamma quantiles when calculating large peaks}
    \item{peak.thresh}{threshold for declaring large peak; see below}
}

\details{
Reads in information from each file created by \code{\link{run.baselines}}, 
calls \code{\link{locate.peaks}} to find potential peaks, and writes the output 
to a file in directory \code{peak.dir}.  The name of each new file is the same 
as the name of the old file with \dQuote{.RData} replaced by 
\dQuote{\_peaks.RData}.  The resulting file contains the data frame 
\code{all.peaks}, which has columns
\tabular{ll}{ \tab \cr
    \code{Center_hat} \tab estimated mass of peak \cr
    \code{Max_hat} \tab estimated intensity of peak \cr
    \code{Width_hat} \tab estimated width of peak \cr
}    
and is ready to be used by \code{\link{run.lrg.peaks}}.
  
The parameters \code{gengamma.quantiles} and \code{peak.thresh} are relevant 
only if \code{calc.all.peaks = FALSE}.  In that case, if 
\code{gengamma.quantiles == TRUE}, then \code{peak.thresh} is interpreted as a 
multiplier for the baseline.  Anything larger than \code{peak.thresh} times the 
estimated baseline is declared to be a real peak.  
If \code{gengamma.quantiles == FALSE}, then \code{peak.thresh} is interpreted as 
two-thirds of the value of \eqn{K} used in a Tukey's biweight estimation of 
center and scale (so roughly equal to the number of standard deviations above 
the mean for iid normal data).  Anything with weight zero in the calculation is 
then declared to be a real peak.
}

\value{
No value returned; the files are simply created.
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
If \code{use.par.file == TRUE} and other parameters are entered into the function
call, then the parameters entered in the function call overwrite those read in
from the file.  Note that this is opposite from the behavior for
\pkg{\link{FTICRMS}} versions 0.7 and earlier.

\code{peak.method} and \code{trans.method} can be abbreviated.

Using \code{calc.all.peaks == FALSE} can speed up computation time immensely,
but will affect the final result.  It probably won't affect it much, but
\emph{caveat emptor}.
}

\seealso{
\code{\link{run.baselines}}, \code{\link{run.lrg.peaks}}, \code{\link{locate.peaks}}
}

\examples{}
