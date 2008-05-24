\name{run.peaks}
\alias{run.peaks}
\title{Locate Potential Peaks in FT-ICR MS Spectra}
\description{Takes baseline-corrected data and locates potential peaks in the spectra.}
\usage{
run.peaks(add.par = 10, trans.method = "shiftedlog", root.dir = ".", 
          base.dir, peak.dir, overwrite = FALSE, use.par.file = FALSE, 
          par.file = "parameters.RData", num.pts = 5, R2.thresh = 0.98, 
          oneside.min = 1, peak.method = "parabola", 
          calc.all.peaks = TRUE, numsds = 4)
}
\arguments{
  \item{add.par}{additive parameter for \code{"shiftedlog"} or \code{"glog"} options for \code{trans.method}}
  \item{trans.method}{type of transformation to use on spectra before statistical analysis; 
      currently, only \code{"shiftedlog"}, \code{"glog"}, and \code{"none"} are supported}
  \item{root.dir}{string containing location of raw data directory}
  \item{base.dir}{directory for baseline-corrected files; default is \code{paste(root.dir, "/Baseline_Corrected", sep = "")}}
  \item{peak.dir}{directory for peak location files; default is \code{paste(root.dir, "/All_Peaks", sep = "")}}
  \item{overwrite}{logical; whether to replace existing files with new ones}
  \item{use.par.file}{logical; if \code{TRUE}, then parameters are read from 
      \code{par.file} in directory \code{root.dir}}
  \item{par.file}{string containing name of parameters file}
  \item{num.pts}{minimum number of points needed to have a peak}
  \item{R2.thresh}{minimum \eqn{R^2} needed to have a peak}
  \item{oneside.min}{minimum number of points needed on each side of the local maximum}
  \item{peak.method}{method to use for locating peaks; currently only \code{"parabola"} and \code{"locmax"} are supported}
  \item{calc.all.peaks}{logical; whether or not to locate all peaks or only those which are sufficiently large}
  \item{numsds}{only relevant if \code{calc.all.peaks = FALSE}; the number of standard deviations above the mean required 
     for a peak to be considered}
}
\details{Reads in information from each file created by \code{\link{run.baselines}}, calls \code{\link{locate.peaks}} to
find potential peaks, and writes the output to a file in directory \code{peak.dir}.  The name of each new file
is the same as the name of the old file with \dQuote{.RData} replaced by \dQuote{\_peaks.RData}.  The resulting file contains
the data frame \code{all.peaks}, which has columns
\tabular{ll}{ \tab \cr
  \code{Center_hat}\tab estimated mass of peak\cr
  \code{Max_hat}\tab estimated intensity of peak\cr
  \code{Width_hat}\tab estimated width of peak\cr
}
  and is ready to be used by \code{\link{run.lrg.peaks}}.
}
\value{
No value returned; the files are simply created.
}
\references{
Barkauskas, D.A. \emph{et al}. (2008) \dQuote{Detecting glycan cancer biomarkers in serum 
samples using MALDI FT-ICR mass spectrometry data}.  Submitted to \emph{Bioinformatics}
}
\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}
\note{If \code{use.par.file = TRUE}, then the parameters read in from the file overwrite any arguments entered in the
function call.

Using \code{calc.all.peaks = FALSE} will speed up computation time immensely, but will affect the final result.  
It probably won't affect it much, but \emph{caveat emptor}.
}
\seealso{\code{\link{run.baselines}}, \code{\link{run.lrg.peaks}}, \code{\link{locate.peaks}} }
\examples{}