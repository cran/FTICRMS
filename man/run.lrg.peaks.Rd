\name{run.lrg.peaks}
\alias{run.lrg.peaks}
\title{Extract "Large" Peaks from Files}
\description{
Takes the files output by \code{\link{run.peaks}}, extracts \dQuote{large} peaks, combines them
into a single data frame, and writes the data frame to a file.
}
\usage{
run.lrg.peaks(trans.method = c("shiftedlog", "glog", "none"), 
              add.par = 0, subtract.base = FALSE,
              root.dir = ".", peak.dir, base.dir, lrg.dir,
              lrg.file = lrg_peaks.RData, overwrite = FALSE,
              use.par.file = FALSE, par.file = "parameters.RData",
              calc.all.peaks = FALSE, gengamma.quantiles = TRUE,
              peak.thresh = 3.798194, subs)
}
\arguments{
    \item{trans.method}{type of transformation to use on spectra before statistical analysis}
    \item{add.par}{additive parameter for \code{"shiftedlog"} or \code{"glog"} options for \code{trans.method}}
    \item{subtract.base}{logical; whether to subtract calculated baseline from spectrum}
    \item{root.dir}{directory for parameters file and raw data}
    \item{peak.dir}{directory for peak location files; default is \code{paste(root.dir, "/All_Peaks", sep = "")}}
    \item{base.dir}{directory for baseline files; default is \code{paste(root.dir, "/Baselines", sep = "")}}
    \item{lrg.dir}{directory for large peaks file; default is \code{paste(root.dir, "/Large_Peaks", sep = "")}}
    \item{lrg.file}{name of file to store large peaks in}
    \item{overwrite}{logical; whether to replace existing files with new ones}
    \item{use.par.file}{logical; if \code{TRUE}, then parameters are read from \code{par.file} in directory \code{root.dir}}
    \item{par.file}{string containing name of parameters file}
    \item{calc.all.peaks}{logical; whether to calculate all possible peaks or only sufficiently large ones}
    \item{gengamma.quantiles}{logical; whether to use generalized gamma quantiles when calculating large peaks}
    \item{peak.thresh}{threshold for declaring large peak; see below}
    \item{subs}{subset of spectra to use for analysis; see below}
}
\details{
Reads in information from each file created by \code{\link{run.peaks}},
extracts peaks which are \dQuote{large} (see below), and creates the file
\code{lrg.file} in \code{lrg.dir}.  The resulting file contains the data frame
\code{lrg.peaks}, which has columns
\tabular{ll}{ \tab \cr
    \code{Center_hat} \tab estimated mass of peak \cr
    \code{Max_hat} \tab estimated intensity of peak \cr
    \code{Width_hat} \tab estimated width of peak \cr
    \code{File} \tab name of file the peak was extracted from, with \dQuote{\_peaks.RData} deleted \cr
}    
and is ready to be used by \code{\link{run.strong.peaks}}.
}

\value{
No value returned; the file is simply created.
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
If \code{use.par.file == TRUE} and other parameters are entered into the 
function call, then the parameters entered in the function call overwrite those 
read in from the file.  This is opposite from the behavior for 
\pkg{\link{FTICRMS}} versions 0.7 and earlier.

\code{trans.method} can be abbreviated.

If \code{gengamma.quantiles == TRUE}, then a peak is \dQuote{large} if it is at
least \code{peak.thresh} times as large as the estimated baseline at that point.

If \code{gengamma.quantiles == FALSE}, then a peak is \dQuote{large} if it has
zero weight in the data generated by \code{\link{run.peaks}} for the spectrum it
comes from when using Tukey's biweight with parameter 
\code{K = 1.5 * peak.thresh} to estimate center and scale.

If \code{subs} is not defined, the algorithm finds large peaks for all files in
\code{peak.dir}.  If it is defined, \code{subs} can be logical or numeric or
character; if it is defined, then the algorithm finds large peaks for all
entries in \code{subs} (character) or \code{list.files(peak.dir)[subs]} (logical
or numeric).
}

\seealso{
\code{\link{run.peaks}}, \code{\link{run.cluster.matrix}}
}

\examples{}
