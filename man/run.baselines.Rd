\name{run.baselines}
\alias{run.baselines}
\title{Create Baseline-Corrected Spectra for FT-ICR MS}
\description{
Takes the spectra from files in \code{raw.dir}, calculates and subtracts the baselines from them, 
and writes the results in the directory \code{base.dir}.}
\usage{
run.baselines(root.dir = ".", raw.dir, base.dir, overwrite = FALSE, 
              use.par.file = FALSE, par.file = "parameters.RData", 
              sm.par = 1.1E-9, neg.pen = sqrt(pi/2), k.biweight = 6,
              max.iter = 30, frac.changed = 0.001)
}
\arguments{
  \item{root.dir}{string containing location of raw data directory}
  \item{raw.dir}{directory for raw data files; default is \code{paste(root.dir, "/Raw_Data", sep = "")}}
  \item{base.dir}{directory for baseline-corrected files; default is \code{paste(root.dir, "/Baseline_Corrected", sep = "")}}
  \item{overwrite}{logical; whether to replace existing files with new ones}
  \item{use.par.file}{logical; if \code{TRUE}, then parameters are read from 
      \code{par.file} in directory \code{root.dir}}
  \item{par.file}{string containing name of parameters file}
  \item{sm.par}{smoothing parameter for baseline calculation}
  \item{neg.pen}{negativity penalty in baseline calculation}
  \item{k.biweight}{parameter for robust center and scale estimation using Tukey's biweight}
  \item{max.iter}{convergence criterion in baseline calculation}
  \item{frac.changed}{convergence criterion in baseline calculation}
}
\details{Goes through the entire directory \code{raw.dir} file-by-file and computes each baseline using \code{\link{baseline}}, 
subtracts the baseline from the spectrum, and writes the resulting data frame to a file in directory \code{base.dir}.  
The name of the new file is the same as the name of the old file with \dQuote{.txt} replaced by \dQuote{.RData}, and the 
new file is ready to be used by \code{\link{run.peaks}}.}
\value{
No value returned; the files are simply created.
}
\references{
Barkauskas, D.A. \emph{et al}. (2008) \dQuote{Detecting glycan cancer biomarkers in serum 
samples using MALDI FT-ICR mass spectrometry data}.  Submitted to \emph{Bioinformatics}

Xi, Y. and Rocke, D.M. (2008) \dQuote{Baseline Correction for NMR Spectroscopic Metabolomics Data Analysis}.  \emph{BMC Bioinformatics}, \bold{9}:324.
}
\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}
\note{If \code{use.par.file = TRUE}, then the parameters read in from the file overwrite any arguments entered in the
function call.}
\seealso{\code{\link{baseline}}, \code{\link{run.peaks}}}
\examples{}
