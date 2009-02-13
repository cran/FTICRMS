\name{run.all}
\alias{run.all}
\title{Complete Analysis of FT-ICR MS Data}
\description{
A wrapper that calls all six functions needed for a full analysis.
}
\usage{
run.all(par.file = "parameters.RData", root.dir = ".")
}
\arguments{
  \item{par.file}{string containing the name of the parameters file}
  \item{root.dir}{string containing location of raw data directory and parameters file}
}
\details{Requires \code{par.file} to be in place before starting---for example by creating it with
\code{\link{make.par.file}}.

Calls (in order) \code{\link{run.baselines}}, \code{\link{run.peaks}}, \code{\link{run.lrg.peaks}}, 
\code{\link{run.strong.peaks}}, \code{\link{run.cluster.matrix}}, and \code{\link{run.analysis}}.}

\references{
Barkauskas, D.A. (2009) \dQuote{Statistical Analysis of Matrix-Assisted Laser Desorption/Ionization
Fourier Transform Ion Cyclotron Resonance Mass Spectrometry Data with Applications to Cancer
Biomarker Detection}.  Ph.D. dissertation, University of California at Davis.

Barkauskas, D.A. \emph{et al}. (2009) \dQuote{Detecting glycan cancer biomarkers in serum
samples using MALDI FT-ICR mass spectrometry data}.  \emph{Bioinformatics}, \bold{25}:2, 251--257.

Benjamini, Y. and Hochberg, Y.  (1995) \dQuote{Controlling the false discovery rate: a practical and powerful 
approach to multiple testing.}  \emph{J. Roy. Statist. Soc. Ser. B}, \bold{57}:1, 289--300.

Xi, Y. and Rocke, D.M. (2008) \dQuote{Baseline Correction for NMR Spectroscopic Metabolomics Data Analysis}.  \emph{BMC Bioinformatics}, \bold{9}:324.
}
\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}
\note{
The analysis described in Barkauskas \emph{et al}. (2008) can be reproduced using the following parameter
values instead of the defaults:
\tabular{ll}{
    \tab \code{add.par = 10} \cr
    \tab \code{calc.all.peaks = TRUE} \cr
    \tab \code{gengamma.quantiles = FALSE} \cr
    \tab \code{neg.norm.by = "constant"} \cr
    \tab \code{peak.thresh = 4} \cr
    \tab \code{sm.norm.by = "constant"} \cr
    \tab \code{subtract.base = TRUE}
}
}
\seealso{\code{\link{make.par.file}}, \code{\link{run.baselines}}, \code{\link{run.peaks}}, 
\code{\link{run.lrg.peaks}}, \code{\link{run.strong.peaks}}, \code{\link{run.cluster.matrix}}, 
\code{\link{run.analysis}}}
\examples{}
