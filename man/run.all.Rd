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
Barkauskas, D.A. \emph{et al}. (2008) \dQuote{Detecting glycan cancer biomarkers in serum 
samples using MALDI FT-ICR mass spectrometry data}.  Submitted to \emph{Bioinformatics}

Xi, Y. and Rocke, D.M. (2008) \dQuote{Baseline Correction for NMR Spectroscopic Metabolomics Data Analysis}.  \emph{BMC Bioinformatics}, \bold{9}:324.

Benjamini, Y. and Hochberg, Y.  (1995) \dQuote{Controlling the false discovery rate: a practical and powerful 
approach to multiple testing.}  \emph{J. Roy. Statist. Soc. Ser. B}, \bold{57}:1, 289--300.
}
\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}
\note{}
\seealso{\code{\link{make.par.file}}, \code{\link{run.baselines}}, \code{\link{run.peaks}}, 
\code{\link{run.lrg.peaks}}, \code{\link{run.strong.peaks}}, \code{\link{run.cluster.matrix}}, 
\code{\link{run.analysis}}}
\examples{}
