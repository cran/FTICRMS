\name{extract.pars}
\alias{extract.pars}
\title{Extract Parameters from File}
\description{
Extracts the parameters in the file specified by \code{par.file} and returns them in list form.
}
\usage{
extract.pars(root.dir = ".", par.file = "parameters.RData")
}
\arguments{
  \item{root.dir}{string containing location of parameters file to be extracted from}
  \item{par.file}{string containing name of parameters file}
}
\details{Used by \code{\link{run.analysis}} to record all the parameter choices in an analysis for future reference.}
\value{
   A list with the following components:
        \item{add.norm}{logical; whether to normalize additively or multiplicatively on the log scale}
        \item{add.par}{additive parameter for \code{"shiftedlog"} or \code{"glog"} options for \code{trans.method}}
        \item{align.method}{alignment algorithm for peaks}
        \item{base.dir}{directory for baseline-corrected files}
        \item{binsize}{size of bins used in estimating noise variance in baseline calculation}
        \item{cor.thresh}{threshhold correlation for declaring isotopes}
        \item{covariates}{data frame containing covariates used in analysis}
        \item{FDR}{False Discovery Rate in Benjamini-Hochberg test}
        \item{form}{formula used in \code{\link{t.test}} or \code{\link{lm}}}
        \item{frac.changed}{convergence criterion in baseline calculation}
        \item{isotope.dist}{maximum distance for declaring isotopes}
        \item{masses}{specific masses to test}
        \item{max.iter}{convergence criterion in baseline calculation}
        \item{neg.pen}{negativity penalty in baseline calculation}
        \item{normalization}{type of normalization to use on spectra before statistical analysis}
        \item{num.pts}{number of points needed for peak fitting}
        \item{numsds}{number of standard deviations above the mean to be declared non-noise}
        \item{oneside.min}{minimum number of points on each side of local maximum for peak fitting}
        \item{overwrite}{whether to replace existing files with new ones}
        \item{peak.dir}{directory for peak location files}
        \item{peak.method}{method for locating peaks}
        \code{pre.align}{shifts to apply before running \code{\link{run.strong.peaks}}}
        \item{R2.thresh}{\eqn{R^2} value needed for peak fitting}
        \item{raw.dir}{directory for raw data files}
        \item{repl.method}{how to deal with replicates}
        \item{res.dir}{directory for result file}
        \item{root.dir}{directory for parameters file and raw data directory}
        \item{lrg.dir}{directory for significant peaks file} 
        \item{lrg.only}{whether to consider only peaks that have at least one peak \dQuote{large}; i.e.,
            identified by \code{\link{run.lrg.peaks}}}
        \item{sm.fac}{smoothing factor for baseline calculation}
        \item{trans.method}{data transformation method}
        \item{use.t.test}{whether to use a \emph{t}-test or an ANOVA/regression model}
}
\references{
Barkauskas, D.A. \emph{et al}. (2008) \dQuote{Detecting glycan cancer biomarkers in serum 
samples using MALDI FT-ICR mass spectrometry data}.  Submitted to \emph{Bioinformatics}

Xi, Y. and Rocke, D.M. (2006) \dQuote{Baseline Correction for NMR Spectroscopic Metabolomics Data Analysis}, unpublished.
}
\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}
\note{\code{do.call(make.par.file, extract.pars())} recreates the original parameter file}
\seealso{\code{\link{make.par.file}}, \code{\link{run.analysis}}}
\examples{}