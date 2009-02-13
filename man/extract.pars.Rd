\name{extract.pars}
\alias{extract.pars}
\title{Extract Parameters from File}
\description{
Extracts the parameters in the file specified by \code{par.file} and returns them in list form.
}
\usage{
extract.pars(par.file = "parameters.RData", root.dir = ".")
}
\arguments{
    \item{par.file}{string containing name of parameters file}
    \item{root.dir}{string containing location of parameters file to be extracted from}
}
\details{Used by \code{\link{run.analysis}} to record all the parameter choices in an analysis for future reference.}
\value{
A list with the following components:
    \item{add.norm}{logical; whether to normalize additively or multiplicatively on the log scale}
    \item{add.par}{additive parameter for \code{"shiftedlog"} or \code{"glog"} options for \code{trans.method}}
    \item{align.method}{alignment algorithm for peaks}
    \item{base.dir}{directory for baseline-corrected files}
    \item{calc.all.peaks}{whether to calculate all possible peaks or only sufficiently large ones}
    \item{cluster.constant}{parameter used in running \code{cluster.method}}
    \item{cluster.method}{method for determining when two peaks from different spectra are the same}
    \item{cor.thresh}{threshhold correlation for declaring isotopes}
    \item{covariates}{data frame containing covariates used in analysis}
    \item{FDR}{False Discovery Rate in Benjamini-Hochberg test}
    \item{form}{formula used in \code{\link{t.test}} or \code{\link{lm}}}
    \item{gengamma.quantiles}{whether to use generalized gamma quantiles when calculating large peaks}
    \item{halve.search}{whether to use a halving-line search if step leads to smaller value of function}
    \item{isotope.dist}{maximum distance for declaring isotopes}
    \item{lrg.dir}{directory for significant peaks file}
    \item{lrg.file}{name of file for storing large peaks}
    \item{lrg.only}{whether to consider only peaks that have at least one peak \dQuote{large}; i.e., identified by \code{\link{run.lrg.peaks}}}
    \item{masses}{specific masses to test}
    \item{max.iter}{convergence criterion in baseline calculation}
    \item{neg.div}{negativity divisor in baseline calculation}
    \item{neg.norm.by}{method for negativity penalty in baseline analysis}
    \item{normalization}{type of normalization to use on spectra before statistical analysis}
    \item{num.pts}{number of points needed for peak fitting}
    \item{oneside.min}{minimum number of points on each side of local maximum for peak fitting}
    \item{overwrite}{whether to replace existing files with new ones}
    \item{par.file}{string containing name of parameters file}
    \item{peak.dir}{directory for peak location files}
    \item{peak.method}{method for locating peaks}
    \item{peak.thresh}{threshold for declaring large peak}
    \item{pre.align}{shifts to apply before running \code{\link{run.strong.peaks}}}
    \item{pval.fcn}{function to calculate \emph{p}-values if \code{use.t.test = FALSE}}
    \item{R2.thresh}{\eqn{R^2} value needed for peak fitting}
    \item{raw.dir}{directory for raw data files}
    \item{rel.conv.crit}{whether convergence criterion should be relatiev to current baseline estimate}
    \item{repl.method}{how to deal with replicates}
    \item{res.dir}{directory for result file}
    \item{res.file}{name for results file}
    \item{root.dir}{directory for parameters file and raw data directory}
    \item{sm.div}{smoothness divisor in baseline calculation}
    \item{sm.norm.by}{method for smoothness penalty in baseline analysis}
    \item{sm.ord}{order of derivative to kill in baseline analysis}
    \item{sm.par}{smoothing parameter for baseline calculation}
    \item{subtract.base}{whether to subtract calculated baseline from spectrum}
    \item{tol}{convergence criterion in baseline calculation}
    \item{trans.method}{data transformation method}
    \item{use.t.test}{whether to use a \emph{t}-test to calculate \emph{p}-values}
    \item{zero.rm}{whether to replace zeros in spectra with average of surrounding values}
}
\references{
Barkauskas, D.A. (2009) \dQuote{Statistical Analysis of Matrix-Assisted Laser Desorption/Ionization
Fourier Transform Ion Cyclotron Resonance Mass Spectrometry Data with Applications to Cancer
Biomarker Detection}.  Ph.D. dissertation, University of California at Davis.

Barkauskas, D.A. \emph{et al}. (2009) \dQuote{Detecting glycan cancer biomarkers in serum
samples using MALDI FT-ICR mass spectrometry data}.  \emph{Bioinformatics}, \bold{25}:2, 251--257.

Xi, Y. and Rocke, D.M. (2008) \dQuote{Baseline Correction for NMR Spectroscopic Metabolomics Data Analysis}.  \emph{BMC Bioinformatics}, \bold{9}:324.
}
\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}
\note{\code{do.call(make.par.file, extract.pars())} recreates the original parameter file}
\seealso{\code{\link{make.par.file}}, \code{\link{run.analysis}}}
\examples{}
