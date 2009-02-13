\name{make.par.file}
\alias{make.par.file}
\title{Create Parameter File for FT-ICR MS Analysis}
\description{
Creates a file of parameters that can be read by the functions in the \code{FTICRMS} package
}
\usage{
make.par.file(covariates, form, par.file = "parameters.RData", root.dir = ".", \dots)
}
\arguments{
    \item{covariates}{data frame with rownames given by raw data files}
    \item{form}{object of class \dQuote{\code{\link{formula}}} to be used for testing using \code{covariates}}
    \item{par.file}{string containing name of file}
    \item{root.dir}{string containing location for file}
    \item{\dots}{parameters whose default values are to be overwritten (see below)}
}
\details{
Creates a file with name given by \code{par.file} in directory given by \code{root.dir} which
contains values for all of the parameters used in the programs in the \code{FTICRMS} package.
The possible parameters that can be included in \code{\dots}, their default values, their
descriptions, and the program(s) in which they are used are as follows:
\tabular{lll}{ \tab \tab \cr
    \code{add.norm = TRUE}\tab logical; whether to normalize additively or multiplicatively on the log scale\tab \code{\link{run.analysis}}\cr
    \code{add.par = 0}\tab additive parameter for \code{"shiftedlog"} or \code{"glog"} options for \code{trans.method}\tab \code{\link{run.cluster.matrix}}, \code{\link{run.lrg.peaks}}, \code{\link{run.peaks}}\cr
    \code{align.method = "spline"}\tab alignment algorithm for peaks\tab \code{\link{run.cluster.matrix}}\cr
    \code{base.dir = paste(root.dir, "/Baselines", sep = "")}\tab directory for baseline files\tab \code{\link{run.baselines}}, \code{\link{run.cluster.matrix}}, \code{\link{run.lrg.peaks}}, \code{\link{run.peaks}}\cr
    \code{calc.all.peaks = FALSE}\tab logical; whether to calculate all possible peaks or only sufficiently large ones\tab \code{\link{run.cluster.matrix}}, \code{\link{run.lrg.peaks}}, \code{\link{run.peaks}}\cr
    \code{cluster.constant = 10}\tab NA\tab \code{\link{run.cluster.matrix}}\cr
    \code{cluster.method = "ppm"}\tab NA\tab \code{\link{run.cluster.matrix}}\cr
    \code{cor.thresh = 0.8}\tab threshold correlation for declaring isotopes\tab \code{\link{run.strong.peaks}}\cr
    \code{FDR = 0.1}\tab False Discovery Rate in Benjamini-Hochberg test\tab \code{\link{run.analysis}}\cr
    \code{gengamma.quantiles = TRUE}\tab logical; whether to use generalized gamma quantiles when calculating large peaks\tab \code{\link{run.lrg.peaks}}, \code{\link{run.peaks}}\cr
    \code{halve.search = FALSE}\tab logical; whether to use a halving-line search if step leads to smaller value of function\tab \code{\link{run.baselines}}\cr
    \code{isotope.dist = 7}\tab maximum distance for declaring isotopes\tab \code{\link{run.analysis}}, \code{\link{run.cluster.matrix}}, \code{\link{run.strong.peaks}}\cr
    \code{lrg.dir = paste(root.dir, "/Large_Peaks", sep = "")}\tab directory for large peaks file\tab \code{\link{run.analysis}}, \code{\link{run.cluster.matrix}}, \code{\link{run.lrg.peaks}}, \code{\link{run.strong.peaks}}\cr
    \code{lrg.file = "lrg_peaks.RData"}\tab name of file for storing large peaks\tab \code{\link{run.analysis}}, \code{\link{run.cluster.matrix}}, \code{\link{run.lrg.peaks}}, \code{\link{run.strong.peaks}}\cr
    \code{lrg.only = TRUE}\tab logical; whether to consider only peaks that have at least one \dQuote{large}peak; i.e., identified by \code{run.lrg.peaks}\tab \code{\link{run.analysis}}, \code{\link{run.cluster.matrix}}\cr
    \code{masses = NULL}\tab specific masses to test\tab \code{\link{run.analysis}}, \code{\link{run.cluster.matrix}}\cr
    \code{max.iter = 40}\tab convergence criterion in baseline calculation\tab \code{\link{run.baselines}}\cr
    \code{neg.div = NA}\tab negativity divisor in baseline calculation\tab \code{\link{run.baselines}}\cr
    \code{neg.norm.by = c("baseline", "overestimate", "constant")}\tab method for negativity penalty in baseline analysis\tab \code{\link{run.baselines}}\cr
    \code{normalization = "common"}\tab type of normalization to use on spectra before statistical anlaysis\tab \code{\link{run.analysis}}\cr
    \code{num.pts = 5}\tab number of consecutive points needed for peak fitting\tab \code{\link{run.cluster.matrix}}, \code{\link{run.peaks}}\cr
    \code{oneside.min = 1}\tab minimum number of points on each side of local maximum for peak fitting\tab \code{\link{run.cluster.matrix}}, \code{\link{run.peaks}}\cr
    \code{overwrite = FALSE}\tab whether to replace exisiting files with new ones\tab All six programs\cr
    \code{par.file = "parameters.RData"}\tab string containing name of parameters file\tab All six programs\cr
    \code{peak.dir = paste(root.dir, "/All_Peaks", sep = "")}\tab directory for peak location files\tab \code{\link{run.cluster.matrix}}, \code{\link{run.lrg.peaks}}, \code{\link{run.peaks}}\cr
    \code{peak.method = "parabola"}\tab method for locating peaks\tab \code{\link{run.cluster.matrix}}, \code{\link{run.peaks}}\cr
    \code{peak.thresh = 3.798194}\tab threshold for declaring large peak\tab \code{\link{run.lrg.peaks}}, \code{\link{run.peaks}}\cr
    \code{pre.align = FALSE}\tab shifts to apply before running \code{\link{run.strong.peaks}}\tab \code{\link{run.cluster.matrix}}, \code{\link{run.strong.peaks}}\cr
    \code{pval.fcn = "default"}\tab function to calculate \emph{p}-values if \code{use.t.test = FALSE}\tab \code{\link{run.analysis}}\cr
    \code{R2.thresh = 0.98}\tab \eqn{R^2} value needed for peak fitting\tab \code{\link{run.cluster.matrix}}, \code{\link{run.peaks}}\cr
    \code{raw.dir = paste(root.dir, "/Raw_Data", sep = "")}\tab directory for raw data files\tab \code{\link{run.baselines}}\cr
    \code{rel.conv.crit = TRUE}\tab whether convergence criterion should be relatiev to current baseline estimate\tab \code{\link{run.baselines}}\cr
    \code{repl.method = max}\tab how to deal with replicates\tab \code{\link{run.analysis}}\cr
    \code{res.dir = paste(root.dir, "/Results", sep = "")}\tab directory for result file\tab \code{\link{run.analysis}}\cr
    \code{res.file = "analyzed.RData"}\tab name for results file\tab \code{\link{run.analysis}}\cr
    \code{root.dir = "."}\tab directory for parameters file and raw data\tab All six programs\cr
    \code{sm.div = NA}\tab smoothness divisor in baseline calculation\tab \code{\link{run.baselines}}\cr
    \code{sm.norm.by = c("baseline", "overestimate", "constant")}\tab method for smoothness penalty in baseline analysis\tab \code{\link{run.baselines}}\cr
    \code{sm.ord = 2}\tab order of derivative to kill in baseline analysis\tab \code{\link{run.baselines}}\cr
    \code{sm.par = 1e-11}\tab smoothing parameter for baseline calculation\tab \code{\link{run.baselines}}\cr
    \code{subtract.base = FALSE}\tab logical; whether to subtract calculated baseline from spectrum\tab \code{\link{run.cluster.matrix}}, \code{\link{run.lrg.peaks}}, \code{\link{run.peaks}}\cr
    \code{tol = 5e-8}\tab convergence criterion in baseline calculation\tab \code{\link{run.baselines}}\cr
    \code{trans.method = "shiftedlog"}\tab data transformation method\tab \code{\link{run.cluster.matrix}}, \code{\link{run.lrg.peaks}}, \code{\link{run.peaks}}\cr
    \code{use.t.test = FALSE}\tab whether to use a \emph{t}-test to calculate \emph{p}-values\tab \code{\link{run.analysis}}\cr
    \code{zero.rm = TRUE}\tab whether to replace zeros in spectra with average of surrounding values\tab \code{\link{run.baselines}}\cr
}
}
\value{
No value returned; the file \code{par.file} is simply created in \code{root.dir}.
}
\references{
Barkauskas, D.A. (2009) \dQuote{Statistical Analysis of Matrix-Assisted Laser Desorption/Ionization
Fourier Transform Ion Cyclotron Resonance Mass Spectrometry Data with Applications to Cancer
Biomarker Detection}.  Ph.D. dissertation, University of California at Davis.

Barkauskas, D.A. \emph{et al}. (2009) \dQuote{Detecting glycan cancer biomarkers in serum
samples using MALDI FT-ICR mass spectrometry data}.  \emph{Bioinformatics}, \bold{25}:2, 251--257.

Xi, Y. and Rocke, D.M. (2008) \dQuote{Baseline Correction for NMR Spectroscopic Metabolomics Data Analysis}.  \emph{BMC Bioiniformatics}, \bold{9}:324.
}
\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}
\note{
\code{do.call(make.par.file, extract.pars())} recreates the original parameter file.

See the individual function help pages for each function for more detailed descriptions of the above parameters.
}
\seealso{\code{\link{extract.pars}}}
\examples{}
