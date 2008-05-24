\name{make.par.file}
\alias{make.par.file}
\title{Create Parameter File for FT-ICR MS Analysis}
\description{
Creates a file of parameters that can be read by the functions in the \code{FTICRMS} package
}
\usage{
make.par.file(covariates, form, root.dir = ".", par.file = "parameters.RData", \dots)
}
\arguments{
  \item{covariates}{data frame with rownames given by raw data files}
  \item{form}{object of class \dQuote{\code{\link{formula}}} to be used for testing using \code{covariates}}
  \item{root.dir}{string containing location for file}
  \item{par.file}{string containing name of file}
  \item{\dots}{parameters whose default values are to be overwritten (see below)}
}
\details{
Creates a file called \code{parameters.RData} in directory \code{root.dir} which contains
values for all of the parameters used in the programs in the \code{FTICRMS} package.  The 
possible parameters, their default values, and their descriptions are as follows:
\tabular{ll}{ \tab \cr
        \code{add.norm = TRUE}\tab logical; whether to normalize additively or multiplicatively on the log scale\cr
        \code{add.par = 10}\tab additive parameter for \code{"shiftedlog"} or \code{"glog"} options for \code{trans.method}\cr
        \code{align.method = "spline"}\tab alignment algorithm for peaks\cr
        \code{base.dir = paste(root.dir, "/Baseline_Corrected", sep="")}\tab directory for baseline-corrected files\cr
        \code{binsize = 64}\tab size of bins used in estimating noise variance in baseline calculation\cr
        \code{cor.thresh = 0.8}\tab threshhold correlation for declaring isotopes\cr
        \code{FDR = 0.1}\tab False Discovery Rate in Benjamini-Hochberg test\cr
        \code{frac.changed = 0.001}\tab convergence criterion in baseline calculation\cr
        \code{isotope.dist = 7}\tab maximum distance for declaring isotopes\cr
        \code{masses = NULL}\tab specific masses to test\cr
        \code{max.iter = 30}\tab convergence criterion in baseline calculation\cr
        \code{neg.pen = sqrt(pi/2)}\tab negativity penalty in baseline calculation\cr
        \code{normalization = "common"}\tab type of normalization to use on spectra before statistical anlaysis\cr
        \code{num.pts = 5}\tab number of points needed for peak fitting\cr
        \code{numsds = 4}\tab number of standard deviations above the mean to be declared non-noise\cr
        \code{oneside.min = 1}\tab minimum number of points on each side of local maximum for peak fitting\cr
        \code{overwrite = FALSE}\tab whether to replace exisiting files with new ones\cr
        \code{peak.dir = paste(root.dir, "/All_Peaks", sep="")}\tab directory for peak location files\cr
        \code{peak.method = "parabola"}\tab method for locating peaks\cr
        \code{pre.align = 0}\tab shifts to apply before running \code{\link{run.strong.peaks}}\cr
        \code{R2.thresh = 0.98}\tab \eqn{R^2} value needed for peak fitting\cr
        \code{raw.dir = paste(root.dir, "/Raw_Data", sep="")}\tab directory for raw data files\cr
        \code{repl.method = max}\tab how to deal with replicates\cr
        \code{res.dir = paste(root.dir, "/Results", sep="")}\tab directory for result file\cr
        \code{root.dir = "."}\tab directory for parameters file and raw data\cr
        \code{lrg.dir = paste(root.dir, "/Large_Peaks", sep="")}\tab directory for significant peaks file\cr 
        \code{lrg.only = TRUE}\tab whether to consider only peaks that have at least one peak \dQuote{significant}; i.e.,
            identified by \code{run.lrg.peaks} \cr
        \code{sm.fac = 10^15}\tab smoothing factor for baseline calculation\cr
        \code{trans.method = "shiftedlog"}\tab data transformation method\cr
        \code{use.t.test = TRUE}\tab whether to use a \emph{t}-test or an ANOVA/regression model \cr
}
}
\value{
No value returned; the file \code{par.file} is simply created in \code{root.dir}.
}
\references{
Barkauskas, D.A. \emph{et al}. (2008) \dQuote{Detecting glycan cancer biomarkers in serum 
samples using MALDI FT-ICR mass spectrometry data}.  Submitted to \emph{Bioinformatics}

Xi, Y. and Rocke, D.M. (2006) \dQuote{Baseline Correction for NMR Spectroscopic Metabolomics Data Analysis}, unpublished.
}
\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}
\note{
\code{do.call(make.par.file, extract.pars())} recreates the original parameter file.

See the individual function help pages for each function for more detailed descriptions of the above parameters.
}
\seealso{\code{\link{extract.pars}}}
\examples{}
