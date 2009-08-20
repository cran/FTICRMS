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
    \item{root.dir}{string containing directory of parameters file to be extracted from}
}

\details{Used by \code{\link{run.analysis}} to record all the parameter choices in an analysis for future reference.}
\value{
A list with the following components:
    \item{add.norm}{logical; whether to normalize additively or multiplicatively on the log scale}
    \item{add.par}{additive parameter for \code{"shiftedlog"} or \code{"glog"} options for \code{trans.method}}
    \item{align.fcn}{function (and inverse) to apply to masses before (and after) applying \code{align.method}}
    \item{align.method}{alignment algorithm for peaks}
    \item{base.dir}{directory for baseline files}
    \item{bhbysubj}{logical; whether to look for number of large peaks by subject (i.e., combining replicates) or by spectrum}
    \item{calc.all.peaks}{whether to calculate all possible peaks or only sufficiently large ones}
    \item{cluster.constant}{parameter used in running \code{cluster.method}}
    \item{cluster.method}{method for determining when two peaks from different spectra are the same}
    \item{cor.thresh}{threshhold correlation for declaring isotopes}
    \item{covariates}{data frame containing covariates used in analysis}
    \item{FDR}{False Discovery Rate in Benjamini-Hochberg test}
    \item{FTICRMS.version}{Version of \pkg{\link{FTICRMS}} that created file}
    \item{form}{formula used in \code{use.model}}
    \item{gengamma.quantiles}{whether to use generalized gamma quantiles when calculating large peaks}
    \item{halve.search}{whether to use a halving-line search if step leads to smaller value of function}
    \item{isotope.dist}{maximum distance for declaring isotopes}
    \item{lrg.dir}{directory for significant peaks file}
    \item{lrg.file}{name of file for storing large peaks}
    \item{lrg.only}{whether to consider only peaks that have at least one \dQuote{large} peak; i.e., identified by \code{\link{run.lrg.peaks}}}
    \item{masses}{specific masses to test}
    \item{max.iter}{convergence criterion in baseline calculation}
    \item{min.spect}{minimum number of spectra necessary for peak to be used in \code{\link{run.analysis}}}
    \item{neg.div}{negativity divisor in baseline calculation}
    \item{neg.norm.by}{method for negativity penalty in baseline analysis}
    \item{norm.peaks}{which peaks to use in normalization}
    \item{norm.post.repl}{logical; whether to normalize after combining replicates}
    \item{normalization}{type of normalization to use on spectra before statistical analysis}
    \item{num.pts}{number of points needed for peak fitting}
    \item{oneside.min}{minimum number of points on each side of local maximum for peak fitting}
    \item{overwrite}{whether to replace existing files with new ones}
    \item{par.file}{string containing name of parameters file}
    \item{peak.dir}{directory for peak location files}
    \item{peak.method}{method for locating peaks}
    \item{peak.thresh}{threshold for declaring large peak}
    \item{pre.align}{shifts to apply before running \code{\link{run.strong.peaks}}}
    \item{pval.fcn}{function to calculate \emph{p}-values}
    \item{R2.thresh}{\eqn{R^2} value needed for peak fitting}
    \item{raw.dir}{directory for raw data files}
    \item{rel.conv.crit}{whether convergence criterion should be relative to size of current baseline estimate}
    \item{repl.method}{how to deal with replicates}
    \item{res.dir}{directory for result file}
    \item{res.file}{name for results file}
    \item{root.dir}{directory for parameters file and raw data directory}
    \item{sm.div}{smoothness divisor in baseline calculation}
    \item{sm.norm.by}{method for smoothness penalty in baseline analysis}
    \item{sm.ord}{order of derivative to penalize in baseline analysis}
    \item{sm.par}{smoothing parameter for baseline calculation}
    \item{subs}{subset of spectra to use for analysis}
    \item{subtract.base}{whether to subtract calculated baseline from spectrum}
    \item{tol}{convergence criterion in baseline calculation}
    \item{trans.method}{data transformation method}
    \item{use.model}{what model to apply to data}
    \item{zero.rm}{whether to replace zeros in spectra with average of surrounding values}
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

Benjamini, Y. and Hochberg, Y.  (1995) \dQuote{Controlling the false discovery 
rate: a practical and powerful approach to multiple testing.}  \emph{J. Roy. 
Statist. Soc. Ser. B}, \bold{57}:1, 289--300.

Xi, Y. and Rocke, D.M. (2008) \dQuote{Baseline Correction for NMR Spectroscopic 
Metabolomics Data Analysis}.  \emph{BMC Bioinformatics}, \bold{9}:324.
}

\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}

\note{
\code{do.call(make.par.file, extract.pars())} recreates the original parameter file

\code{align.method}, \code{cluster.method}, \code{neg.norm.by}, \code{normalization}, 
\code{peak.method}, \code{sm.norm.by}, and \code{trans.method} can be abbreviated.

See \code{\link{make.par.file}} for a summary of which programs use each of the
parameters in the list.
}

\seealso{
\code{\link{make.par.file}}, \code{\link{run.analysis}}
}

\examples{}
