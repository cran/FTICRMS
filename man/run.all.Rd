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
The analysis described in Barkauskas \emph{et al}. (2008) can be (approximately)
reproduced using the following parameter values instead of the defaults:
\tabular{l}{
    \code{add.par = 10}\cr
    \code{calc.all.peaks = TRUE}\cr
    \code{gengamma.quantiles = FALSE}\cr
    \code{max.iter = 30}\cr
    \code{neg.norm.by = "constant"}\cr
    \code{peak.thresh = 4}\cr
    \code{pval.fcn = function(x){anova(x)[2,5]}}\cr
    \code{rel.conv.crit = FALSE}\cr
    \code{sm.norm.by = "constant"}\cr
    \code{subtract.base = TRUE}\cr
    \code{zero.rm = FALSE}\cr
}
(It is only an approximate reproduction because the stopping criterion for
baseline calculation used in the article turned out to be a poor one and is no
longer supported in the package.  This shouldn't make a very large difference,
however.)
}

\seealso{
\code{\link{make.par.file}}, \code{\link{run.baselines}},
\code{\link{run.peaks}}, \code{\link{run.lrg.peaks}},
\code{\link{run.strong.peaks}}, \code{\link{run.cluster.matrix}},
\code{\link{run.analysis}}
}

\examples{}
