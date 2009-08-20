\name{run.baselines}
\alias{run.baselines}
\title{Calculate and Store Baselines for Spectroscopic Data}

\description{
Takes the spectra from files in \code{raw.dir}, calculates the baselines from them,
and writes the results in the directory \code{base.dir}.
}

\usage{
run.baselines(root.dir = ".", raw.dir, base.dir, overwrite = FALSE,
              use.par.file = FALSE, par.file = "parameters.RData",
              sm.par = 1e-11, sm.ord = 2, max.iter = 20, tol = 5e-8,
              sm.div = NA, sm.norm.by = c("baseline", "overestimate", "constant"),
              neg.div = NA, neg.norm.by = c("baseline", "overestimate", "constant"),
              rel.conv.crit = TRUE, zero.rm = TRUE, halve.search = FALSE)
}

\arguments{
    \item{root.dir}{directory for parameters file and raw data}
    \item{raw.dir}{directory for raw data files; default is \code{paste(root.dir, "/Raw_Data", sep = "")}}
    \item{base.dir}{directory for baseline files; default is \code{paste(root.dir, "/Baselines", sep = "")}}
    \item{overwrite}{logical; whether to replace existing files with new ones}
    \item{use.par.file}{logical; if \code{TRUE}, then parameters are read from \code{par.file} in directory \code{root.dir}}
    \item{par.file}{string containing name of parameters file}
    \item{sm.par}{smoothing parameter for baseline calculation}
    \item{sm.ord}{order of derivative to penalize in baseline analysis}
    \item{max.iter}{convergence criterion in baseline calculation}
    \item{tol}{convergence criterion}
    \item{sm.div}{smoothness divisor in baseline calculation}
    \item{sm.norm.by}{method for smoothness penalty in baseline analysis}
    \item{neg.div}{negativity divisor in baseline calculation}
    \item{neg.norm.by}{method for negativity penalty in baseline analysis}
    \item{rel.conv.crit}{logical; whether convergence criterion should be relative to size of current baseline estimate}
    \item{zero.rm}{logical; whether to replace zeros with average of surrounding values}
    \item{halve.search}{logical; whether to use a halving-line search if step leads to smaller value of function}
}

\details{
Goes through the entire directory \code{raw.dir} file-by-file and computes each 
baseline using \code{\link{baseline}}, then writes the spectrum and the baseline 
to a file in directory \code{base.dir}.  The name of the new file is the same as 
the name of the old file with \dQuote{.txt} replaced by \dQuote{.RData}, and the 
new file is ready to be used by \code{\link{run.peaks}}.

The files in \code{raw.dir} must be in a specific format (future versions of the 
package will allow for more flexibility).  The files should be two-column text 
files with mass in the first column and spectrum intensity in the second column.  
There should be no header row (just start the file with the first data point).  
The columns can be either comma-separated or whitespace-separated and the 
program will automatically detect which each file is.  The decimal separator 
should be \code{"."}, as using \code{","} will cause errrors in reading the 
files.

See \code{\link{baseline}} for details of all the parameters after 
\code{par.file}.
}

\value{
No value returned; the files are simply created.
}

\references{
Barkauskas, D.A. (2009) \dQuote{Statistical Analysis of Matrix-Assisted Laser Desorption/Ionization
Fourier Transform Ion Cyclotron Resonance Mass Spectrometry Data with Applications to Cancer
Biomarker Detection}.  Ph.D. dissertation, University of California at Davis.

Barkauskas, D.A. \emph{et al}. (2009) \dQuote{Detecting glycan cancer biomarkers in serum
samples using MALDI FT-ICR mass spectrometry data}.  \emph{Bioinformatics}, \bold{25}:2, 251--257.

Xi, Y. and Rocke, D.M. (2008) \dQuote{Baseline Correction for NMR Spectroscopic Metabolomics Data Analysis}.
\emph{BMC Bioinformatics}, \bold{9}:324.
}

\author{Don Barkauskas (\email{barkda@wald.ucdavis.edu})}

\note{
If \code{use.par.file == TRUE} and other parameters are entered into the function
call, then the parameters entered in the function call overwrite those read in
from the file.  Note that this is opposite from the behavior for
\pkg{\link{FTICRMS}} versions 0.7 and earlier.

The values of \code{sm.norm.by} and \code{neg.norm.by} can be abbreviated and
both have default value \code{"baseline"}.
}

\seealso{
\code{\link{baseline}}, \code{\link{run.peaks}}
}

\examples{}
