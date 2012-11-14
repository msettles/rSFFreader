\name{readSffHeader}
\alias{readSffHeader}

\title{Read in sff header}

\description{
  \code{readSffHeader} Processes an sff file or set of sff files, and store the header information in a \code{\linkS4class{SffHeader}} object.
}

\usage{
readSffHeader(filenames, verbose=TRUE)
}


\arguments{
  \item{filenames}{File or files to be read in.}
  
  \item{verbose}{\code{logical(1)} indicating whether to return verbose output.}
  
}

\value{A \code{\linkS4class{SffHeader}} object.}


\author{Matt Settles <msettles@uidaho.edu>}

\examples{
sffhead <- readSffHeader(system.file("extdata","Small454Test.sff",package="rSFFreader"))
sffhead@header
}