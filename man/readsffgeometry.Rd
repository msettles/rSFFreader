\name{readSffGeometry}
\alias{readSffGeometry}

\title{Read in basic statistics about an sff file}

\description{
  \code{readSffGeometry} processes a sff file or set of sff files, counting the number of reads and the length of each read.
}

\usage{
readSffGeometry(filenames)
}


\arguments{
  \item{filenames}{File or files to be read in.}
}

\value{A \code{\linkS4class{list}} with two elements, \code{nReads} and \code{Read_Widths} }

\author{Matt Settles <msettles@uidaho.edu>}

\examples{
sffgeom <- readSffGeometry(system.file("extdata","Small454Test.sff",package="rSFFreader"))
sffgeom$nReads
sffgeom$Read_Widths
}