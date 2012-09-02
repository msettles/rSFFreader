\name{readsffgeometry}
\alias{readsffgeometry}

\title{Read in basic statistics about an sff file}

\description{
  \code{readsffgeometry} processes a sff file or set of sff files, counting the number of reads and the length of each read.
}

\usage{
readsffgeometry(filenames)
}


\arguments{
  \item{filenames}{File or files to be read in.}
}

\value{A \code{\linkS4class{list}} with two elements, \code{nReads} and \code{Read_Widths} }

\author{Matt Settles <msettles@uidaho.edu>}

\examples{
sffstats <- readsffgeometry(system.file("extdata","Small454Test.sff",package="rSFFreader"))
sffstats$nReads
sffstats$Read_Widths
}