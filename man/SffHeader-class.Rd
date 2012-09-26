\name{SffHeader-class}
\Rdversion{1.1}
\docType{class}

%Class
\alias{class:SffHeader}
\alias{SffHeader-class}
\alias{SffHeader}

\alias{show,SffHeader-method}
\alias{detail,SffHeader-method}

%Methods
\alias{header,SffHeader-method}
\alias{header}

\title{Class \code{"SffHeader"} contains top header information from the SFF files read in.}

\description{
This class contains a list-like structure with header descriptions
derived from one or more SFF files.

Objects from this class are the result of \code{\link{readSffHeader}},
or from the result of \code{\link{readSff}}.

}

\section{Slots}{
  \describe{
    \item{\code{header}:}{Object of class \code{"list"}, containing
      data frames or lists of data frames summarizing a description
      of the SFF files.}
  }
}

\section{Methods}{
  \describe{
    \item{header}{\code{signature(object = "SffHeader")}: access the 
    header slot of \code{object}, returning a list object.}
   }
}
\author{Matt Settles <msettles@uidaho.edu>}

\seealso{
  \code{\linkS4class{SffReads}},\code{\linkS4class{SffReadsQ}}.
}
\examples{
showClass("SffHeader")
header <- readSffHeader(system.file("extdata","Small454Test.sff",package="rSFFreader"))
header
}
\keyword{classes}
