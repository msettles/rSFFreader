\name{SffHeader-class}
\Rdversion{1.1}
\docType{class}

%Class
\alias{class:SffHeader}
\alias{SffHeader-class}
\alias{SffHeader}

%Methods
\alias{header,SffHeader-method}
\alias{header}

\title{Class \code{"SffHeader"} contains top header information from the SFF files read in.}

\description{
Objects can be created by calls of the form \code{new("SffHeader", ...)}.
This class contains a list-like structure with header descriptions
derived from one or more SFF files.

Objects from this class are the result of \code{\link{readsffheader}},
or from the result of \code{\link{readsff}}.

}

\section{Slots}{
  \describe{
    \item{\code{header}:}{Object of class \code{"list"}, containing
      data frames or lists of data frames summarizing a description
      of the SFF files.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{.SFFBase}"}, directly.
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
}
\keyword{classes}
