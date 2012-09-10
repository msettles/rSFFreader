\name{SffHeader-class}
\docType{class}
\alias{SffHeader-class}
\alias{header,SffHeader-method}
\alias{show,SffHeader-method}
\alias{detail,SffHeader-method}

\title{SFF Header information}

\description{

  This class contains a list-like structure with header descriptions
  derived from one or more SFF files.

}

\section{Objects from the Class}{

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
Class \code{"\linkS4class{.SffBase}"}, directly.
}

\section{Methods}{
  \describe{
  \item{header}{inherited from \code{signature(object = "ANY")}:
    access the header slot of \code{object}.}
  }
}
\author{Matt Settles <msettles@uidaho.edu>}

\seealso{
  \code{\link{SffReads}}.
}
\examples{
showClass("SffHeader")
}
\keyword{classes}
