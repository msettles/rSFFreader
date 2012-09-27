\name{SffReads-class}
\Rdversion{1.1}
\docType{class}
%Class
\alias{class:SffReads}
\alias{SffReads-class}
\alias{SffReads}

%Methods
\alias{show,SffReads-method}
\alias{detail,SffReads-method}
\alias{[,SffReads,ANY,ANY-method}
\alias{[,SffReads,ANY,missing-method}
\alias{[,SffReads,missing,ANY-method}
\alias{[,SffReads,missing,missing-method}
\alias{adapterClip,SffReads-method}
\alias{adapterClip<-,SffReads-method}
\alias{clipMode,SffReads-method}
\alias{clipMode<-,SffReads-method}
\alias{customClip,SffReads-method}
\alias{customClip<-,SffReads-method}
\alias{fullClip,SffReads-method}
\alias{length,SffReads-method}
\alias{names,SffReads-method}
\alias{names<-,SffReads,ANY-method}
\alias{qualityClip,SffReads-method}
\alias{qualityClip<-,SffReads-method}
\alias{rawClip,SffReads-method}
%\alias{sread,ANY-method}
\alias{sread}
\alias{id,SffReads-method}
\alias{name,SffReads-method}
\alias{append,SffReads,SffReads,missing-method}
\alias{width,SffReads-method}
\alias{writeFasta,SffReads-method}
\alias{solveSffSEW}
\alias{adapterClip}
\alias{adapterClip<-}
\alias{clipMode}
\alias{clipMode<-}
\alias{customClip}
\alias{customClip<-}
\alias{fullClip}
\alias{qualityClip}
\alias{qualityClip<-}
\alias{rawClip}



\title{Class \code{"SffReads"}. A class for next generation reads stored in
SFF files (Roche 454 and Life Sciences Ion Torrent)}

\description{
  This class provides a way to store and manipulate, in a coordinated
  fashion, next generation reads stored in SFF files and their identifiers.

  Objects from this class are created by \code{readSff}, or by calls to the
  constructor \code{SffReads}, as outlined below.
}
\section{Slots}{

    Slots \code{header} are inherited from \code{\linkS4class{SffHeader}}. 
    Additional slots defined in this class are:

  \describe{
    \item{\code{sread}:}{Object of class \code{"DNAStringSet"} containing
    IUPAC-standard, variable-length DNA strings representing sequence reads. }
    \item{\code{qualityIR}:}{Object of class \code{"IRanges"} contains the 
    clip points associated with quality clip locations specified in the sff 
    file read header. Clip locations are auto filled to those specified in
    the sff file. }
    \item{\code{adapterIR}:}{Object of class \code{"IRanges"} contains the 
    clip points associated with adapter clip locations specified in the sff 
    file read header. Clip locations are auto filled to those specified in 
    the sff file. }
    \item{\code{customIR}:}{Object of class \code{"IRanges"} contains the 
    clip points associated with custom clip locations as specified by the 
    user. Initial object is set to empty. }
    \item{\code{clipMode}:}{Object of class \code{"character"} contains the
    current active clipMode, one of "raw", "adapter", "quality", "custom", 
    "full". See explaination below for an expanded desription of clip modes.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SffHeader}"}, directly.
}

\section{Contructors}{
  \describe{
      \item{SffReads}{\code{sread = "DNAStringSet", id = "BStringSet"}:
      Create a \code{SffReads} object from reads and their
      identifiers. The length of \code{id} must match that of \code{sread}.}
    \item{SffReads}{\code{sread = "DNAStringSet", id = "missing"}:
      Create a \code{SffReads} object from reads, creating empty identifiers.}
    \item{SffReads}{\code{sread = "missing", id = "missing", ...}:
      Create an empty \code{SffReads} object.}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "SffReads", i = "ANY", j = "missing")}:
      This method creates a new \code{SffReads} object containing only
      those reads indexed by \code{i}. Additional methods on
      \sQuote{[,SffReads} do not provide additional functionality, but
      are present to limit appropriate use.}
    \item{adapterClip}{\code{signature(object = "SffReads")}: Returns an IRanges
      object of clip points when clipMode is set to "adapter"}
    \item{adapterClip<-}{\code{signature(object = "SffReads")}: Provide new adapter
      clip points with an IRanges object.}
    \item{clipMode}{\code{signature(object = "SffReads")}: returns the current
    active clipMode being used.}
    \item{clipMode<-}{\code{signature(object = "SffReads")}: set the current clip mode
    to one of "raw", "adapter", "quality", "full", or "custom"}
    \item{customClip}{\code{signature(object = "SffReads")}: Returns an IRanges
      object of clip points when clipMode is set to "custom" }
    \item{customClip<-}{\code{signature(object = "SffReads")}:Provide new custom
      clip points with an IRanges object.}
    \item{fullClip}{\code{signature(object = "SffReads")}: Returns an IRanges
      object of clip points when clipMode is set to "full" }
    \item{length}{\code{signature(x = "SffReads")}: returns a
      \code{integer(1)} vector describing the number of reads in this object.}
    \item{width}{\code{signature(x = "SffReads")}: returns an
      \code{integer()} vector of the widths of each read in this object.}
    \item{names}{\code{signature(x = "SffReads")}: access the
      id slot of the object returning a character vector of read names.  }
    \item{names<-}{\code{signature(x = "SffReads", value = "ANY")}: Provide new read
      ids with a \code{character()} vector of same length as the object.}
    \item{qualityClip}{\code{signature(object = "SffReads")}: Returns an IRanges
      object of clip points when clipMode is set to "quality". }
    \item{qualityClip<-}{\code{signature(object = "SffReads")}: Provide new quality
      clip points with an IRanges object.}
    \item{rawClip}{\code{signature(object = "SffReads")}: Returns an IRanges
      object of clip points when clipMode is set to "full". }
    \item{sread}{\code{signature(object = "ANY")}: access the
      sread slot of \code{object}, returning a DNAStringSet object of reads.}
    \item{id}{\code{signature(object = "SffReads")}: this function is included to maintain compatibilty with 
      the \linkS4class{ShortRead} class, returning a BStringSet object of read names for the \code{object}.}
    \item{append}{\code{signature(x = "SffReads", values = "SffReads", length = "missing")}:
      append the \code{sread} and \code{id} slots of \code{values} after
      the corresponding fields of \code{x}.}
    \item{writeFasta}{\code{signature(object="SffReads", file, ...)}: write
      \code{object} to \code{file} in fasta format. See
      \code{\link{writeFasta}} for \code{...} argument values.}
  }
}

\author{Matt Settles with lots of code 'borrowed' from Martin Morgan's ShortRead Package}
\seealso{
  \code{\linkS4class{SffReadsQ} and \link{availableClipModes}}
}
\examples{
showClass("SffReads")
showMethods(class="SffReads", where=getNamespace("rSFFreader"),
            inherit=FALSE)
            
sff <- readSff(system.file("extdata","Small454Test.sff",package="rSFFreader"), use.qualities=FALSE)
detail(sff)

}
\keyword{classes}
