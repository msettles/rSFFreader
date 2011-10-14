\name{SffReadsQ-class}
\docType{class}

\alias{SffReadsQ-class}
\alias{SffReadsQ}
\alias{SffReadsQ,DNAStringSet,QualityScore,BStringSet-method}
\alias{SffReadsQ,DNAStringSet,QualityScore,missing-method}
\alias{SffReadsQ,DNAStringSet,BStringSet,BStringSet-method}
\alias{SffReadsQ,DNAStringSet,BStringSet,missing-method}
\alias{SffReadsQ,missing,missing,missing-method}
\alias{coerce,SffReadsQ,QualityScaledDNAStringSet-method}
\alias{writeFastq,SffReadsQ-method}
\alias{alphabetByCycle,SffReadsQ-method}
\alias{alphabetScore,SffReadsQ-method}
\alias{[,SffReadsQ,ANY,ANY-method}
\alias{[,SffReadsQ,ANY,missing-method}
\alias{[,SffReadsQ,missing,ANY-method}
\alias{[,SffReadsQ,missing,missing-method}
\alias{[,SffReadsQ,missing,missing,ANY-method}
\alias{[,SffReadsQ,missing,ANY,ANY-method}
\alias{[,SffReadsQ,ANY,ANY,ANY-method}
\alias{[,SffReadsQ,ANY,missing,ANY-method}
\alias{append,SffReadsQ,SffReadsQ,missing-method}
\alias{narrow,SffReadsQ-method}
\alias{trimTailw,SffReadsQ-method}
\alias{trimTails,SffReadsQ-method}
\alias{qa,SffReadsQ-method}
\alias{detail,SffReadsQ-method}

\title{"SffReadsQ" class for short reads and their quality scores}

\description{

  This class provides a way to store and manipulate, in a coordinated
  fashion, the reads, identifiers, and quality scores of uniform-length
  short reads.

}

\section{Objects from the Class}{

  Objects from this class are the result of \code{\link{readFastq}}, or
  can be constructed from \code{DNAStringSet}, \code{QualityScore}, and
  \code{BStringSet} objects, as described below.

}
\section{Slots}{

  Slots \code{sread} and \code{id} are inherited from
  \code{\linkS4class{SffReads}}. An additional slot defined in this
  class is:
  \describe{
    \item{\code{quality}:}{Object of class \code{"BStringSet"}
      representing a quality score (see \code{\link{readFastq}} for some
      discussion of quality score).}
  }
}
\section{Extends}{
Class \code{"\linkS4class{SffReads}"}, directly.
Class \code{"\linkS4class{.SffReadsBase}"}, by class "SffReads", distance 2.
}
\section{Methods}{

  Constructors include:

  \describe{

    \item{SffReadsQ}{\code{signature(sread = "DNAStringSet", quality = "QualityScore", id = "BStringSet")}:}
    \item{SffReadsQ}{\code{signature(sread = "DNAStringSet", quality = "BStringSet", id = "BStringSet")}:

      Create a \code{SffReadsQ} object from reads, their quality
      scores, and identifiers. When \code{quality} is of class
      \code{BStringSet}, the type of encoded quality score is inferred
      from the letters used in the scores. The length of \code{id} and
      \code{quality} must match that of \code{sread}.}

    \item{SffReadsQ}{\code{signature(sread = "DNAStringSet", quality = "QualityScore", id = "missing")}:}
    \item{SffReadsQ}{\code{signature(sread = "DNAStringSet", quality = "BStringSet", id = "missing")}:

      Create a \code{SffReadsQ} object from reads and their quality
      scores, creating empty identifiers. When \code{quality} is of
      class \code{BStringSet}, the type of encoded quality score is
      inferred from the letters used in the scores. }

    \item{SffReadsQ}{\code{signature(sread = "missing", quality = "missing", id = "missing", ...)}:
      Create an empty \code{SffReadsQ} object.}

  }

  See \code{\link{accessors}} for additional functions to access slot
  content, and \code{\linkS4class{SffReads}} for inherited
  methods. Additional methods include:

  \describe{

    \item{quality}{inherited from \code{signature(object = "ANY")}:
      access the quality slot of \code{object}.}


    \item{coerce}{\code{signature(from = "SFastqQuality", to = "QualityScaledDNAStringSet")}:
      
      (Use \code{as(from, "QualityScaledDNAStringSet")}) coerce objects
      of class \code{from} to class \code{to}, using the quality
      encoding implied by \code{quality(from)}. See
      \code{\linkS4class{QualityScore}} for supported quality
      classes and their coerced counterparts.}

    \item{writeFastq}{\code{signature(object = "SffReadsQ", file = "character", mode="character", ...)}:
      Write \code{object} to \code{file} in fastq format. \code{mode}
      defaults to \sQuote{w}. This creates a new file, or fails if
      \code{file} already exists. Use \code{mode="a"} to append to an
      existing file. \code{file} is expanded using \code{path.expand}.}

    \item{[}{\code{signature(x = "SffReadsQ", i = "ANY", j = "missing")}: 
      This method creates a new \code{SffReadsQ} object containing only
      those reads indexed by \code{i}. Additional methods on
      \sQuote{[,SffReads} do not provide additional functionality, but
      are present to limit inappropriate use.}

    \item{append}{\code{signature(x = "SffReadsQ", values = "SffReads", length = "missing")}:
      append the \code{sread}, \code{quality} and \code{id} slots of
      \code{values} after the corresponding fields of \code{x}.}

    \item{narrow}{\code{signature(x = "SffReadsQ", start = NA, end =
      NA, width = NA, use.names = TRUE)}: narrow \code{sread} and
      \code{quality} so that sequences are between \code{start} and
      \code{end} bases, according to
      \code{\link[IRanges:Ranges-utils]{narrow}} in the \code{IRanges}
      package.}

    \item{trimTailw}{\code{signature(object="SffReadsQ", k="integer",
	a="character", halfwidth="integer", ..., ranges=FALSE)}:
      trim trailing nucleotides when a window of width 2 * halfwidth + 1
      contains \code{k} or more quality scores falling at or below
      \code{a}.}

    \item{trimTails}{\code{signature(object="SffReadsQ", k="integer",
	a="character", successive=FALSE, ..., ranges=FALSE)}: trim
      trailing nucleotides if \code{k} nucleotides fall below the
      quality encoded by \code{a}. If \code{successive=FALSE}, the k'th
      failing nucleotide and all subseqent nucleotides are trimmed. If
      \code{successive=TRUE}, failing nucleotides must occur
      successively; the sequence is trimmed from the first of the
      successive failing nucleotides.}

    \item{alphabetByCycle}{\code{signature(stringSet = "SffReadsQ")}:
      Apply \code{\link{alphabetByCycle}} to the \code{sread} component,
      the \code{quality} component, and the combination of these two
      components of \code{stringSet}, returning a list of matrices with
      three elements: \code{"sread"}, \code{"quality"}, and
      \code{"both"}.}

    \item{alphabetScore}{\code{signature(object = "SffReadsQ")}:
      See \code{\link{alphabetScore}} for details.}

    \item{qa}{\code{signature(dirPath = "SffReadsQ", lane="character",
	..., verbose=FALSE)}:
      Perform quality assessment on the \code{SffReadsQ} object using
      \code{lane} to identify the object and returning an instance of
      \code{\linkS4class{SffReadsQQA}}. See \code{\link{qa}}}

    \item{detail}{\code{signature(x = "SffReadsQ")}: display the
      first and last entries of each of \code{sread}, \code{id}, and
      \code{quality} entries of \code{object}.}
  }
}

\author{Matt Settles with lots of code 'borrowed' from Martin Morgan's ShortRead Package}

\seealso{

  \code{\link{readFastq}} for creation of objects of this class from
  fastq-format files.

}
\examples{
showClass("SffReadsQ")
showMethods(class="SffReadsQ", where=getNamespace("SffReads"),
            inherit=FALSE)
showMethods(class="SffReads", where=getNamespace("SffReads"),
            inherit=FALSE)

sp <- SolexaPath(system.file('extdata', package='SffReads'))
rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
quality(rfq)
quality(trimTails(rfq, 2, "H", successive=TRUE))
}
\keyword{classes}
