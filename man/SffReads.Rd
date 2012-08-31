\name{SffReads-class}
\docType{class}
\alias{SffReads-class}
\alias{SffReads}
\alias{SffReads,DNAStringSet,BStringSet-method}
\alias{SffReads,DNAStringSet,missing-method}
\alias{SffReads,missing,missing-method}
\alias{srrank,SffReads-method}
\alias{srorder,SffReads-method}
\alias{srsort,SffReads-method}
\alias{srduplicated,SffReads-method}
\alias{srdistance,SffReads,ANY-method}
\alias{trimLRPatterns,SffReads-method}
\alias{alphabetByCycle,SffReads-method}
\alias{tables,SffReads-method}
\alias{clean,SffReads-method}
\alias{[,SffReads,ANY,ANY-method}
\alias{[,SffReads,ANY,missing-method}
\alias{[,SffReads,missing,ANY-method}
\alias{[,SffReads,missing,missing-method}
\alias{[,SffReads,missing,missing,ANY-method}
\alias{[,SffReads,missing,ANY,ANY-method}
\alias{[,SffReads,ANY,ANY,ANY-method}
\alias{[,SffReads,ANY,missing,ANY-method}
\alias{append,SffReads,SffReads,missing-method}
\alias{narrow,SffReads-method}
\alias{detail,SffReads-method}
\alias{id,SffReads-method}
\alias{width,SffReads-method}
\alias{length,SffReads-method}
\alias{writeFasta,SffReads-method}
\alias{show,SffReads-method}

\title{"SffReads" class for short reads}

\description{

  This class provides a way to store and manipulate, in a coordinated
  fashion, uniform-length short reads and their identifiers.

}

\section{Objects from the Class}{

  Objects from this class are created by \code{readFasta}, or by
  calls to the constructor \code{SffReads}, as outlined below.

}
\section{Slots}{
  \describe{
    \item{\code{sread}:}{Object of class \code{"DNAStringSet"}
      containing IUPAC-standard, uniform-length DNA strings represent
      short sequence reads.}
    \item{\code{id}:}{Object of class \code{"BStringSet"} containing
      identifiers, one for each short read.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{.SffReadsBase}"}, directly.
}
\section{Methods}{

  Constructors include:

  \describe{

    \item{SffReads}{\code{signature(sread = "DNAStringSet", id = "BStringSet")}:
      Create a \code{SffReads} object from reads and their
      identifiers. The length of \code{id} must match that of \code{sread}.}

    \item{SffReads}{\code{signature(sread = "DNAStringSet", id = "missing")}:
      Create a \code{SffReads} object from reads, creating empty identifiers.}

    \item{SffReads}{\code{signature(sread = "missing", id = "missing", ...)}:
      Create an empty \code{SffReads} object.}

  }

  Methods include:

  \describe{

    \item{sread}{\code{signature(object = "AlignedRead")}: access the
      sread slot of \code{object}.}

    \item{id}{\code{signature(object = "AlignedRead")}: access the
      id slot of \code{object}.}

    \item{[}{\code{signature(x = "SffReads", i = "ANY", j = "missing")}: 
      This method creates a new \code{SffReads} object containing only
      those reads indexed by \code{i}. Additional methods on
      \sQuote{[,SffReads} do not provide additional functionality, but
      are present to limit inappropriate use.}

    \item{append}{\code{signature(x = "SffReads", values = "SffReads", length = "missing")}:
      append the \code{sread} and \code{id} slots of \code{values} after
      the corresponding fields of \code{x}.}

    \item{narrow}{\code{signature(x = "SffReads", start = NA, end = NA, width = NA, use.names = TRUE)}:
      \sQuote{narrow} \code{sread} so that sequences are between
      \code{start} and \code{end} bases, according to
      \code{\link[IRanges:Ranges-utils]{narrow}} in the \code{IRanges}
      package.
    }

    \item{length}{\code{signature(x = "SffReads")}: returns a
      \code{integer(1)} vector describing the number of reads in this
      object.}

    \item{width}{\code{signature(x = "SffReads")}: returns an
      \code{integer()} vector of the widths of each read in this
      object.}

    \item{srorder}{\code{signature(x = "SffReads")}:}
    \item{srrank}{\code{signature(x = "SffReads")}:}
    \item{srsort}{\code{signature(x = "SffReads")}:}
    \item{srduplicated}{\code{signature(x = "SffReads")}:
      Order, rank, sort, and find duplicates in \code{SffReads} objects
      based on \code{sread(x)}, analogous to the corresponding functions
      \code{order}, \code{rank}, \code{sort}, and \code{duplicated},
      ordering nucleotides in the order \code{ACGT}.}

    \item{srdistance}{\code{signature(pattern="SffReads", subject="ANY")}:
      Find the edit distance between each read in \code{pattern} and the
      (short) sequences in \code{subject}. See \code{\link{srdistance}}
      for allowable values for \code{subject}, and for additional
      details.}

    \item{trimLRPatterns}{\code{signature(Lpattern = "", Rpattern = "", subject = "SffReads", max.Lmismatch = 0, max.Rmismatch = 0, with.Lindels = FALSE, with.Rindels = FALSE, Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)}:

      Remove left and / or right flanking patterns from
      \code{sread(subject)}, as described in
      \code{\link[Biostrings:trimLRPatterns]{trimLRPatterns}}. Classes
      derived from \code{SffReads} (e.g., \code{\link{SffReadsQ}},
      \code{\link{AlignedRead}}) have corresponding base quality scores
      trimmed, too. The class of the return object is the same as the
      class of \code{subject}, except when \code{ranges=TRUE} when the
      return value is the ranges to use to trim 'subject'.}

    \item{alphabetByCycle}{\code{signature(stringSet = "SffReads")}:
      Apply \code{\link{alphabetByCycle}} to the \code{sread} component
      of \code{stringSet}, returning a matrix as described in
      \code{\link{alphabetByCycle}}.}

    \item{tables}{\code{signature(x= "SffReads", n = 50)}:
      Apply \code{\link{tables}} to the \code{sread} component of
      \code{x}, returning a list summarizing frequency of reads in
      \code{x}.}

    \item{clean}{\code{signature(object="SffReads")}: Remove all reads
      containing non-nucleotide (\code{"N", "-"}) symbols.}

    \item{show}{\code{signature(object = "SffReads")}: provides a brief
      summary of the object, including its class, length and width.}

    \item{detail}{\code{signature(x = "SffReads")}: provides a
      more extensive summary of this object, displaying the first and
      last entries of \code{sread} and \code{id}.}

    \item{writeFasta}{\code{signature(object, file, ...)}: write
      \code{object} to \code{file} in fasta format. See
      \code{\link{writeFASTA}} for \code{...} argument values.}

  }
}
\author{Matt Settles with lots of code 'borrowed' from Martin Morgan's ShortRead Package}
\seealso{
  \code{\linkS4class{SffReadsQ}}
}
\examples{
showClass("SffReads")
showMethods(class="SffReads", where=getNamespace("rSFFreader"),
            inherit=FALSE)

showMethods(class="SffReadsQ", where=getNamespace("rSFFreader"),
            inherit=FALSE)
            
reads <- readsff(system.file("extdata","Small454Test.sff",package="rSFFreader"), use.qualities=FALSE)
detail(reads)

}
\keyword{classes}
