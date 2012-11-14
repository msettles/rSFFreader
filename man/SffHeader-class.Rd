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

\title{\code{SffHeader}}

\description{
Class \code{SffHeader} contains meta-data stored in the header of the SFF files read in.

Objects from this class are the result of \code{\link{readSffHeader}},
or from the result of \code{\link{readSff}}. The resulting object will
contain a \code{header} slot which is a list. If multiple sff files 
were processed by either of the above functions, this list will contain
meta-data about each of these files in corresponding positions in the list.

Meta-data included is defined in the SFF file specifications and include:
\describe{
    \item{filename:}{The name of the file that was read in.}
    \item{magic_number:}{779314790, which encodes the string ".sff"}
    \item{version:}{Version number}
    \item{index_offset:}{An optional field which indicates the position of a read index within the file.}
    \item{index_length:}{An optional field which indicates the length of a read index within the file.}
    \item{number_of_reads:}{Stores the number of reads in the file.}
    \item{header_length:}{The number of bytes required by header fields.}
    \item{key_length:}{The length of the key sequence used for these reads.}
    \item{number_of_flows_per_read:}{The number of flows carried out during sequencing of the reads.}
    \item{flowgram_format_code:}{Indicates the format of the flowgram encoding. Currently "1" is the only valid value.}
    \item{flow_chars:}{Indicates the nucleotide bases used for each sequencing flow.}
    \item{key_sequence:}{Nucleotide sequence use for these reads.}
  }
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

## The readSffHeader can be used to extract header information from one or more sff files:
sffFiles = c(system.file("extdata", "SmallTorrentTest.sff", package = "rSFFreader"), 
             system.file("extdata","Small454Test.sff",package="rSFFreader"))
header <- readSffHeader(sffFiles)
header
header(header)[[1]]$number_of_reads
header(header)[[2]]$number_of_reads

## Header information is also retrieved when using readSff:
sff <- readSff(sffFiles)

## Number of flows for dataset 1 and 2:
header(sff)[[1]]$number_of_flows_per_read
header(sff)[[2]]$number_of_flows_per_read

}
\keyword{classes}
