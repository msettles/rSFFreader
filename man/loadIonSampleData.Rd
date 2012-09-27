\name{loadIonSampleData}
\alias{loadIonSampleData}

\title{Load the example Ion Torrent dataset}

\description{
Load a small, 1000 record example sff file.
}

\usage{
load454SampleData()
}

\arguments{
None
}

\value{Loads a small example Life Sciences Ion Torrent dataset (sff file) into a \code{\linkS4class{SffReadsQ}} object containing 1000 reads.}

\author{Matt Settles <msettles@uidaho.edu>}

\examples{
##Read in the included Ion Torrent sample data set.
sff <- loadIonSampleData()
##Generate summary statistics for read lengths.
summary(width(sff))

##Generate some QA plots:
##Read length histograms:
par(mfrow=c(2,2))
clipMode(sff) <- "raw"
hist(width(sff),breaks=500,col="grey",xlab="Read Length",main="Raw Read Length")

clipMode(sff) <- "full"
hist(width(sff),breaks=500,col="grey",xlab="Read Length",main="Clipped Read Length")

## Base by position plots:
clipMode(sff) <- "raw"
ac <- alphabetByCycle(sread(sff),alphabet=c("A","C","T","G","N"))
ac.reads <- apply(ac,2,sum)
acf <- sweep(ac,MARGIN=2,FUN="/",STATS=apply(ac,2,sum))

matplot(cbind(t(acf),ac.reads/ac.reads[1]),col=c("green","blue","black","red","darkgrey","purple"),
          type="l",lty=1,xlab="Base Position",ylab="Base Frequency",
          main="Base by position")
cols <- c("green","blue","black","red","darkgrey","purple")
leg <-  c("A","C","T","G","N","\%reads")
legend("topright", col=cols, legend=leg, pch=18, cex=.8)

clipMode(sff) <- "full"
ac <- alphabetByCycle(sread(sff),alphabet=c("A","C","T","G","N"))
ac.reads <- apply(ac,2,sum)
acf <- sweep(ac,MARGIN=2,FUN="/",STATS=apply(ac,2,sum))

matplot(cbind(t(acf),ac.reads/ac.reads[1]),col=c("green","blue","black","red","darkgrey","purple"),
          type="l",lty=1,xlab="Base Position",ylab="Base Frequency",
          main="Base by position")
legend("topright", col=cols, legend=leg, pch=18, cex=.8)
}