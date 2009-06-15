\name{LIMbarents}
\alias{LIMbarents}
\docType{data}
\title{
  Linear Inverse Model input for the Barents Sea
}

\description{
  Linear inverse model specification for the Barents Sea

  Hier wat meer uitleg  over dit experiment, maar niet te veel zodat
  "ze" er niet mee weglopen...
}

\usage{
LIMbarents
}

\format{
  a list of matrices, vectors, names and values that specify the
  linear inverse model problem.

  see the return value of \code{\link[LIM]{Setup}} for more information
  about this list

  A more complete description of this structures is in vignette("LIM") from
  package \code{LIM}
}

\author{
  Frederik de Laender <f.delaender@nioo.knaw.nl>

  Karline Soetaert <k.soetaert@nioo.knaw.nl>

}

\examples{

plotweb(Flowmatrix(LIMbarents),main="Barents Sea Food Web",
  sub="mgC/m2/day")

pm <- par(mfrow=c(1,2))
# ranges of flows
Plotranges(LIMbarents,lab.cex=0.7,xlab="mgC/m2/d",
   main="Flows")
# ranges of variables
Plotranges(LIMbarents,type="V",lab.cex=0.7,xlab="mgC/m2/d",
   main="variables")
mtext(outer=TRUE,"Barents Sea food web",side=3,line=-1,cex=1.5)
par(mfrow=pm)

}

\references{

  De Laender, F., Van Oevelen, D., Middelburg, J.J. and Soetaert, K., 2009.
  MEPS VERHAAL

}

\seealso{
  \code{\link{LIMlake}},  \code{\link{LIMlakeFish}}
}

\keyword{datasets}

