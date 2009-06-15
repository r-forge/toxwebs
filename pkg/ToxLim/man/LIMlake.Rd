\name{LIMlake}
\alias{LIMlake}
\alias{LIMlakeFish}
\docType{data}
\title{
  Linear Inverse Model input for two
  Experimental freshwater lake enclosure (Lac Croche, Canada)
}

\description{
  Linear inverse model specification for two freshwater lake enclosures
  (Lac Croche, Canada), used as demonstration example in Laender et al., 2009.

  It comprises the enclosure without \code{LIMlake} and with \code{LIMlakeFish}
  fish.

  Hier wat meer uitleg  over dit experiment...
}

\usage{
LIMlake
LIMlakeFish
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
pm <- par(mfrow=c(1,2))
lake <- Flowmatrix(LIMlake)
plotweb(lake,main="Lake enclosure food web",sub="mgC/m2/day")

plotweb(Flowmatrix(LIMlakeFish),main="Lake enclosure food web + FISH",
  sub="mgC/m2/day")

# ranges of flows
Plotranges(LIMlake,lab.cex=0.7,xlab="mgC/m2/d",
   main="Flows")
# ranges of variables
Plotranges(LIMlake,type="V",lab.cex=0.7,xlab="mgC/m2/d",
   main="variables")
mtext(outer=TRUE,"Lake enclosure food web",side=3,line=-1,cex=1.5)
par(mfrow=pm)

}

\references{
  De Laender, F., Van Oevelen, D., Middelburg, J.J. and Soetaert, K., 2009.
  Incorporating Ecological Data and Associated Uncertainty in
  Bioaccumulation Modeling: Methodology Development and Case Study.
  Environ. Sci. Technol., 2009, 43 (7), 2620-2626.

  Ridal, J. J., Mazumder, A., Lean, D.R.S., 2001. Effects of nutrient loading
  and planktivory on the accumulation of organochlorine pesticides
  in aquatic food chains. Environ. Toxicol. Chem. 20, 1312 - 1319.
}

\seealso{
  \code{\link{LimOmega}}
}

\keyword{datasets}

