\name{LIMbarents}
\alias{LIMbarents}
\docType{data}
\title{
  Linear Inverse Model Input For the Barents Sea
}

\description{
  Linear inverse model specification for the Southern Barents Sea

  Food web compartments for the area are dissolved organic carbon (DOC), 
  detritus, bacteria, heterotrophic flagellates and ciliates, 
  phytoplankton (pico- and nanoplankton, diatoms and Phaeocystis sp.), 
  mesozooplankton (copepods), macrozooplankton (krill and chaetognaths), 
  cod (Gadus morhua), herring (Clupea harengus) and capelin (Mallotus villosus).
  Adult cod (> 3 yrs) and young cod (< 3 yrs) were considered as two different populations.  
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
  Frederik de Laender <frederik.delaender@ugent.be>

  Karline Soetaert <k.soetaert@nioz.nl>

}

\examples{

plotweb(Flowmatrix(LIMbarents),main="Barents Sea Food Web",
  sub="gC/m2/day")

pm <- par(mfrow=c(1,2))
# ranges of flows
Plotranges(LIMbarents,lab.cex=0.7,xlab="gC/m2/d",
   main="Flows")
# ranges of variables
Plotranges(LIMbarents,type="V",lab.cex=0.7,xlab="gC/m2/d",
   main="variables")
mtext(outer=TRUE,"Southern Barents Sea food web",side=3,line=-1,cex=1.5)
par(mfrow=pm)

}

\references{

  De Laender, F., Van Oevelen, D., Middelburg, J.J. and Soetaert, K., 2010. 
  Carbon transfer in herbivore- and microbial loop-dominated pelagic 
  food webs in the southern Barents Sea during spring and summer. 
  Marine Ecology Progress Series 398: 93 -- 107.

  De Laender, F., Van Oevelen, D., Frantzen, S., Middelburg, J.J. and Soetaert, K. 
  Seasonal PCB Bioaccumulation in an Arctic Marine Ecosystem: A
  Model Analysis Incorporating Lipid Dynamics, Food-Web Productivity
  and Migration. Environ. Sci. Technol. 2010, 44, 356 -- 361.
}

\note{
An application of \code{LIMbarents} with \link{LimOmega} can be found in 
De Laender et al. (2010).
}
\seealso{
  \code{\link{LIMlake}},  \code{\link{LIMlakeFish}}
}

\keyword{datasets}

