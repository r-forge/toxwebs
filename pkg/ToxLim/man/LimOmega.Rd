\name{LimOmega}
\alias{LimOmega}
\title{
  LIM OMEGA bioaccumulation model
}
\description{
  Solves a LIMOMEGA model

  Uncertainty ranges of specific ingestion and production rates, and of assimilation efficiency are estimated with a LIM, 
  and used to set Kup,food and Kout+dil in OMEGA, as explained in detail in De Laender et al (2009). 
  If the LIM is solved using the Markow Chain Monte Carlo approach developed by Van den Meersche et al (2009), 
  uncertainty ranges for Kup,food and Kout+dil and thus for the internal concentrations can be estimated.
}

\usage{
LimOmega (lim=NULL, flowmatrix=NULL,
  INST, KINE, EXPO, DOC=NULL,
  DET, DIC, WW_KINE=7.596535e-08,
  SS_KINE=0.079508, Growth=0, k=0.25, Q=1,
  logKow=6,  Cwater=10,  Koc_Kow=0.41,
  LIPID_INST=0.0053,  LIPID_KINE=0.02,
  OC =0.028,  rH2O = 1.1e-5,  rH2O_0 = 2.8e-3,
  rCH2 = 68,  g_0 = 200  )
}

\arguments{
  \item{lim }{a list that contains the linear inverse model
    specification, when \code{NULL} then \code{flowmatrix} should be specified.
  }
  \item{flowmatrix }{a flow matrix, with flows *from* as rows, flows
    *to* as columns. The matrix should be square; its column and row names
    should contain the names of the food web components; they should be equal.
    If \code{NULL}, then created from the \code{lim}
  }
  \item{INST }{the names of those food web compartments in \code{lim} that are
    assumed to be in rapid equilibrium with water.
  }
  \item{KINE }{the names of those food web compartments in \code{lim} for which
    uptake/loss kinetics are explicitly modelled.
  }
  \item{EXPO }{the names of all \code{lim} externals that do not represent
    dissolved inorganic carbon.
  }
  \item{DOC }{the name of the \code{lim} food web compartment that represents
    dissolved organic carbon; the default is no such compartment.
  }
  \item{DET }{the name of the \code{lim} food web compartment that represents
    detritus.
  }
  \item{DIC }{the name of \code{lim} food web compartment that represents
    dissolved inorganic carbon.
  }
  \item{WW_KINE }{wet weights of \code{KINE} individuals, in same order as in
    \code{KINE}: [kg wet weight].
  }
  \item{SS_KINE }{standing stocks of \code{KINE} individuals, in same order as
    in \code{KINE}: [g C/m2].
  }
  \item{Growth }{growth rates of \code{KINE} individuals, essentially the
    rate of change of the KINE food web compartments: [/d].
  }
  \item{k }{shape factor; relates rates of uptake/loss of carbon to
    uptake/loss of chemical: [-].
  }
  \item{Q }{correction factor; default (1) is for cold blooded,
    use 10 for warm blooded animals: [-].
  }
  \item{logKow }{10-base logarithm of octanol water partition coefficient:
    [log10 (L/kg)].
  }
  \item{Cwater }{dissolved concentration: [microg/L].
  }
  \item{Koc_Kow }{ratio of Koc (organic carbon water partitioning coefficient)
    over Kow
  }
  \item{LIPID_INST }{proportion of wet weight that is lipid for the \code{INST}
    food web compartments, either one value or one value per \code{INST}: [-].
  }
  \item{LIPID_KINE }{proportion of wet weight that is lipid for the \code{KINE}
    food web compartments, either one value or one value per \code{KINE}: [-].
  }
  \item{OC }{proportion of wet weight that is organic carbon in \code{INST}
    groups: [-].
  }
  \item{rH2O }{water layer diffusion resistance for uptake (loss) of chemical
    from food (through egestion): [\eqn{kg^{-k}d }{d kg^(-k)}].
  }
  \item{rH2O_0 }{water layer diffusion resistance for absorption (excretion)
    of chemical from (to) water: [\eqn{kg^{-k}d}{d kg^(-k)}].
  }
  \item{rCH2 }{lipid layer permeation resistance: [\eqn{kg^{-k}d}{d/kg^k}].
  }
  \item{g_0 }{water absorption - excretion coefficient: [\eqn{kg^k d^{-1}}{kg^(k)/d}].
  }
}
\value{
  a list containing:

  \item{BAF_LC }{The bioaccumulation factor, after lipid normalisation.
    units [mg (kg lipid)-1 / mg L-1]
  }
  \item{BCF_OC }{The bioconcentration factor predicted BCF after organic
    carbon normalisation. units [mg (kg OC)-1 / mg L-1]
    UITLEG OVER WAAROM ALLEMAAL GELIJK
  }
}

\author{
  Frederik de Laender <f.delaender@nioo.knaw.nl>

  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\examples{

#------------------------------------------------------------------
# Three simple examples: the bioaccumulation model for the
# three food webs, estimated once...
#------------------------------------------------------------------

LimOmega (lim = LIMlake,
          INST=c("DET","DOC","BAC","PHY","NAN","CIL","MIZ"),
          KINE= "MEZ", EXPO=c("sed","gr"),
          DOC="DOC", DET="DET", DIC="dic",
          WW_KINE=7.596535e-08, SS_KINE=0.079508)

LimOmega (lim = LIMlakeFish,
          INST=c("DET","DOC","BAC","PHY","NAN","CIL","MIZ"),
          KINE=c("MEZ","FIS"), EXPO=c("sed","gr"),
          DOC="DOC", DET="DET", DIC="dic",
          WW_KINE=c(1.85e-08,1.7e-3),
          SS_KINE=c(0.00514800,0.26))

LimOmega (lim = LIMbarents,
          INST=c("DIA","PHA","AUT","CIL","HNA","DET","BAC"),
          KINE=c("COP","CHA","KRI","CAP","COD","YCO","HER"),
          EXPO=c("SED","GRA","GRO"),
          DOC="DOC", DET="DET", DIC="DIC",
          WW_KINE=c(0.000001,8e-05,0.00006,10e-3,3,1,20e-3),
          LIPID_KINE=c(0.01,0.01,0.01,0.03,0.03,0.03,0.03),
          LIPID_INST=0.04,
          SS_KINE=c(1.79,0.6965,0.003,0.38,0.053,0.006,0.055))

#------------------------------------------------------------------
# Now performing a monte carlo run on food web structure
#------------------------------------------------------------------
# 1. Take niter random samples from all possible solutions using a
#     Markow Chain Monte Carlo approach
  X0       <- Lsei(LIMlake)$X
  niter <- 50
  SolXS    <- Xsample(LIMlake, iter=niter, type = "mirror",
                      jmp=NULL, x0=X0, fulloutput = FALSE)

  BAFlc_all  <- NULL

# 2. For each of these samples: create flowmatrix and run LimOmega
  for (i in 1:niter) {
   flowmat <- Flowmatrix(LIMlake, SolXS[i,])
   LO<- LimOmega (flowmatrix=flowmat,
          INST=c("DET","DOC","BAC","PHY","NAN","CIL","MIZ"),
          KINE=c("MEZ"), EXPO=c("sed","gr"),
          DOC="DOC", DET="DET", DIC="dic",
          WW_KINE=7.596535e-08, SS_KINE=0.079508,
          Growth=0, k=0.25, Q=1)

   BAFlc_all <- c(BAFlc_all,LO$BAF_LC)
  }
# 3. show results
hist(BAFlc_all)

#------------------------------------------------------------------
# Same food web structure, different chemical parameters
#------------------------------------------------------------------
  niter <- 100
  
  # a normally distributed sample of log kow
  lkw <- rnorm(niter,mean=6,sd=0.4)

  BAFlc_all  <- NULL
  BCFoc_all  <- NULL
  
  # the flowmatrix on which this is based
  flowmat <- Flowmatrix(LIMlake)

  for (i in 1:niter) {
   LO<- LimOmega (flowmatrix=flowmat,
          INST=c("DET","DOC","BAC","PHY","NAN","CIL","MIZ"),
          KINE=c("MEZ"), EXPO=c("sed","gr"),
          DOC="DOC", DET="DET", DIC="dic",
          WW_KINE=7.596535e-08, SS_KINE=0.079508,
          Growth=0, k=0.25, Q=1, logKow = lkw[i])

   BAFlc_all <- c(BAFlc_all,LO$BAF_LC)
   BCFoc_all <- c(BCFoc_all,LO$BCF_OC[1])
  }

pm <- par(mfrow=c(2,2))
hist(BAFlc_all,main="BAF_LC")
plot(lkw,BAFlc_all,xlab="log Kow",ylab="BAF_LC")
hist(BCFoc_all,main="BCF_OC")
plot(lkw,BCFoc_all,xlab="log Kow",ylab="BCF_OC")
par(mfrow=pm)

}

\references{
  De Laender, F., Van Oevelen, D., Middelburg, J.J. and Soetaert, K., 2009.
  Incorporating Ecological Data and Associated Uncertainty in
  Bioaccumulation Modeling: Methodology Development and Case Study.
  Environ. Sci. Technol., 2009, 43 (7), 2620-2626.

  Hendriks, A.J., van der Linde, A., Cornelissen, G., Sijm, D., 2001. The
  power of size. 1. Rate constants and equilibrium ratios for
  accumulation of organic substances related to octanol-water
  partition ratio and species weight. Environ. Toxicol. Chem. 20, 1399 - 1420.
  
  Van den Meersche, K., Soetaert, K., Van Oevelen, D., 2009. 
  xsample(): An R Function for Sampling Linear Inverse Problems.
  J. Stat. Soft. 30, 1-15. 
}

\details{
  Solves the following inverse problem:

  DIT HIER VERVANGEN DOOR HET OMEGA MODEL. LaTeX style maar geen echte latex
  \deqn{\min(||AAx-BB||^2)}, the approximate equations
    subject to
  \deqn{Ex=F}, the mass balances
  \deqn{Gx>=H}, the constraints.



}
\seealso{
  \code{LIMlake},  \code{LIMlakeFish},  \code{LIMbarents}, the input food webs.
}

\keyword{algebra}

