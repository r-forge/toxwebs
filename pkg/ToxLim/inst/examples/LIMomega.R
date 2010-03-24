require(ToxLim)

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
