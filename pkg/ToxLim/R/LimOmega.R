
LimOmega <- function(lim=NULL,
                     flowmatrix=NULL,
                     INST, KINE, EXPO,
                     DOC=NULL,  DET,  DIC,
                     WW_KINE=7.596535e-08,
                     SS_KINE=0.079508,
                     Growth=0, k=0.25, Q=1,
                     logKow=6,  Cwater=10,  Koc_Kow=0.41,
                     LIPID_INST=0.0053, LIPID_KINE=0.02,
                     OC =0.028,  rH2O = 1.1e-5,  rH2O_0 = 2.8e-3,
                     rCH2 = 68,  g_0 = 200   )
{

## -----------------------------------------------------------------------------
## Initialisation: check name consistency
## -----------------------------------------------------------------------------
  if (is.null(lim) & is.null(flowmatrix))
    stop("Either lim or flowmatrix should be specified")

# Extract flowmatrix
  if (is.null(flowmatrix)) {
    flowmatrix         <-Flowmatrix(lim)
  }

  LimNull <- is.null(lim)
  limNames <- colnames(flowmatrix)

  # Function to check if an input name is known to the LIM
  checkname <- function(Name,string) {
    ii <- which (! Name  %in% limNames)
    if (length(ii) >0) stop(paste(string,"not known from lim: ",Name[ii]))
  }

  checkname(INST,"INST")
  checkname(EXPO,"EXPO")
  if (! is.null(DOC))
    checkname(DOC,"DOC")
  checkname(DET,"DET")
  checkname(DIC,"DIC")

  #check if number of wet weights and SS equals length of KINE
  if (length(WW_KINE) == 1)
    WW_KINE <- rep(WW_KINE,len=length(KINE))

  if (length(SS_KINE) == 1)
    SS_KINE <- rep(SS_KINE,len=length(KINE))

  if (length(WW_KINE) != length(KINE))
    stop ("check WW_KINE; length should be = 1 or = length of KINE")

  if (length(SS_KINE) != length(KINE))
    stop ("check SS_KINE; length should be = 1 or = length of KINE")

  nCOMPS<- length(KINE)+length(INST)
## -----------------------------------------------------------------------------
## Calculate chemical specifications of experiment
## -----------------------------------------------------------------------------
  Kow <- 10^logKow
  Koc <- Koc_Kow * Kow

  LIPID_INST <- rep(LIPID_INST,len=length(INST))
  LIPID_KINE <- rep(LIPID_KINE,len=length(KINE))

  LIPID <- c(LIPID_INST,LIPID_KINE)
  OC    <- rep(OC,len=length(INST))

## -----------------------------------------------------------------------------
## producing matrices that are needed in LIMOMEGA
## -----------------------------------------------------------------------------

  # diagonal matrix with LIPID fractions
  LIPIDm <- matrix(data=0,ncol=nCOMPS,nrow=nCOMPS)
  diag(LIPIDm) <- LIPID
  colnames(LIPIDm) <- c(INST,KINE)

  # diagonal matrix with 1/(LIPID fractions_diet*(Kow-1)+1) used in "UP_B" calc
  LIPIDmext <- matrix(data=0,ncol=nCOMPS,nrow=nCOMPS)
  diag(LIPIDmext) <- 1/(LIPID*(Kow-1)+1)
  colnames(LIPIDmext) <- c(INST,KINE)

  # diagonal matrix with 1/(LIPID fractions_pred*(Kow-1)+1) used in "OUT_A" calc
  LIPIDmexxt <- matrix(data=0,ncol=length(KINE),nrow=length(KINE))
  diag(LIPIDmexxt) <- 1/(LIPID[which(colnames(LIPIDm)%in%KINE)]*(Kow-1)+1)
  colnames(LIPIDmexxt) <- KINE
  # diagonal matrix with WW^(-k) used in "G" and "UP_B" calc
  WW_KINEm <- matrix(data=0,ncol=length(KINE),nrow=length(KINE))
  diag(WW_KINEm) <- WW_KINE^(-k)
  colnames(WW_KINEm) <- KINE
  #inversion of WW_KINEm
  WW_KINEm_inv <- (WW_KINEm)^(-1)
  WW_KINEm_inv[is.infinite(WW_KINEm_inv)] <- 0

  # diagonal matrix with reciprokals of SS used in "G" and "k_dil" calc
  SS_KINEm <- matrix(data=0,ncol=length(KINE),nrow=length(KINE))
  diag(SS_KINEm) <- 1/SS_KINE
  colnames(SS_KINEm) <- KINE

## -----------------------------------------------------------------------------
## Predict the Bioaccumulation factors for the
## KINE populations (µg kg ww-1 / µg L-1) = internal concentrations divided
## by concentration in water
## -----------------------------------------------------------------------------

## re-arrange flow matrix: INST KINE EXPO
  Fl               <- flowmatrix[c(INST,KINE,EXPO),c(INST,KINE,EXPO)]

## take submatrix from Fl for KINE calculations
  Fl_KIN           <- matrix(data=Fl[,which(colnames(Fl) %in% KINE)],
                               ncol=length(KINE))
  colnames(Fl_KIN) <- KINE
  rownames(Fl_KIN) <- rownames(Fl)

## remove externals
  Fl_KINE          <- matrix(data=Fl_KIN[-(which(rownames(Fl_KIN) %in% EXPO)),],
                             ncol=length(KINE))
  colnames(Fl_KINE) <- KINE
  rownames(Fl_KINE) <- c(INST,KINE)

## -----------------------------------------------------------------------------
## calculation of chemical in those food web compartments that are assumed to
## be in rapid equilibrium with water (INST)
## -----------------------------------------------------------------------------

  X_INST             <- matrix(data=(Cwater * Koc * OC), ncol=1,
                               nrow=length(INST))

## -----------------------------------------------------------------------------
## calculation of chemical in those food web compartments for which
## uptake/loss kinetics are explicitly modelled (KINE)
## -----------------------------------------------------------------------------

  #food uptake
  Uptake             <- c(colSums(Fl_KINE))

  #diagonal matrix of uptakes by KINEs
  Uptake_KINEm       <- matrix(data=0,ncol=length(KINE),nrow=length(KINE))
  diag(Uptake_KINEm) <- 1/Uptake

  #diet matrix: proportions of different food sources
  Diet               <- Fl_KINE %*% Uptake_KINEm
  colnames(Diet)     <- KINE

  #egestion
  DETLoss            <- c(Fl[which(rownames(Fl) %in% KINE),
                             which(colnames(Fl) == DET)])
  DETLossm           <- matrix(data=DETLoss,ncol=length(KINE),nrow=1)

  #assimilation efficiency
  p                  <- (Uptake - DETLossm) / Uptake

  # 1-assimilation efficiency for KINEs in a diagonal matrix
  pm             <- matrix(data=0,ncol=length(KINE),nrow=length(KINE))
  colnames(pm)   <- KINE
  diag(pm)       <- (1-p)

  #diagonal matrix of assimilation efficiency /(1-assimilation efficiency)
  palt           <- matrix((data=0),ncol=length(KINE),nrow=length(KINE))
  colnames(palt) <- KINE
  diag(palt)     <- p/(1-p)

## -----------------------------------------------------------------------------
## BIOACCUMULATION
## -----------------------------------------------------------------------------

  #chemical uptake from water
  k_uptake_w     <- matrix(data=(WW_KINE^(-k) / (rH2O_0 + rCH2/Kow)),
                           nrow=length(KINE),ncol=1)

  #chemical excretion to water - rate constant calculation
  k_excretion    <- 1/(LIPID[which(colnames(Fl) %in% KINE)]*(Kow-1)+1) *
                        WW_KINE^(-k)/(rH2O_0+rCH2/Kow)

  #chemical uptake from food - ingestion coefficient
  G              <- (Fl_KINE %*% SS_KINEm) %*% WW_KINEm_inv

  #chemical dilution over body: a result of growth/production
  k_dil          <- matrix(data=Fl[which(rownames(Fl) %in% KINE),
                                   which(colnames(Fl) %in% c(KINE,EXPO))],
                           nrow=length(KINE),ncol=length(KINE)+length(EXPO))
  k_dil          <- rowSums(k_dil)
  k_dil          <- c(SS_KINEm %*% k_dil)

  #chemical dilution over body: a result of net growth (rate of change: is zero in steady state LIM)
  k_grdil        <- c(Growth/SS_KINE)

  #chemical uptake from food - rate constant calculation
  MC1            <- matrix(data=c(rH2O + rCH2/(Kow*Q)),
                           nrow=nCOMPS,ncol=length(KINE),byrow=TRUE)
  UP_B           <- (1/(MC1 + 1/(((LIPIDm %*% G) %*% pm) * Kow * Q)))%*%WW_KINEm
  UP_A           <- LIPIDmext
  k_uptake       <- UP_A %*% (UP_B %*% palt)  #d-1

  rownames(k_uptake)  <- rownames(Fl_KINE)
  colnames(k_uptake)  <- KINE
  k_uptake            <- t(k_uptake)
  rowK <-

  #uptake through feeding on KINE groups
  k_uptake_KINE   <- matrix(data=(k_uptake[which(rownames(k_uptake) %in% KINE),
                                           which(colnames(k_uptake) %in% KINE)]),
                            nrow=length(KINE),ncol=length(KINE))
  #uptake through feeding on INST groups
  k_uptake_INST   <- matrix(data=(k_uptake[which(rownames(k_uptake) %in% KINE),
                                           which(colnames(k_uptake) %in% INST)]),
                            nrow=length(KINE),ncol=length(INST))

  #chemical egestion to faeces - rate constant calculation
  OUT_A               <- LIPIDmexxt
  OUT_B               <- matrix(data=(colSums(UP_B)),nrow=1,ncol=length(KINE))
  k_egestion          <- OUT_B %*% OUT_A
  colnames(k_egestion)<- KINE
  k_egestion          <- c(colSums(k_egestion))

  #all loss processes of chemical - rate constant calculation
  k_loss              <- matrix(data=0,nrow=length(KINE),ncol=length(KINE))
  diag(k_loss)        <- k_egestion + k_dil + k_excretion +k_grdil

  #net rate of accumulation through eating KINE compartments and loosing
  # chemical by egestion, dilution, excretion and growth dilution - rate calculation
  k_net               <- k_uptake_KINE- k_loss

  #chemical uptake from eating INST compartments
  INS   <- k_uptake_INST %*% X_INST

  #chemical uptake from water
  W     <- k_uptake_w * Cwater

  #constant term "B" in subsequent equation
  Ct    <- INS + W

## -----------------------------------------------------------------------------
##solution of Xdot = AX + B with
## X the vector of internal concentrations and Xdot equal to zero
## -----------------------------------------------------------------------------

  BCF                          <- X_INST / Cwater                                                      # Bioconcentration factor (µg kg ww-1 / µg L-1) for INST groups = internal concentartions divided by concentration in water

  BAF                          <- c(solve(k_net,-Ct)) / Cwater                                         # Bioaccumulation factor (µg kg ww-1 / µg L-1) for KINE groups = internal concentartions divided by concentration in water

  BCF                          <- t(BCF)

  BAF                          <- t(BAF)

## -----------------------------------------------------------------------------
## Save results
## -----------------------------------------------------------------------------

## store predicted BAF after lipid normalisation: (µg kg lipid-1 / µg L-1)
  BAF_LC          <- as.vector(BAF/LIPID_KINE)
  names(BAF_LC)   <- KINE

## store predicted BCF after organic carbon normalisation: (µg kg OC-1 / µg L-1)
  BCF_OC          <- as.vector(BCF/OC)
  names(BCF_OC)   <- INST

  # return as a list
  return(list (BAF_LC=BAF_LC, BCF_OC=BCF_OC))
}
