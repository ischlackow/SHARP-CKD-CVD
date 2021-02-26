###################################################################
#                       libraries              	                  #
###################################################################

library(data.table)
library(plyr)  
library(scales)
library(ggplot2)
#library(reshape2)
library(xtable)
library(foreach)
library(snowfall)
library(doSNOW)

###################################################################
#          functions from "extrapolation, preliminary" file       #
###################################################################

get_state_combined <- function(N_CV, N_CKD) {
  return(paste(N_CV, N_CKD, sep = "_"))
}

### possible transitions ###

.get_next_N_CV <- function(N_CV) {
  if (N_CV == 1)
    Ns_CV1 <- c(1, 2, 5) else
      if (N_CV %in% c(2, 3, 4))
        Ns_CV1 <- c(2, min(N_CV + 1, 4)) else
          if (N_CV == 5)
            Ns_CV1 <- c(2, 5, 6) else
              if (N_CV %in% c(6, 7))
                Ns_CV1 <- c(2, 5, min(N_CV + 1, 7))
  return(Ns_CV1)
}

.get_next_N_CKD <- function(N_CKD) {
  if (N_CKD %in% c(1, 2, 3))
    Ns_CKD1 <- c(1, 2, 3, 4, 8) else {
      tr <- N_CKD < 8
      Ns_CKD1 <- c(min(N_CKD + 1 - 4 * (1 - tr), 7), min(N_CKD + 1 + 4 * tr, 11))
    }
  return(Ns_CKD1)
}

.get_rs <- function(N_CKD0, N_CKD1) {
  # pre-ESRD
  if (N_CKD1 %in% c(1, 2, 3))
    rs <- N_CKD1 else
      if (N_CKD1 %in% c(4, 5, 6, 7)) { # transplant
        if (!(N_CKD0 %in% c(4, 5, 6, 7)))
          rs <- 4 else
            rs <- 5
        } else 
          # dialysis
          if (N_CKD1 %in% c(8, 9, 10, 11)) { # dialysis
            if (!(N_CKD0 %in% c(8, 9, 10, 11)))
              rs <- 6 else
                rs <- 7
          }
  return(rs)
}

### combined ###

.get_stateN <- function(states_info, BVD, N_CV, N_CKD) {
  if (N_CV == 1) {
    DH_CV <- if (BVD) "no NFMVE, baseline VD" else "no NFMVE, no baseline VD"  
    N <- which(sapply(states_info, "[", "N_CV") == N_CV 
               & sapply(states_info, "[", "N_CKD") == N_CKD
               & sapply(states_info, "[", "DH_CV") == DH_CV)
  } else
    N <- which(sapply(states_info, "[", "N_CV") == N_CV 
               & sapply(states_info, "[", "N_CKD") == N_CKD)
  names(N) <- NULL
  return(N)
}

# indicator of initial CKD state
.get_N_CKD_0 <- function(stageRand_2, ESRDDuration_T) {
  if (stageRand_2 == "1-3b")
    N_0 <- 1
  if (stageRand_2 == "4")
    N_0 <- 2
  if (stageRand_2 == "5")
    N_0 <- 3
  if (stageRand_2 == "transplant" & ESRDDuration_T == 0)
    N_0 <- 4
  if (stageRand_2 == "transplant" & ESRDDuration_T == 1)
    N_0 <- 5
  if (stageRand_2 == "transplant" & ESRDDuration_T == 2)
    N_0 <- 6
  if (stageRand_2 == "transplant" & ESRDDuration_T >= 3)
    N_0 <- 7
  if (stageRand_2 == "dialysis" & ESRDDuration_T == 0)
    N_0 <- 8
  if (stageRand_2 == "dialysis" & ESRDDuration_T == 1)
    N_0 <- 9
  if (stageRand_2 == "dialysis" & ESRDDuration_T == 2)
    N_0 <- 10
  if (stageRand_2 == "dialysis" & ESRDDuration_T >= 3)
    N_0 <- 11
  return(N_0)
}

# indicator of initial CV state
.get_N_CV_0 <- function(CVRand) {
  if (CVRand %in% c("None", "Another cardiovascular event"))
    N_0 <- 1
  if (CVRand == "Major atherosclerotic event in the last year")
    N_0 <- 2
  if (CVRand == "Major atherosclerotic event 1-2 years ago")
    N_0 <- 3
  if (CVRand == "Major atherosclerotic event >2 years ago")
    N_0 <- 4
  if (CVRand == "No MAE, but haemorrhagic stroke in the last year")
    N_0 <- 5
  if (CVRand == "No MAE, but haemorrhagic stroke 1-2 years ago")
    N_0 <- 6
  if (CVRand == "No MAE, but haemorrhagic stroke >2 years ago")
    N_0 <- 7
  return(N_0)
}

### possible next states

### possible next CV states
.get_Ns_CV <- function(N_0, cycle) {
  if (N_0 == 1) {
    if (cycle == 1)
      Ns <- 1
    if (cycle == 2)
      Ns <- c(1, 2, 5)
    if (cycle == 3)
      Ns <- c(1, 2, 3, 5, 6)
    if (cycle >= 4)
      Ns <- c(1, 2, 3, 4, 5, 6, 7)
  } else if (N_0 == 2) {
    if (cycle == 1)
      Ns <- 2
    if (cycle == 2)
      Ns <- c(2, 3)
    if (cycle >= 3)
      Ns <- c(2, 3, 4)
  } else if (N_0 == 3) {
    if (cycle == 1)
      Ns <- 3
    if (cycle == 2)
      Ns <- c(2, 4)
    if (cycle >= 3)
      Ns <- c(2, 3, 4)
  } else if (N_0 == 4) {
    if (cycle == 1)
      Ns <- 4
    if (cycle == 2)
      Ns <- c(2, 4)
    if (cycle >= 3)
      Ns <- c(2, 3, 4)
    } else if (N_0 == 5) {
    if (cycle == 1)
      Ns <- 5
    if (cycle == 2)
      Ns <- c(2, 5, 6)
    if (cycle == 3)
      Ns <- c(2, 3, 5, 6, 7)
    if (cycle >= 4)
      Ns <- c(2, 3, 4, 5, 6, 7)
  } else if (N_0 == 6) {
    if (cycle == 1)
      Ns <- 6
    if (cycle == 2)
      Ns <- c(2, 5, 7)
    if (cycle == 3)
      Ns <- c(2, 3, 5, 6, 7)
    if (cycle >= 4)
      Ns <- c(2, 3, 4, 5, 6, 7)
  } else if (N_0 == 7) {
    if (cycle == 1)
      Ns <- 7
    if (cycle == 2)
      Ns <- c(2, 5, 7)
    if (cycle == 3)
      Ns <- c(2, 3, 5, 6, 7)
    if (cycle >= 4)
      Ns <- c(2, 3, 4, 5, 6, 7)
  }
  
  return(Ns)
}

### possible next CKD states

.get_Ns_CKD_preESRD <- function(N_0, cycle) {
  if (cycle == 1)
    Ns <- c(N_0)
  if (cycle == 2)
    Ns <- c(1, 2, 3, 4, 8)
  if (cycle == 3)
    Ns <- c(1, 2, 3, 4, 5, 8, 9)
  if (cycle == 4)
    Ns <- c(1, 2, 3, 4, 5, 6, 8, 9, 10)
  if (cycle >= 5)
    Ns <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
  return(Ns) 
}

.get_Ns_CKD_ESRD <- function(N_0, cycle) {
  if (cycle == 1)
    Ns <- c(N_0)
  N_temp <- if (N_0 >= 8) N_0 - 4 else N_0  # transplant status corresponding to the same time as N_0
  if (cycle == 2)
    Ns <- c(min(N_temp + 1, 7), min(N_temp + 5, 11))
  if (cycle == 3)
    Ns <- c(min(N_temp + 2, 7), min(N_temp + 6, 11))
  if (cycle >= 4)
    Ns <- c(7, 11)
  return(Ns) 
}

.get_Ns_CKD <- function(N_0, cycle) {
  Ns <- if (N_0 <= 3) .get_Ns_CKD_preESRD(N_0, cycle) else .get_Ns_CKD_ESRD(N_0, cycle)
  return(Ns)
}

###################################################################
#                        preliminary calculations        		      #
###################################################################

### risk equations

prelimRiskEquations <- function(vars_0, dfBaselineMM_0,
                                NVDExternal, 
                                coeff_VD, coeff_NFMAE, coeff_NFMVE, coeff_dial_into_tr, coeff_from_nonESRD) {
  
  ### split coefficients into baseline and time-updated

  vars <- c("coeff_VD", "coeff_NFMAE", "coeff_NFMVE")
  
  for (var in vars) {
    names_0 <- intersect(vars_0, names(get(var)))
    names_T <- setdiff(names(get(var)), names_0)
    assign(paste(var), list(#coeff = get(var), 
                            coeff_0 = get(var)[names_0],
                            coeff_T = get(var)[names_T]))
  
  }

  for (var in c("coeff_dial_into_tr")){
    names_0 <- intersect(vars_0, names(get(var)))
    names_T <- setdiff(names(get(var)), names_0)
    assign(paste(var), list(#coeff = get(var), 
                            coeff_0 = get(var)[names_0],
                            coeff_T = get(var)[names_T]))
  
  }
   
  for (var in c("coeff_from_nonESRD")){
    names_0 <- intersect(vars_0, rownames(get(var)))
    names_T <- setdiff(rownames(get(var)), names_0)
    assign(paste(var), list(#coeff = get(var), 
                            coeff_0 = get(var)[names_0, ],
                            coeff_T = get(var)[names_T, ]))
    
  }
  
  ### produce baseline lambdas
  
  lam_VD_0_All <- exp(.lamM(beta = coeff_VD$coeff_0, X = dfBaselineMM_0))
  lam_NFMAE_0_All <- exp(.lamM(beta = coeff_NFMAE$coeff_0, X = dfBaselineMM_0))
  lam_NFMVE_0_All <- exp(.lamM(beta = coeff_NFMVE$coeff_0, X = dfBaselineMM_0))
  LP_from_nonESRD_0_All <- list()
  for (var in c("1-3b", "4", "5", "transplant", "dialysis"))
    LP_from_nonESRD_0_All[[var]] <- exp(.lamM(beta = coeff_from_nonESRD$coeff_0[, var], X = dfBaselineMM_0))
  LP_dial_into_tr_0_All <- exp(.lamM(beta = coeff_dial_into_tr$coeff_0, X = dfBaselineMM_0))
  
  return(list(
    lambdas = list(lam_VD_0_All = lam_VD_0_All, 
              lam_NFMAE_0_All = lam_NFMAE_0_All, lam_NFMVE_0_All = lam_NFMVE_0_All,
              LP_dial_into_tr_0_All = LP_dial_into_tr_0_All, LP_from_nonESRD_0_All = LP_from_nonESRD_0_All),
    coeff = list(coeff_VD = coeff_VD, coeff_NFMAE = coeff_NFMAE, coeff_NFMVE = coeff_NFMVE,
                  coeff_dial_into_tr = coeff_dial_into_tr, coeff_from_nonESRD = coeff_from_nonESRD)))


}

### QoL

prelimQoL <- function(vars_0, dfBaselineMM_0, coeff_QoL) {
  
  names_0 <- intersect(vars_0, names(coeff_QoL))
  names_T <- setdiff(names(coeff_QoL), names_0)
  
  coeff_QoL <- list(coeff_0 = coeff_QoL[names_0], coeff_T = coeff_QoL[names_T])
  
  lam_QoL_0 <- .lamM(beta = coeff_QoL$coeff_0, X = dfBaselineMM_0)
  
  return(list(lambdas = list(lam_QoL_0 = lam_QoL_0),
              coeff = list(coeff_QoL = coeff_QoL)))
}
 
### Costs

prelimCost <- function(vars_0, dfBaselineMM_0, coeff_cost) {
    
  names_0 <- intersect(vars_0, names(coeff_cost))
  names_T <- setdiff(names(coeff_cost), names_0)
  
  coeff_cost <- list(coeff_0 = coeff_cost[names_0], coeff_T = coeff_cost[names_T])
  
  lam_cost_0 <- .lamM(beta = coeff_cost$coeff_0, X = dfBaselineMM_0)
  
  return(list(lambdas = list(lam_cost_0 = lam_cost_0),
              coeff = list(coeff_cost = coeff_cost)))
}

###################################################################
#                     State conversion          			            #
###################################################################

get_state_combined <- function(N_CV, N_CKD) {
  return(paste(N_CV, N_CKD, sep = "_"))
}

###################################################################
#                  extracting probabilities            		        #
###################################################################

# calculate \beta*X
# beta is a vector and X is a matrix
# X must contain all columns corresponding to beta
.lamM <- function(beta, X) {
  X <- X[, names(beta)]  # TODO: possibly do this outside the main cycle??
  if (length(beta) == 1)
    X <- as.matrix(X, ncol = 1)
  return(X %*% beta)
}

# beta is a vector and X is a vector
.lamVec <- function(beta, X) {
  n <- intersect(names(beta), names(X))
  return(X[n] %*% beta[n])
}

# cumulative hazard for cox ph distribution
H_coxph <- function(lam_0, beta, X, bh0, bh1) {
  lam <- lam_0 * exp(.lamVec(beta, X))
  return((bh0 - bh1) * lam)
}
  
# cumulative hazard for exponential distribution
H_e <- function(lam_0, beta, X, cycleLengthInDays) {
  lam <- lam_0 * exp(.lamVec(beta, X))
  return(lam * (-cycleLengthInDays))
}

# cumulative hazard for Weibull distribution
H_w <- function(lam_0, beta, X, p, t, cycleLengthInDays) {
  lam <- lam_0 * exp(.lamVec(beta, X))
  return(lam * ((t - cycleLengthInDays)^p - t^p))
}

# cumulative hazard for Gompertz distribution
H_g <- function(lam_0, beta, X, g, t, cycleLengthInDays) {
  lam <- lam_0 * exp(.lamVec(beta, X))
  return((lam / g) * exp(g * t) * (- 1 + exp(-g * cycleLengthInDays))) 
}

# Calculate the transition probability from a survival model 

tp_survival_coxph <- function(lam_0, X, beta, anc, t, cycleLengthInDays) {
  t <- t / cycleLengthInDays
  bh0 <- anc[t]
  bh1 <- anc[t + 1]
  H <- H_coxph(lam_0 = lam_0, beta = beta, X = X, bh0 = bh0, bh1 = bh1) 
  return(1 - exp(H))
}
  
tp_survival_e <- function(lam_0, X, beta, anc = NULL, t, cycleLengthInDays) {
  H <- H_e(lam_0 = lam_0, beta = beta, X = X, cycleLengthInDays = cycleLengthInDays) 
  return(1 - exp(H))
}

tp_survival_w <- function(lam_0, X, beta, anc, t, cycleLengthInDays) {
  H <- H_w(lam_0 = lam_0, beta = beta, X = X, p = anc, t = t, cycleLengthInDays = cycleLengthInDays) 
  return(1 - exp(H))
}

tp_survival_g <- function(lam_0, X, beta, anc, t, cycleLengthInDays) {
  H <- H_g(lam_0 = lam_0, beta = beta, X = X, g = anc, t = t, cycleLengthInDays = cycleLengthInDays)
  return(1 - exp(H))
}

# Calculate the transition probability from a logistic model 
tp_logistic <- function(LP_0, X, beta) {
  # linear predictor
  LP <- LP_0 * exp(.lamVec(beta = beta, X = X))
  # transition probability
  tp <- LP / (1 + LP)
  return(tp)
}

# Calculate the transition probability from a multinomial model 
tp_multinomial <- function(LP_0, X, beta) {
  n <- intersect(names(X), rownames(beta))
  X <- X[n]
  beta <- beta[n, ]
  if (is.vector(beta))
    beta <- t(beta)
  # linear predictors
  LP <- c()
  for (var in colnames(beta))    
    LP[var] <- LP_0[var] * exp(X %*% beta[, var])
  L <- sum(LP)
  tp <- c()
  for (var in names(LP))
    tp[var] <- LP[var] / L
  return(tp)
}

# get cost from the standard linear regression model
c_lm <- function(lam_0, X, beta) {
  # identity link function
  return(lam_0 + .lamVec(beta = beta, X = X))
}

get_pNVD_Ext <- function(ratesNVD, age, N_CKD) {
  return(ratesNVD[ratesNVD[, 1] == N_CKD & ratesNVD[, 2] <= age & age < ratesNVD[, 3], 4])  
}

get_pNFMAE <- function(pNFMAEorVD, pVD) {
  return(max(0, pNFMAEorVD - pVD))
}

get_pNFMVEnotNFMAE <- function(pNFMVEorVD, pNFMAEorVD, pVD) {
  pNFMVEnotNFMAE <- pNFMVEorVD - max(pNFMAEorVD, pVD)
  return(max(0, pNFMVEnotNFMAE))
}

get_pcond <- function(p, pNVD) {
  return(p * (1 - pNVD))
}

########## internal functions for generating the vector of covariates ##########

.getX_CV_CKD <- function(N_CV, vX_state, vX_T) {
  
  vX_T <- c(vX_state, vX_T)
  
  return(vX_T)
}


.getX_QALY <- function(vX_state, vX_T, fatal = FALSE) {
  vX_T <- c(vX_T, age_T2 = as.numeric(vX_T[1]) + 1 * fatal)
  vX_T <- c(vX_state, vX_T)
  return(vX_T)
}

gen_nextCycle <- function(age_T, CKDDuration_T, vP_0, vP_1, 
                          states_info, 
                          states0_i_j,
                          cycleLengthInDays, 
                          endpts_CV, endpts_CKD, endpts_F,  
                          ratesNVD,
                          tp_survival_VD, coeff_VD_T, anc_VD, lam_VD_0,
                          tp_survival_NFMAE, coeff_NFMAE_T, anc_NFMAE, lam_NFMAE_0,
                          tp_survival_NFMVE, coeff_NFMVE_T, anc_NFMVE, lam_NFMVE_0,
                          coeff_from_nonESRD_T, LP_from_nonESRD_0,
                          coeff_dial_into_tr_T, coeff_tr_into_dial, LP_dial_into_tr_0,
                          lamQoL_0, coeff_QoL_T,
                          lamCost_0, coeff_cost_T) {
  
  ########## update relevant covariates ##########
  
  ### update output probability vectors
  
  # 1st output vector - after the CKD model is run
  cycleTempCount <- vP_0[2]
  cycleTemp <- cycleTempCount + 1
  vP_1[2] <- cycleTemp
  
  # 2nd output vector - after the CV model is run
  vP_2 <- vP_1
  
  # generate vector with time-updated covariates
  vX_T <- c(age_T + cycleTempCount, CKDDuration_T + cycleTempCount)
  names(vX_T) = c("age_T", "CKDDuration_T")
  
  ########## run CKD model ##########
  
  ### extract possible start states and outcomes 

  states0_num <- states0_i_j
  
  l0 <- length(states0_num)
  states0 <- states_info[unlist(states0_num)]
  states0_lab <- sapply(states0, '[[', 1)
  
  states1_num <- unique(unlist(lapply(states0, '[[', 15), use.names = FALSE ))  # Ns1
  l1 <- length(states1_num)
  states1 <- states_info[unlist(states1_num)]
  states1_lab <-  sapply(states1, '[[', 1)
  
  states2_num <- unique(unlist(lapply(states0, '[[', 16), use.names = FALSE))  # Ns2
  l2 <- length(states2_num)
  states2 <- states_info[unlist(states2_num)]
  states2_lab <-  sapply(states2, '[[', 1)
  
  ### generate probability matrix
  MCKDrownames <- states0_lab 
  MCKDcolnames <- states1_lab  
  MCKD <- matrix(nrow = length(MCKDrownames), ncol = length(MCKDcolnames),
              dimnames = list(MCKDrownames, MCKDcolnames))
  MCKD[, ] <- 0  
  
  # separate probability matrix for the endpoints
  MCKD_endptsrownames <- states0_lab 
  MCKD_endptscolnames <- endpts_CKD
  MCKD_endpts <- matrix(nrow = length(MCKD_endptsrownames), ncol = length(MCKD_endptscolnames),
                        dimnames = list(MCKD_endptsrownames, MCKD_endptscolnames))
  MCKD_endpts[, ] <- 0  
  
  for (s in 1 : l0) {
    
    state0 <- states0[[s]]
    
    N_CV0 <- state0[[2]]
    N_CKD0 <- state0[[8]]
    vX_state0 <- state0[[12]]
    
    X <- .getX_CV_CKD(N_CV0, vX_state0, vX_T)   
    
    if (N_CKD0 %in% c(1, 2, 3)) {  # from non-ESRD
      # calculate
      tp <- tp_multinomial(LP_0 = LP_from_nonESRD_0, X = X, beta = coeff_from_nonESRD_T)    
      # update probability matrix: final state pre-ESRD
      for (N_CKD1 in c(1, 2, 3)) {
        MCKD[s, get_state_combined(N_CV0, N_CKD1)] <- tp[N_CKD1]
      }
      # update probability matrix: final state ESRD
      pTr <- tp[4]
      MCKD[s, get_state_combined(N_CV0, 4)] <- pTr
      pDial <- tp[5]
      MCKD[s, get_state_combined(N_CV0, 8)] <- pDial
      MCKD_endpts[s, "endpt_CKD_ESRD_first"] <- pTr + pDial
    } else if (N_CKD0 %in% c(8, 9, 10, 11)) {  # from dialysis
        # modelled probability, ie p(dialysis -> transplant)
        pTr <- tp_logistic(LP_0 = LP_dial_into_tr_0, X = X, beta = coeff_dial_into_tr_T)
        N_next <- min(N_CKD0 + 1, 11)
        MCKD[s, get_state_combined(N_CV0, N_next - 4)] <- pTr
        # residual probability, ie p(dialysis -> dialysis)
        MCKD[s, get_state_combined(N_CV0, N_next)] <- 1 - pTr
      } else if (N_CKD0 %in% c(4, 5, 6, 7)) {  # from transplant
          # modelled probability, ie p(transplant -> dialysis)
          pDial <- coeff_tr_into_dial
          N_next <- min(N_CKD0 + 1, 7)
          MCKD[s, get_state_combined(N_CV0, N_next + 4)] <- pDial
          # residual probability, ie p(transplant -> transplant)
          MCKD[s, get_state_combined(N_CV0, N_next)] <- 1 - pDial
      }
  }
  
  ### update v with cumulative probabilities
  
  # only update outcomes that are dealt with in M
  vP_0_temp <- vP_0[MCKDrownames]
  for (outcome in MCKDcolnames){
      vP_1[outcome] <- vP_0_temp %*% MCKD[, outcome]
  }
  for (outcome in MCKD_endptscolnames){
    vP_1[outcome] <- vP_0_temp %*% MCKD_endpts[, outcome]
  }  

  ########## run CV model ##########
  
  ### extract possible end states and outcomes 
  
  ### generate probability matrix
  MCVrownames <- states1_lab  
  MCVcolnames <- c(endpts_CV, states2_lab)
  MCV <- matrix(nrow = length(MCVrownames), ncol = length(MCVcolnames),
              dimnames = list(MCVrownames, MCVcolnames))
  MCV[, ] <- 0
  
  ### calculate T0 and T1 for cumulative hazard   
  
  t <- cycleTemp * cycleLengthInDays 
  #s <- 1
  for (s in 1 : l1) {
    
    state1 <- states1[[s]]
    
    N_CV1 <- state1[[2]]
    N_CKD1 <- state1[[8]]
    N_CKD1_2 <- state1[[11]]
    vX_state1 <- state1[[12]]
    
    X <- .getX_CV_CKD(N_CV1, vX_state1, vX_T)
    
    ### NVD
    pNVD <- get_pNVD_Ext(ratesNVD = ratesNVD, age = vX_T[1], N_CKD = N_CKD1_2)
    
    # update the probability matrix
    MCV[s, "endpt_CV_NVD"] <- pNVD
    # record the number of patients who experienced NVD from a noNFMVE state
    if (N_CV1 == 1)
      MCV[s, "endpt_CV_NVD_noNFMVE"] <- pNVD
    
    ### VD
    
    pVD <- do.call(tp_survival_VD, args = list(lam_0 = lam_VD_0, X = X, beta = coeff_VD_T, anc = anc_VD, 
                                               t = t, cycleLengthInDays = cycleLengthInDays))
    
    p <- get_pcond(p = pVD, pNVD = pNVD)
    # update probability matrix
    MCV[s, "endpt_CV_VD"] <- p
    
    ### NFMAE
    
    # NFMAE or VD
    pNFMAEorVD <- do.call(tp_survival_NFMAE, args = list(lam_0 = lam_NFMAE_0, X = X, beta = coeff_NFMAE_T, anc = anc_NFMAE, 
                                                         t = t, cycleLengthInDays = cycleLengthInDays))
    #p <- get_pcond(p = pNFMAEorVD, pNVD = pNVD)
    
    # NFMAE
    pNFMAE <- get_pNFMAE(pNFMAEorVD = pNFMAEorVD, pVD = pVD)
    p <- get_pcond(p = pNFMAE, pNVD = pNVD)
    MCV[s, get_state_combined(2, N_CKD1)] <- p
  
    ### NFMVE not NFMAE
  
    # NFMVE or VD
    pNFMVEorVD <- do.call(tp_survival_NFMVE, args = list(lam_0 = lam_NFMVE_0, X = X, beta = coeff_NFMVE_T, anc = anc_NFMVE, 
                                                           t = t, cycleLengthInDays = cycleLengthInDays))
      
    # NFMVE not NFMAE
    pNFMVEnotNFMAE <- get_pNFMVEnotNFMAE(pNFMVEorVD = pNFMVEorVD, pNFMAEorVD = pNFMAEorVD, pVD = pVD)
    p <- get_pcond(p = pNFMVEnotNFMAE, pNVD = pNVD)
    #MCV[s, "endpt_CV_NFMVEnotNFMAE"] <- p 
    if (N_CV1 %in% c(1, 5, 6, 7))
      MCV[s, get_state_combined(5, N_CKD1)] <- p  
    
    ### changed on 05/07
    ### probability of first NFMVE or VD
    
    # original code
    #if (N_CV1 == 1) {
    #  p <- get_pcond(p = pNFMVEorVD, pNVD = pNVD)
    #  MCV[s, "endpt_CV_NFMVEorVD_first"] <- p
    #}
    
    if (N_CV1 == 1) {
      p <- get_pcond(p = max(pVD, pNFMAEorVD, pNFMVEorVD), pNVD = pNVD)
      MCV[s, "endpt_CV_NFMVEorVD_first"] <- p
    }
    
    ### end of change
    
    ### no NFMVE
    
    # end state
    if (N_CV1 == 1)  
      N_CV2 <- 1 else 
        if (N_CV1 %in% c(2, 3, 4))
          N_CV2 <- min(N_CV1 + 1, 4) else 
            if (N_CV1 %in% c(5, 6, 7))
              N_CV2 <- min(N_CV1 + 1, 7)
    # residual probability
    MCV[s, get_state_combined(N_CV2, N_CKD1)] <- 1 - sum(MCV[s, c(states2_lab, endpts_F)])
  }
  
  ### update v with cumulative probabilities
    
  # only update outcomes that are dealt with in MCV
  vP_1_temp <- vP_1[MCVrownames]
  for (outcome in MCVcolnames){
    vP_2[outcome] <- vP_1_temp %*% MCV[, outcome]
  }
  for (outcome in MCKD_endptscolnames){
    vP_2[outcome] <- vP_1[outcome]
  } 
  
  ########## transitional probabilities from vP_0 to vP_1 ##########
  
  M <- MCKD %*% MCV
  
  ########## QALYs ##########
  
  ### for a non-fatal state
  
  QALY_NF <- 0
  
  for (s in 1 : l2) {
      
    state2 <- states2[[s]]
    
    state2_lab <- state2[[1]]
    N_CV2 <- state2[[2]]
    N_CKD2 <- state2[[8]]
    vX_state2 <- state2[[12]]
    
    X <- .getX_QALY(vX_state2, vX_T)
    y <- c_lm(lam_0 = lamQoL_0, X = X, beta = coeff_QoL_T)
      
    # convert into QALY and calculate total contribution
    QALY_NF <- QALY_NF + y * vP_2[state2_lab]  # * 1
  }
  vP_2["QALY_NF"] <- QALY_NF
  
  ### for a fatal state
  
  QALY_F <- 0
  
  for (s in 1 : l1) {   
      
    state1 <- states1[[s]]
    
    state1_lab <- state1[[1]]
    N_CV1 <- state1[[2]]
    N_CKD1 <- state1[[8]]
    vX_state1 <- state1[[12]]
  
    X <- .getX_QALY(vX_state1, vX_T, fatal = TRUE)
    y <- c_lm(lam_0 = lamQoL_0, X = X, beta = coeff_QoL_T)
      
    # convert into QALY and calculate total contribution
    # (probability of ending up in this non-fatal state at the end of the CKD model) x (probability of dying when transitioning from this state)
    #multiply by 1/2 assuming the patient died half-way through    
    QALY_F <- QALY_F + y * vP_1[state1_lab] * (MCV[state1_lab, "endpt_CV_NVD"] + MCV[state1_lab, "endpt_CV_VD"]) 
  }
  vP_2["QALY_F"] <- QALY_F / 2
  
  ########## costs ##########
  
  cost <- 0  
  
  for (s0 in 1 : l0) {
      
    state0 <- states0[[s0]]
    
    state0_lab <- state0[[1]]
    
    N_CV0 <- state0[[2]]
    N_CKD0 <- state0[[8]]
    vX_state0 <- state0[[12]]
    
    # non-fatal outcome
    for (s2 in 1 : length(state0$s2)){
   
      alpha <- state0[[14]][[s2]]
      state2 <- states_info[[alpha[[1]]]]
      rs <- alpha[[2]]
      
      state2_lab <- state2[[1]]
      N_CV2 <- state2[[2]]
      N_CKD2 <- state2[[8]]
      vX_state2 <- state2[[12]]        
      
      X <- c(vX_T, vX_state2, age_T2 = as.numeric(vX_T[1]))
      X <- c(X, 1)
      names(X)[length(X)] <- paste("rs", rs, sep = "")
      
      # cost
      cTemp <- c_lm(lam_0 = lamCost_0, X = X, beta = coeff_cost_T)
      cost <- cost + cTemp * vP_0[state0_lab] * M[state0_lab, state2_lab]                      
    }   
    # fatal outcome
    for (s1 in 1 : length(state0$s1)) {

      alpha <- state0[[13]][[s1]]
      state1 <- states_info[[alpha[[1]]]]
      rs <- alpha[[2]]
      rs_lab <- paste("rs", rs, sep = "")
      
      state1_lab <- state1[[1]]
      N_CV1 <- state1[[2]]
      N_CKD1 <- state1[[8]]
      vX_state1 <- state1[[12]]
      
      # vascular death
      X <- c(vX_T, age_T2 = as.numeric(vX_T[1]))
      X <- c(X, rep(1, 2))
      labs <- c("cvd2", rs_lab)
      l <-  length(X)
      names(X)[(l - 1) : l] <- labs
      if (N_CKD1 %in% c(8, 9, 10, 11))
        X <- c(X, death_dialTRUE = 1)
      
      # cost
      cTemp <- c_lm(lam_0 = lamCost_0, X = X, beta = coeff_cost_T)
      cost <- cost + cTemp * vP_0[state0_lab] * MCKD[state0_lab, state1_lab] * MCV[state1_lab, "endpt_CV_VD"] 
      
      # non-vascular death
      X <- c(vX_T, age_T2 = as.numeric(vX_T[1]))
      X <- c(X, rep(1, 2))
      labs <- c("cvd3", rs_lab)
      l <- length(X)
      names(X)[(l - 1) : l] <- labs
      if (N_CKD1 %in% c(8, 9, 10, 11))
        X <- c(X, death_dialTRUE = 1)
      
      cTemp <- c_lm(lam_0 = lamCost_0, X = X, beta = coeff_cost_T)
      cost <- cost + cTemp * vP_0[state0_lab] * MCKD[state0_lab, state1_lab] * MCV[state1_lab, "endpt_CV_NVD"]    
    }
  } 
  
  vP_2["cost_hosp"] <- cost
  
  return(vP_2)
}

###################################################################
# clean-up functions
###################################################################

.cleanup_dfOutput <- function(dfOutput, disc_events, disc_costs, compl, Tx_price) {
  
  # death
  dfOutput$endpt_CV_D <- dfOutput$endpt_CV_NVD + dfOutput$endpt_CV_VD
  
  ### FU times

  # D
  dfOutput$aliveAtEnd <- eval(parse(text = paste("dfOutput[, '", as.vector(outer(c(1 : 7), c(1 : 11), function(x, y) paste(x, y, sep = "_"))), "']", sep = "", collapse = "+")))
  dfOutput$aliveAtStart <- dfOutput$aliveAtEnd + dfOutput$endpt_CV_D
  # VD
  dfOutput$FU_VD <- dfOutput$aliveAtStart - dfOutput$endpt_CV_NVD
  # first NFMVE or VD
  dfOutput$noNFMVEAtEnd <- eval(parse(text = paste("dfOutput[, '", as.vector(outer(c(1), c(1 : 11), function(x, y) paste(x, y, sep = "_"))), "']", sep = "", collapse = "+")))
  dfOutput$noNFMVEAtStart <- c(1, dfOutput$noNFMVEAtEnd[1 : (nrow(dfOutput) - 1)])
  dfOutput$noNFMVEAtStart[dfOutput$cycle == 1] <- 1
  dfOutput$FU_NFMVEorVD_first <- dfOutput$noNFMVEAtStart - dfOutput$endpt_CV_NVD_noNFMVE
  # first ESRD
  dfOutput$preESRDAtEnd <- eval(parse(text = paste("dfOutput[, '", as.vector(outer(c(1 : 7), c(1 : 3), function(x, y) paste(x, y, sep = "_"))), "']", sep = "", collapse = "+")))
  # TODO: better to pass an argument whether the patient is pre-ESRD & replace FU_ESRD_first with 0 or 1 accordingly
  dfOutput$FU_ESRD_first <- c(1, dfOutput$preESRDAtEnd[1 : (nrow(dfOutput) - 1)])
  # remove cycle 0
  dfOutput <- subset(dfOutput, cycle > 0)
  
  # Life expectancy
  dfOutput$LY <- dfOutput$aliveAtEnd + 1/2 * (dfOutput$endpt_CV_NVD + dfOutput$endpt_CV_VD)
  dfOutput$QALY <- dfOutput$QALY_NF + dfOutput$QALY_F

  # discounted values
  for (var in c("LY", "QALY")) 
    dfOutput[, paste(var, "disc", sep = "_")] <- dfOutput[, var] / (1 + disc_events)^dfOutput$cycle
  for (var in c("cost_hosp")) 
    dfOutput[, paste(var, "disc", sep = "_")] <- dfOutput[, var] / (1 + disc_costs)^dfOutput$cycle
  
  # medication cost 
  dfOutput$cost_tx <- dfOutput$LY * compl * Tx_price * 365
  dfOutput$cost_tx_disc <- dfOutput$LY_disc * compl * Tx_price * 365
  
  return(dfOutput)
}

# cleanup dfOutput to produce rates 
.get_rates <- function(dfOutput, byVars, all = FALSE, lifetime = FALSE) {
  
  ### events

  x <- c(byVars, "cycle")
  df_events <- ddply(dfOutput, x, summarise, 
                     NFMVEorVD_first_prod = 1 - sum(endpt_CV_NFMVEorVD_first) / sum(FU_NFMVEorVD_first),
                     ESRD_first_prod = 1 - sum(endpt_CKD_ESRD_first) / sum(FU_ESRD_first),
                     VD_prod = 1 - sum(endpt_CV_VD) / sum(FU_VD),
                     D_prod = 1 - sum(endpt_CV_D) / sum(aliveAtStart))
  df_events <- ddply(df_events, byVars, transform, 
                     NFMVEorVD_first = 1 - cumprod(NFMVEorVD_first_prod),
                     ESRD_first = 1 - cumprod(ESRD_first_prod),
                     VD = 1 - cumprod(VD_prod),
                     D = 1 - cumprod(D_prod))
  # remove .id column if present
  df_events$.id <- NULL
  
  events <- c("NFMVEorVD_first", "ESRD_first", "VD", "D")
  df_events <- subset(df_events, select = c(byVars, "cycle", events))
  # add "id" column if required
  if (all)
    df_events$id <- "all"
  
  # summarise: at 5 & 10 years
  df_events_5 <- subset(df_events, cycle == 5)
  df_events_5 <- df_events_5[, -which(colnames(df_events_5) == "cycle")]
  colnames(df_events_5)[which(colnames(df_events_5) %in% events)] <- 
    paste(colnames(df_events_5)[which(colnames(df_events_5) %in% events)], 5, sep = "_")
  df_events_10 <- subset(df_events, cycle == 10)
  df_events_10 <- df_events_10[, -which(colnames(df_events_10) == "cycle")]
  colnames(df_events_10)[which(colnames(df_events_10) %in% events)] <- 
    paste(colnames(df_events_10)[which(colnames(df_events_10) %in% events)], 10, sep = "_")
  # combine into one dataset; NA produced for those patients with no data at 5/10 years
  df_events_combined <- merge(df_events_5, df_events_10, all.x = TRUE)
  if (!(lifetime)) {
  # if simulation over a fixed number of years, add a summary over this number of years
    max_cycle <- ddply(df_events, byVars, summarise, cycle = max(cycle))
    max_cycle$.id <- NULL
    df_events_LE <- merge(df_events, max_cycle)
    df_events_LE <- df_events_LE[, -which(colnames(df_events_LE) == "cycle")]
    colnames(df_events_LE)[which(colnames(df_events_LE) %in% events)] <- 
      paste(colnames(df_events_LE)[which(colnames(df_events_LE) %in% events)], "all", sep = "_")
    df_events_combined <- merge(df_events_combined, df_events_LE, all.y = TRUE)
  }
  
  ### life expectancy and costs
  
  df_LE <- subset(dfOutput, 
                  select = union(byVars, c("id", "cycle", 
                                           "cost_hosp", "cost_hosp_disc", "cost_tx", "cost_tx_disc",
                                           "LY", "LY_disc", "QALY", "QALY_disc")))
  
  df_LE <- ddply(df_LE, .(id), summarise, cycle = cycle, 
                 cost_hosp = cumsum(cost_hosp), cost_hosp_disc = cumsum(cost_hosp_disc), 
                 cost_tx = cumsum(cost_tx), cost_tx_disc = cumsum(cost_tx_disc),
                 LY = cumsum(LY), LY_disc = cumsum(LY_disc), 
                 QALY = cumsum(QALY), QALY_disc = cumsum(QALY_disc))

  max_cycle <- ddply(df_LE, .(id), summarise, cycle = max(cycle))
  df_LE <- merge(df_LE, max_cycle)
  
  df_LE <- df_LE[, -which(colnames(df_LE) == "cycle")]
  # calculate average LE if summary
  if (all) {
    df_LE <- ddply(df_LE, byVars, summarise, 
                   LY = mean(LY), LY_disc = mean(LY_disc),
                   QALY = mean(QALY), QALY_disc = mean(QALY_disc),
                   cost_hosp = mean(cost_hosp), cost_hosp_disc = mean(cost_hosp_disc),
                   cost_tx = mean(cost_tx), cost_tx_disc = mean(cost_tx_disc))
    # remove .id column
    df_LE$.id <- NULL
    # add "id" column
    df_LE$id <- "all"
  }
  
  ### combine and save
  
  df_all <- merge(df_events_combined, df_LE)
  return(df_all)
}

# add incremental costs, LYs, QALYs and ICERs to a dataset produced by the master function
.add_inc <- function(df) {
  
  # incremental costs, LYs and QALYs
  for (var in c("LY", "QALY", "cost_hosp", "cost_tx")){
    varNew <- paste(var, "inc", sep = "_")
    df[, varNew] <- 
      df[, paste(var, "T", sep = "_")] - df[, paste(var, "C", sep = "_")]
    df[, paste(varNew, "disc", sep = "_")] <- 
      df[, paste(var, "disc", "T", sep = "_")] - df[, paste(var, "disc", "C", sep = "_")]
  }
  
  # total incremental costs
  df[, "cost_total_inc"] <- df[, "cost_hosp_inc"] + df[, "cost_tx_inc"]
  df[, "cost_total_inc_disc"] <- df[, "cost_hosp_inc_disc"] + df[, "cost_tx_inc_disc"]
  
  # ICERs
  df[, "cost_LY"] <- df[, "cost_total_inc"] / df[, "LY_inc"]
  df[, "cost_QALY"] <- df[, "cost_total_inc"] / df[, "QALY_inc"]
  
  df[, "cost_LY_disc"] <- df[, "cost_total_inc_disc"] / df[, "LY_inc_disc"]
  df[, "cost_QALY_disc"] <- df[, "cost_total_inc_disc"] / df[, "QALY_inc_disc"]
  
  return(df)
  
}

# CI for events
.PSA_CI_events <- function(vec, level = 0.95) {
  N <- length(vec)
  vec <- sort(vec)
  temp <- (1 - level) / 2
  lower <- vec[max(floor(N * temp), 1)]
  upper <- vec[min(ceiling(N * (1 - temp)), N)]
  return(list(l = lower, u = upper))
}

# CI for ICERs
.PSA_CI_ICERs <- function(vec_IC, vec_IE, level = 0.95) {
  
  # create the dataset with IC & IE
  df <- data.frame(IC = vec_IC, IE = vec_IE)
  
  # generate ICER
  df$ICER <- df$IC / df$IE
  #df$IE <- as.vector(df$IE)
  #df$IC <- as.vector(df$IC)
  
  # determine quadrant
  df$Q <- "NW"
  df$Q[df$IC > 0 & df$IE > 0] <- "NE"
  df$Q[df$IC < 0 & df$IE > 0] <- "SE"
  df$Q[df$IC < 0 & df$IE < 0] <- "SW"
  
  ### points in the south-west quadrant (cheaper & less effective) - exclude from analyses
  N_SW <- nrow(subset(df, Q == "SW"))
  
  ### points in the south-east quadrant (cheaper & more effective) - best outcome
  df_SE <- subset(df, Q == "SE")
  # order by increasing ICER
  df_SE <- df_SE[order(df_SE$ICER), ]
  #head(df_SE)
  
  ### points in the north-east quadrant (more expensive & more effective) - second best outcome
  df_NE <- subset(df, Q == "NE")
  # order by increasing ICER
  df_NE <- df_NE[order(df_NE$ICER), ]
  #head(df_NE)
  
  ### points in the north-west quadrant (more expensive & less effective) - third best outcome
  df_NW <- subset(df, Q == "NW")
  # order by decreasing ICER (with the same IE, more expensive = smaller negative ICER = worse)
  df_NW <- df_NW[order(-df_NW$ICER), ]
  #head(df_NW)
  
  ### combine
  df_all <- rbind(df_SE, df_NE, df_NW)
  
  ### TODO: what if all point in SW quadrant?
  
  ###  CI
  N <- nrow(df_all)
  temp <- (1 - level) / 2
  N_lower <- max(floor(N * temp), 1) 
  N_upper <- min(ceiling(N * (1 - temp)), N)
  
  # values corresponding to these rows
  ICER_lower <- df_all$ICER[N_lower]
  IE_lower <- df_all$IE[N_lower]
  ICER_upper <- df_all$ICER[N_upper]
  IE_upper <- df_all$IE[N_upper]
  
  # TODO: what exactly should the output be?
  
  #### return value
  #retval <- list(
  #  # number of points in the south-west quadrant
  #  N_SW = N_SW, 
  #  # lower end of CI, together with corresponding IE, to determine which quadrant it comes from
  #  lower = list(ICER = ICER_lower, IE = IE_lower),
  #  # upper end of CI, together with corresponding IE, to determine which quadrant it comes from
  #  upper = list(ICER = ICER_upper, IE = IE_upper)
  #retval <- paste("(", 
  #                round(ICER_lower, digits = digits), ", ", 
  #                round(ICER_upper, digits = digits), ")", sep = "")
  
  return(list(l = ICER_lower, u = ICER_upper))
}

# CI for ICERs
.PSA_CEAC <- function(df, vec_IC_lab, vec_IE_lab, vec_R, N) {
  
  df <- df[, c("id", vec_IC_lab, vec_IE_lab)]
  
  # generate ICER
  df$ICER <- df[, vec_IC_lab] / df[, vec_IE_lab]
  
  temp <- df
  temp$CE <- 0
  
  retval <- list()
  for (R in vec_R) {
    temp2 <- temp
    temp2$CE[temp2[, vec_IE_lab] < 0 & temp2$ICER >= R] <- 1
    temp2$CE[temp2[, vec_IE_lab] >= 0 & temp2$ICER <= R] <- 1
    output <- temp2[, c(1, 5 : ncol(temp2))]
    output <- ddply(output, .(id), numcolwise(sum))
    output$R <- R
    retval <- append(retval, list(output))
  }
  retval <- as.data.frame(do.call("rbind", retval))
  retval$CE <- retval$CE / N
  
  return(retval)
}

###################################################################
# wrappers: master function
###################################################################

master <- function(N_cores, cycleLengthInDays,
                   ids, 
                   dfBaselineMM_0, dfBaselineMM_T, P_0, 
                   stageRand_2, stageRand_2_num, 
                   BVD, CVRand_num,
                   states0, states_info, endpts_CV, endpts_CKD, endpts_F,
                   ratesNVD_M, ratesNVD_M_ESRD, ratesNVD_F, ratesNVD_F_ESRD,
                   tp_survival_VD, anc_VD, coeff_VD,
                   tp_survival_NFMAE, anc_NFMAE, coeff_NFMAE,
                   tp_survival_NFMVE, anc_NFMVE, coeff_NFMVE,
                   coeff_from_nonESRD, coeff_dial_into_tr, coeff_tr_into_dial,  
                   coeff_QoL, coeff_cost, 
                   compl, Tx_price,
                   disc_events, disc_costs,
                   lifetime, years) {
  
  ### add compliance variable to dfBaselineMM_0
  
  dfBaselineMM_0 <- cbind(dfBaselineMM_0, complAll = compl)
  
  ### beta * X for baseline characteristics
  
  vars_0 <- colnames(dfBaselineMM_0)
  lamRE <- prelimRiskEquations(vars_0 = vars_0,
                               dfBaselineMM_0 = dfBaselineMM_0,
                               coeff_VD = coeff_VD, coeff_NFMAE = coeff_NFMAE, coeff_NFMVE = coeff_NFMVE,
                               coeff_dial_into_tr = coeff_dial_into_tr, coeff_from_nonESRD = coeff_from_nonESRD)
  
  lamRE_lambdas <- lamRE$lambdas
  lamRE_coeff <- lamRE$coeff
  
  coeff_VD_T <- lamRE_coeff$coeff_VD$coeff_T
  coeff_NFMAE_T <- lamRE_coeff$coeff_NFMAE$coeff_T
  coeff_NFMVE_T <- lamRE_coeff$coeff_NFMVE$coeff_T
  coeff_from_nonESRD_T <- lamRE_coeff$coeff_from_nonESRD$coeff_T
  coeff_dial_into_tr_T <- lamRE_coeff$coeff_dial_into_tr$coeff_T
  
  lamQoL <- prelimQoL(vars_0 = vars_0, dfBaselineMM_0 = dfBaselineMM_0, coeff_QoL = coeff_QoL)
  lamQoL_lambdas <- lamQoL$lambdas
  lamQoL_coeff <- lamQoL$coeff$coeff_QoL
  coeff_QoL_T <- lamQoL_coeff$coeff_T
  
  lamCost <- prelimCost(vars_0 = vars_0, dfBaselineMM_0 = dfBaselineMM_0, coeff_cost = coeff_cost)  
  lamCost_lambdas <- lamCost$lambdas
  coeff_cost_T <- lamCost$coeff$coeff_cost$coeff_T
  
  # generate 1st output probability vector
  vP_1 <- rep(0, length(colnames(P_0))) 
  names(vP_1) <- colnames(P_0)
  
  #### Start a snowfall cluster and register it with foreach
  sfInit(cpu = N_cores, parallel = TRUE)
  cl <- sfGetCluster()
  clusterExport(cl, c("H_coxph", "H_e", "H_g", "H_w", "c_lm", "gen_nextCycle", 
                    "get_state_combined", ".get_stateN",
                    "get_pNFMAE", "get_pNFMVEnotNFMAE", "get_pNVD_Ext",
                    "get_pcond",
                    "tp_logistic", "tp_multinomial",
                    "tp_survival_coxph", "tp_survival_e", "tp_survival_g", "tp_survival_w",
                    ".lamVec", ".lamM",
                    ".getX_CV_CKD", ".getX_QALY"))  
  registerDoSNOW(cl)
  
  retval <- foreach (i = ids) %dopar% {
  
#retval <- list()  
#for (i in ids) {
    
    ### extract baseline characteristics
    vX_0 <- dfBaselineMM_0[i, ]
    vX_T <- dfBaselineMM_T[i, ]
    age_T <- vX_T[2]
    CKDDuration_T <- vX_T[4]
    
    # subset the lifetable depending on age, gender & baseline ESRD status
      stageRand <- stageRand_2[i]
      if (vX_0[2] == 1) {
        if (stageRand %in% c("dialysis", "transplant"))
          ratesNVD_2 <- ratesNVD_M_ESRD else
            ratesNVD_2 <- ratesNVD_M
        } else {
              if(stageRand %in% c("dialysis", "transplant"))
              ratesNVD_2 <- ratesNVD_F_ESRD else
                ratesNVD_2 <- ratesNVD_F
            }
      ratesNVD_2 <- ratesNVD_2[ratesNVD_2[, 3] >= vX_T[2], ]
    
    ### extract probability vectors
    
    # baseline probabilities
    
    vP_0 <- P_0[i, ]  
    vP_1[1] <- vP_0[1]
    
    ### generate output matrix
    
    if (lifetime)
      years <- floor(years - vX_T[2] + 1)
    output <- matrix(nrow = years, ncol = length(vP_0))   
   
    ### baseline beta * X for risk equations
  
    lam_VD_0 <- lamRE_lambdas$lam_VD_0_All[i] 
    lam_NFMAE_0 <- lamRE_lambdas$lam_NFMAE_0_All[i]  
    lam_NFMVE_0 <- lamRE_lambdas$lam_NFMVE_0_All[i]
    
    LP_from_nonESRD_0 <- c()
    for (n in names(lamRE_lambdas$LP_from_nonESRD_0))
      LP_from_nonESRD_0[n] <- lamRE_lambdas$LP_from_nonESRD_0[[n]][i]
    
    LP_dial_into_tr_0 <- lamRE_lambdas$LP_dial_into_tr_0_All[i]
    
    ### baseline beta * X for QoL
    lamQoL_0 <- lamQoL_lambdas$lam_QoL_0[i]
    
    ### baseline beta * X for costs
    lamCost_0 <- lamCost_lambdas$lam_cost_0[i]

    # information on initial state
    stageRand_i <- stageRand_2_num[i]
    CVRand_i <- CVRand_num[i]
    BVD_i <- BVD[i]
    
    state_num <- .get_stateN(states_info = states_info, BVD = BVD_i, 
                             N_CV = CVRand_i, N_CKD = stageRand_i)
     
    states0_i <- states0[[state_num]]
    
    for (j in 1 : years) {        
      states0_i_j <- states0_i[[min(j, 5)]]      

      vP_0 <- gen_nextCycle(age_T = age_T, CKDDuration_T = CKDDuration_T,
                            vP_0 = vP_0, vP_1 = vP_1, 
                            states_info = states_info,
                            states0_i_j = states0_i_j,
                            cycleLengthInDays = cycleLengthInDays,
                            endpts_CV = endpts_CV, endpts_CKD = endpts_CKD, endpts_F = endpts_F,
                            ratesNVD = ratesNVD_2, 
                            tp_survival_VD = tp_survival_VD, coeff_VD_T = coeff_VD_T, anc_VD = anc_VD, lam_VD_0 = lam_VD_0,
                            tp_survival_NFMAE = tp_survival_NFMAE, coeff_NFMAE_T = coeff_NFMAE_T, anc_NFMAE = anc_NFMAE, lam_NFMAE_0 = lam_NFMAE_0,
                            tp_survival_NFMVE = tp_survival_NFMVE, coeff_NFMVE_T = coeff_NFMVE_T, anc_NFMVE = anc_NFMVE, lam_NFMVE_0 = lam_NFMVE_0,
                            coeff_from_nonESRD_T = coeff_from_nonESRD_T, LP_from_nonESRD_0 = LP_from_nonESRD_0,
                            coeff_dial_into_tr_T = coeff_dial_into_tr_T, LP_dial_into_tr_0 = LP_dial_into_tr_0,
                            coeff_tr_into_dial = coeff_tr_into_dial, 
                            lamQoL_0 = lamQoL_0, coeff_QoL_T = coeff_QoL_T,
                            lamCost_0 = lamCost_0, coeff_cost_T = coeff_cost_T
      )
      
      output[j, ] <- vP_0 
    }
 
    # save the output
    return(output)
    
#retval[[i]] <- output
    
  }
  sfStop()

  ### combine output
  
  dfOutput <- as.data.frame(do.call("rbind", retval))
  colnames(dfOutput) <- colnames(P_0)

  # TODO: I think this is only needed to know whether the patient was pre-ESRD at start
  # an alternative is to pass this in as an argument
  dfOutput <- rbind(dfOutput, P_0[ids, ])
  dfOutput <- dfOutput[order(dfOutput$id, dfOutput$cycle), ]
  
  # produce rates
  dfOutput <-.cleanup_dfOutput(dfOutput = dfOutput, 
                               disc_events = disc_events, disc_costs = disc_costs, 
                               compl = compl, Tx_price = Tx_price)
  
  df_ind <- .get_rates(dfOutput = dfOutput, byVars = "id", all = FALSE, lifetime = lifetime)
  df_group <- .get_rates(dfOutput = dfOutput, byVars = c(), all = TRUE, lifetime = lifetime)
  
  # return
  return(list(df_ind = df_ind, df_group = df_group))
  
}

master_nonpar <- function(cycleLengthInDays,
                   ids, 
                   dfBaselineMM_0, dfBaselineMM_T, P_0, 
                   stageRand_2, stageRand_2_num, 
                   BVD, CVRand_num,
                   states0, states_info, endpts_CV, endpts_CKD, endpts_F,
                   ratesNVD_M, ratesNVD_M_ESRD, ratesNVD_F, ratesNVD_F_ESRD,
                   tp_survival_VD, anc_VD, coeff_VD,
                   tp_survival_NFMAE, anc_NFMAE, coeff_NFMAE,
                   tp_survival_NFMVE, anc_NFMVE, coeff_NFMVE,
                   coeff_from_nonESRD, coeff_dial_into_tr, coeff_tr_into_dial,  
                   coeff_QoL, coeff_cost, 
                   compl, Tx_price,
                   disc_events, disc_costs,
                   lifetime, years) {
  
  ### add compliance variable to dfBaselineMM_0
  
  dfBaselineMM_0 <- cbind(dfBaselineMM_0, complAll = compl)
  
  ### beta * X for baseline characteristics
  
  vars_0 <- colnames(dfBaselineMM_0)
  lamRE <- prelimRiskEquations(vars_0 = vars_0,
                               dfBaselineMM_0 = dfBaselineMM_0,
                               coeff_VD = coeff_VD, coeff_NFMAE = coeff_NFMAE, coeff_NFMVE = coeff_NFMVE,
                               coeff_dial_into_tr = coeff_dial_into_tr, coeff_from_nonESRD = coeff_from_nonESRD)
  
  lamRE_lambdas <- lamRE$lambdas
  lamRE_coeff <- lamRE$coeff
  
  coeff_VD_T <- lamRE_coeff$coeff_VD$coeff_T
  coeff_NFMAE_T <- lamRE_coeff$coeff_NFMAE$coeff_T
  coeff_NFMVE_T <- lamRE_coeff$coeff_NFMVE$coeff_T
  coeff_from_nonESRD_T <- lamRE_coeff$coeff_from_nonESRD$coeff_T
  coeff_dial_into_tr_T <- lamRE_coeff$coeff_dial_into_tr$coeff_T
  
  lamQoL <- prelimQoL(vars_0 = vars_0, dfBaselineMM_0 = dfBaselineMM_0, coeff_QoL = coeff_QoL)
  lamQoL_lambdas <- lamQoL$lambdas
  lamQoL_coeff <- lamQoL$coeff$coeff_QoL
  coeff_QoL_T <- lamQoL_coeff$coeff_T
  
  lamCost <- prelimCost(vars_0 = vars_0, dfBaselineMM_0 = dfBaselineMM_0, coeff_cost = coeff_cost)  
  lamCost_lambdas <- lamCost$lambdas
  coeff_cost_T <- lamCost$coeff$coeff_cost$coeff_T
  
  # generate 1st output probability vector
  vP_1 <- rep(0, length(colnames(P_0))) 
  names(vP_1) <- colnames(P_0)
  
  # execute the simulation
  retval <- list()  
  for (i in ids) {
    
    ### extract baseline characteristics
    vX_0 <- dfBaselineMM_0[i, ]
    vX_T <- dfBaselineMM_T[i, ]
    age_T <- vX_T[2]
    CKDDuration_T <- vX_T[4]
    
    # subset the lifetable depending on age, gender & baseline ESRD status
    stageRand <- stageRand_2[i]
    if (vX_0[2] == 1) {
      if (stageRand %in% c("dialysis", "transplant"))
        ratesNVD_2 <- ratesNVD_M_ESRD else
          ratesNVD_2 <- ratesNVD_M
    } else {
      if(stageRand %in% c("dialysis", "transplant"))
        ratesNVD_2 <- ratesNVD_F_ESRD else
          ratesNVD_2 <- ratesNVD_F
    }
    ratesNVD_2 <- ratesNVD_2[ratesNVD_2[, 3] >= vX_T[2], ]
    
    ### extract probability vectors
    
    # baseline probabilities
    
    vP_0 <- P_0[i, ]  
    vP_1[1] <- vP_0[1]
    
    ### generate output matrix
    
    if (lifetime)
      years <- floor(years - vX_T[2] + 1)
    output <- matrix(nrow = years, ncol = length(vP_0))   
    
    ### baseline beta * X for risk equations
    
    lam_VD_0 <- lamRE_lambdas$lam_VD_0_All[i] 
    lam_NFMAE_0 <- lamRE_lambdas$lam_NFMAE_0_All[i]  
    lam_NFMVE_0 <- lamRE_lambdas$lam_NFMVE_0_All[i]
    
    LP_from_nonESRD_0 <- c()
    for (n in names(lamRE_lambdas$LP_from_nonESRD_0))
      LP_from_nonESRD_0[n] <- lamRE_lambdas$LP_from_nonESRD_0[[n]][i]
    
    LP_dial_into_tr_0 <- lamRE_lambdas$LP_dial_into_tr_0_All[i]
    
    ### baseline beta * X for QoL
    lamQoL_0 <- lamQoL_lambdas$lam_QoL_0[i]
    
    ### baseline beta * X for costs
    lamCost_0 <- lamCost_lambdas$lam_cost_0[i]
    
    # information on initial state
    stageRand_i <- stageRand_2_num[i]
    CVRand_i <- CVRand_num[i]
    BVD_i <- BVD[i]
    
    state_num <- .get_stateN(states_info = states_info, BVD = BVD_i, 
                             N_CV = CVRand_num[i], N_CKD = stageRand_i)
    
    states0_i <- states0[[state_num]]
    
    for (j in 1 : years) {        
      states0_i_j <- states0_i[[min(j, 5)]]      
      
      vP_0 <- gen_nextCycle(age_T = age_T, CKDDuration_T = CKDDuration_T,
                            vP_0 = vP_0, vP_1 = vP_1, 
                            states_info = states_info,
                            states0_i_j = states0_i_j,
                            cycleLengthInDays = cycleLengthInDays,
                            endpts_CV = endpts_CV, endpts_CKD = endpts_CKD, endpts_F = endpts_F,
                            ratesNVD = ratesNVD_2, 
                            tp_survival_VD = tp_survival_VD, coeff_VD_T = coeff_VD_T, anc_VD = anc_VD, lam_VD_0 = lam_VD_0,
                            tp_survival_NFMAE = tp_survival_NFMAE, coeff_NFMAE_T = coeff_NFMAE_T, anc_NFMAE = anc_NFMAE, lam_NFMAE_0 = lam_NFMAE_0,
                            tp_survival_NFMVE = tp_survival_NFMVE, coeff_NFMVE_T = coeff_NFMVE_T, anc_NFMVE = anc_NFMVE, lam_NFMVE_0 = lam_NFMVE_0,
                            coeff_from_nonESRD_T = coeff_from_nonESRD_T, LP_from_nonESRD_0 = LP_from_nonESRD_0,
                            coeff_dial_into_tr_T = coeff_dial_into_tr_T, LP_dial_into_tr_0 = LP_dial_into_tr_0,
                            coeff_tr_into_dial = coeff_tr_into_dial, 
                            lamQoL_0 = lamQoL_0, coeff_QoL_T = coeff_QoL_T,
                            lamCost_0 = lamCost_0, coeff_cost_T = coeff_cost_T
      )
      
      output[j, ] <- vP_0 
    }
    
    retval[[i]] <- output
    
  }
  
  ### combine output
  
  dfOutput <- as.data.frame(do.call("rbind", retval))
  colnames(dfOutput) <- colnames(P_0)
  
  # TODO: I think this is only needed to know whether the patient was pre-ESRD at start
  # an alternative is to pass this in as an argument
  dfOutput <- rbind(dfOutput, P_0[ids, ])
  dfOutput <- dfOutput[order(dfOutput$id, dfOutput$cycle), ]
  
  # produce rates
  dfOutput <-.cleanup_dfOutput(dfOutput = dfOutput, 
                               disc_events = disc_events, disc_costs = disc_costs, 
                               compl = compl, Tx_price = Tx_price)
  
  df_ind <- .get_rates(dfOutput = dfOutput, byVars = "id", all = FALSE, lifetime = lifetime)
  df_group <- .get_rates(dfOutput = dfOutput, byVars = c(), all = TRUE, lifetime = lifetime)
  
  # return
  return(list(df_ind = df_ind, df_group = df_group))
  
}


###################################################################
# wrappers: long-term projection
###################################################################

# wrapper for deterministic analysis
wrapper_LE <- function(N_cores,
                       dfBaselineMM_0, dfBaselineMM_T, 
                       states_and_endpoints, stageRand_2, stageRand_2_num, 
                       BVD, CVRand_num,
                       ratesNVD, P_0, coeffs, compl, Tx_price,
                       disc_events, disc_costs,
                       lifetime, years) {
  
  ### clean and prepare the data
  
  ids <- 1: nrow(dfBaselineMM_0)
  N_cores <- min(nrow(dfBaselineMM_0), N_cores)
  
  states0 <- states_and_endpoints$starting_states
  states_info <- states_and_endpoints$states_info
  endpts_CV <- states_and_endpoints$endpts_CV
  endpts_CKD <- states_and_endpoints$endpts_CKD
  endpts_F <- states_and_endpoints$endpts_F
  
  # non-vascular death
  ratesNVD_M <- ratesNVD$ratesNVD_M
  ratesNVD_M_ESRD <- ratesNVD$ratesNVD_M_ESRD
  ratesNVD_F <- ratesNVD$ratesNVD_F
  ratesNVD_F_ESRD <- ratesNVD$ratesNVD_F_ESRD
  
  #other
  cycleLengthInDays <- 365
  
  # survival distributions
  tp_survival_VD <-  "tp_survival_e"
  anc_VD <- NULL
  tp_survival_NFMAE <- "tp_survival_g"
  tp_survival_NFMVE <- "tp_survival_g"
  
  # vascular death
  coeff_VD <- coeffs[["VD"]]
  
  # NFMAE or VD
  temp <- coeffs[["NFMAE_or_VD2"]]
  coeff_NFMAE <- temp[[1]]
  anc_NFMAE <- temp[[2]]
  
  # NFMAE or VD
  temp <- coeffs[["NFMVE_or_VD2"]]
  coeff_NFMVE <- temp[[1]]
  anc_NFMVE <- temp[[2]]
  
  # transition from pre-ESRD
  coeff_from_nonESRD <- coeffs[["nonesrd2"]]
  
  # dialysis into transplant
  coeff_dial_into_tr <- coeffs[["dial2trans"]]
  
  # transplant into dialysis
  coeff_tr_into_dial <- coeffs[["trans2dial"]]
  
  # quality of life
  coeff_QoL <- coeffs[["qol2"]]
  
  # costs
  coeff_cost <- coeffs[["cost2"]]
  
  ### run analyses
  temp <- master(N_cores = N_cores, cycleLengthInDays = cycleLengthInDays,
                   ids = ids, 
                   dfBaselineMM_0 = dfBaselineMM_0, 
                   dfBaselineMM_T = dfBaselineMM_T, P_0 = P_0, 
                   stageRand_2 = stageRand_2, stageRand_2_num = stageRand_2_num, 
                   BVD = BVD, CVRand_num = CVRand_num,
                   states0 = states0, states_info = states_info, 
                   endpts_CV = endpts_CV, endpts_CKD = endpts_CKD, endpts_F = endpts_F,
                   ratesNVD_M = ratesNVD_M, ratesNVD_M_ESRD = ratesNVD_M_ESRD, 
                   ratesNVD_F = ratesNVD_M, ratesNVD_F_ESRD = ratesNVD_F_ESRD,
                   tp_survival_VD = tp_survival_VD, anc_VD = anc_VD, coeff_VD = coeff_VD,
                   tp_survival_NFMAE = tp_survival_NFMAE, anc_NFMAE = anc_NFMAE, coeff_NFMAE = coeff_NFMAE,
                   tp_survival_NFMVE = tp_survival_NFMVE, anc_NFMVE = anc_NFMVE, coeff_NFMVE = coeff_NFMVE,
                   coeff_from_nonESRD = coeff_from_nonESRD, coeff_dial_into_tr = coeff_dial_into_tr, 
                   coeff_tr_into_dial = coeff_tr_into_dial,  
                   coeff_QoL = coeff_QoL, coeff_cost = coeff_cost, 
                   compl = compl, Tx_price = Tx_price,
                   disc_events = disc_events, disc_costs = disc_costs,
                   lifetime = lifetime, years = years)
  
  # summary at the individual level
  output_ind <- temp$df_ind
  # summary at group level
  output_group <- temp$df_group
  
  # return
  return(list(output_ind = output_ind, output_group = output_group))
}

# wrapper for PSA
wrapper_LE_PSA <- function(N_cores, 
                           dfBaselineMM_0, dfBaselineMM_T, 
                           states_and_endpoints, stageRand_2, stageRand_2_num, 
                           BVD = BVD, CVRand_num = CVRand_num,
                           ratesNVD, P_0, 
                           coeffs_PSA_VD,
                           coeffs_PSA_NFMAEorVD, coeffs_PSA_NFMAEorVD_gamma,
                           coeffs_PSA_NFMVEorVD, coeffs_PSA_NFMVEorVD_gamma,
                           coeffs_PSA_nonesrd, coeffs_PSA_dial2trans, coeffs_PSA_trans2dial,
                           coeffs_PSA_cost, coeffs_PSA_qol,
                           compl, Tx_price,
                           sims,
                           disc_events, disc_costs,
                           lifetime, years) {
  
  ### clean and prepare the data
  
  ids <- 1: nrow(dfBaselineMM_0)
  N_cores <- min(length(sims), N_cores)
  
  states0 <- states_and_endpoints$starting_states
  states_info <- states_and_endpoints$states_info
  endpts_CV <- states_and_endpoints$endpts_CV
  endpts_CKD <- states_and_endpoints$endpts_CKD
  endpts_F <- states_and_endpoints$endpts_F
  
  # non-vascular death
  ratesNVD_M <- ratesNVD$ratesNVD_M
  ratesNVD_M_ESRD <- ratesNVD$ratesNVD_M_ESRD
  ratesNVD_F <- ratesNVD$ratesNVD_F
  ratesNVD_F_ESRD <- ratesNVD$ratesNVD_F_ESRD
  
  #other
  cycleLengthInDays <- 365
  
  # survival distributions
  tp_survival_VD <-  "tp_survival_e"
  anc_VD <- NULL
  tp_survival_NFMAE <- "tp_survival_g"
  tp_survival_NFMVE <- "tp_survival_g"
  
  sfInit(cpu = N_cores, parallel = TRUE)
  cl <- sfGetCluster()
  clusterExport(cl, 
                c("H_coxph", "H_e", "H_g", "H_w", "c_lm", "gen_nextCycle", 
                  "get_state_combined", ".get_stateN",
                  "get_pNFMAE", "get_pNFMVEnotNFMAE", "get_pNVD_Ext",
                  "get_pcond",
                  "tp_logistic", "tp_multinomial",
                  "tp_survival_coxph", "tp_survival_e", "tp_survival_g", "tp_survival_w",
                  ".lamVec", ".lamM",
                  ".getX_CV_CKD", ".getX_QALY", 
                  "prelimRiskEquations", "prelimCost", "prelimQoL", 
                  ".cleanup_dfOutput", ".get_rates", "master_nonpar"))  
  registerDoSNOW(cl)
  
  retval <- foreach (n = sims, .packages = c("plyr")) %dopar% {
    
    # vascular death
    coeff_VD <- coeffs_PSA_VD[n, ]
    
    # NFMAE or VD
    coeff_NFMAE <- coeffs_PSA_NFMAEorVD[n, ]
    anc_NFMAE <- coeffs_PSA_NFMAEorVD_gamma[n]
    
    # NFMVE or VD
    coeff_NFMVE <- coeffs_PSA_NFMVEorVD[n, ]
    anc_NFMVE <- coeffs_PSA_NFMVEorVD_gamma[n]
    
    # transition from pre-ESRD
    coeff_from_nonESRD <- coeffs_PSA_nonesrd[[n]]
    
    # dialysis into transplant
    coeff_dial_into_tr <- coeffs_PSA_dial2trans[n, ]
    
    # transplant into dialysis
    coeff_tr_into_dial <- coeffs_PSA_trans2dial[n]
    
    # quality of life
    coeff_QoL <- coeffs_PSA_qol[n, ]
    
    # costs
    coeff_cost <- coeffs_PSA_cost[n, ]
    
    ### run analyses
    
    temp <- master_nonpar(cycleLengthInDays = cycleLengthInDays,
                     ids = ids, 
                     dfBaselineMM_0 = dfBaselineMM_0, 
                     dfBaselineMM_T = dfBaselineMM_T, P_0 = P_0, 
                     stageRand_2 = stageRand_2, stageRand_2_num = stageRand_2_num, 
                     BVD = BVD, CVRand_num = CVRand_num,
                     states0 = states0, states_info = states_info, 
                     endpts_CV = endpts_CV, endpts_CKD = endpts_CKD, endpts_F = endpts_F,
                     ratesNVD_M = ratesNVD_M, ratesNVD_M_ESRD = ratesNVD_M_ESRD, 
                     ratesNVD_F = ratesNVD_M, ratesNVD_F_ESRD = ratesNVD_F_ESRD,
                     tp_survival_VD = tp_survival_VD, anc_VD = anc_VD, coeff_VD = coeff_VD,
                     tp_survival_NFMAE = tp_survival_NFMAE, anc_NFMAE = anc_NFMAE, coeff_NFMAE = coeff_NFMAE,
                     tp_survival_NFMVE = tp_survival_NFMVE, anc_NFMVE = anc_NFMVE, coeff_NFMVE = coeff_NFMVE,
                     coeff_from_nonESRD = coeff_from_nonESRD, coeff_dial_into_tr = coeff_dial_into_tr, 
                     coeff_tr_into_dial = coeff_tr_into_dial,  
                     coeff_QoL = coeff_QoL, coeff_cost = coeff_cost, 
                     compl = compl, Tx_price = Tx_price, 
                     disc_events = disc_events, disc_costs = disc_costs,
                     lifetime = lifetime, years = years)
    
    # summary at the individual level
    output_ind <- temp$df_ind
    # summary at group level
    output_group <- temp$df_group
    # combine
    output <- rbind(output_ind, output_group)
    # add simulation number
    output$sim <- n
    
    return(output)
    #retval[[n]] <- output
  }
  
  sfStop()
  # rbind
  output <- as.data.frame(do.call("rbind", retval))
  
  # calcualte confidence intervals
  
  # CIs for events
  dfCI <- data.frame(setDT(output)[, 
                                   as.list(c(D_5 = .PSA_CI_events(D_5),
                                             VD_5 = .PSA_CI_events(VD_5),
                                             NFMVEorVD_first_5 = .PSA_CI_events(NFMVEorVD_first_5),
                                             ESRD_first_5 = .PSA_CI_events(ESRD_first_5),
                                             D_10 = .PSA_CI_events(D_10),
                                             VD_10 = .PSA_CI_events(VD_10),
                                             NFMVEorVD_first_10 = .PSA_CI_events(NFMVEorVD_first_10),
                                             ESRD_first_10 = .PSA_CI_events(ESRD_first_10),
                                             LY = .PSA_CI_events(LY),
                                             QALY = .PSA_CI_events(QALY),
                                             cost_hosp = .PSA_CI_events(cost_hosp),
                                             cost_tx = .PSA_CI_events(cost_tx))), .(id)])
  if (!(lifetime))
    dfCI <- merge(dfCI, data.frame(setDT(output)[, 
                                     as.list(c(D_all = .PSA_CI_events(D_all),
                                               VD_all = .PSA_CI_events(VD_all),
                                               NFMVEorVD_first_all = .PSA_CI_events(NFMVEorVD_first_all),
                                               ESRD_first_all = .PSA_CI_events(ESRD_first_all))), .(id)]))

  return(list(output = output, dfCI = dfCI))
  
}


###################################################################
# wrappers: cost-effectiveness analysis
###################################################################

# wrapper for deterministic analysis
wrapper <- function(N_cores,
                    dfBaselineMM_0, dfBaselineMM_T, 
                    states_and_endpoints, stageRand_2, stageRand_2_num, 
                    BVD, CVRand_num,
                    ratesNVD, P_0, coeffs,
                    compl_C, compl_T,
                    Tx_price_C, Tx_price_T,
                    disc_events, disc_costs,
                    lifetime, years) {
  
  ### clean and prepare the data
  
  ids <- 1: nrow(dfBaselineMM_0)
  N_cores <- min(nrow(dfBaselineMM_0), N_cores)
  
  states0 <- states_and_endpoints$starting_states
  states_info <- states_and_endpoints$states_info
  endpts_CV <- states_and_endpoints$endpts_CV
  endpts_CKD <- states_and_endpoints$endpts_CKD
  endpts_F <- states_and_endpoints$endpts_F
  
  # non-vascular death
  ratesNVD_M <- ratesNVD$ratesNVD_M
  ratesNVD_M_ESRD <- ratesNVD$ratesNVD_M_ESRD
  ratesNVD_F <- ratesNVD$ratesNVD_F
  ratesNVD_F_ESRD <- ratesNVD$ratesNVD_F_ESRD
  
  #other
  cycleLengthInDays <- 365
  
  # survival distributions
  tp_survival_VD <-  "tp_survival_e"
  anc_VD <- NULL
  tp_survival_NFMAE <- "tp_survival_g"
  tp_survival_NFMVE <- "tp_survival_g"
  
  # vascular death
  coeff_VD <- coeffs[["VD"]]
  
  # NFMAE or VD
  temp <- coeffs[["NFMAE_or_VD2"]]
  coeff_NFMAE <- temp[[1]]
  anc_NFMAE <- temp[[2]]
  
  # NFMAE or VD
  temp <- coeffs[["NFMVE_or_VD2"]]
  coeff_NFMVE <- temp[[1]]
  anc_NFMVE <- temp[[2]]
  
  # transition from pre-ESRD
  coeff_from_nonESRD <- coeffs[["nonesrd2"]]
  
  # dialysis into transplant
  coeff_dial_into_tr <- coeffs[["dial2trans"]]
  
  # transplant into dialysis
  coeff_tr_into_dial <- coeffs[["trans2dial"]]
  
  # quality of life
  coeff_QoL <- coeffs[["qol2"]]
  
  # costs
  coeff_cost <- coeffs[["cost2"]]
  
  ### run analyses
  
  # treatment A
  temp_C <- master(N_cores = N_cores, cycleLengthInDays = cycleLengthInDays,
                  ids = ids, 
                  dfBaselineMM_0 = dfBaselineMM_0, 
                  dfBaselineMM_T = dfBaselineMM_T, P_0 = P_0, 
                  stageRand_2 = stageRand_2, stageRand_2_num = stageRand_2_num, 
                  BVD = BVD, CVRand_num = CVRand_num,
                  states0 = states0, states_info = states_info, 
                  endpts_CV = endpts_CV, endpts_CKD = endpts_CKD, endpts_F = endpts_F,
                  ratesNVD_M = ratesNVD_M, ratesNVD_M_ESRD = ratesNVD_M_ESRD, 
                  ratesNVD_F = ratesNVD_M, ratesNVD_F_ESRD = ratesNVD_F_ESRD,
                  tp_survival_VD = tp_survival_VD, anc_VD = anc_VD, coeff_VD = coeff_VD,
                  tp_survival_NFMAE = tp_survival_NFMAE, anc_NFMAE = anc_NFMAE, coeff_NFMAE = coeff_NFMAE,
                  tp_survival_NFMVE = tp_survival_NFMVE, anc_NFMVE = anc_NFMVE, coeff_NFMVE = coeff_NFMVE,
                  coeff_from_nonESRD = coeff_from_nonESRD, coeff_dial_into_tr = coeff_dial_into_tr, 
                  coeff_tr_into_dial = coeff_tr_into_dial,  
                  coeff_QoL = coeff_QoL, coeff_cost = coeff_cost, 
                  compl = compl_C, Tx_price = Tx_price_C, 
                  disc_events = disc_events, disc_costs = disc_costs,
                  lifetime = lifetime, years = years)
  
  # treatment B
  temp_T <- master(N_cores = N_cores, cycleLengthInDays = cycleLengthInDays,
                  ids = ids, 
                  dfBaselineMM_0 = dfBaselineMM_0, 
                  dfBaselineMM_T = dfBaselineMM_T, P_0 = P_0, 
                  stageRand_2 = stageRand_2, stageRand_2_num = stageRand_2_num, 
                  BVD = BVD, CVRand_num = CVRand_num,
                  states0 = states0, states_info = states_info, 
                  endpts_CV = endpts_CV, endpts_CKD = endpts_CKD, endpts_F = endpts_F,
                  ratesNVD_M = ratesNVD_M, ratesNVD_M_ESRD = ratesNVD_M_ESRD, 
                  ratesNVD_F = ratesNVD_M, ratesNVD_F_ESRD = ratesNVD_F_ESRD,
                  tp_survival_VD = tp_survival_VD, anc_VD = anc_VD, coeff_VD = coeff_VD,
                  tp_survival_NFMAE = tp_survival_NFMAE, anc_NFMAE = anc_NFMAE, coeff_NFMAE = coeff_NFMAE,
                  tp_survival_NFMVE = tp_survival_NFMVE, anc_NFMVE = anc_NFMVE, coeff_NFMVE = coeff_NFMVE,
                  coeff_from_nonESRD = coeff_from_nonESRD, coeff_dial_into_tr = coeff_dial_into_tr, 
                  coeff_tr_into_dial = coeff_tr_into_dial,  
                  coeff_QoL = coeff_QoL, coeff_cost = coeff_cost, 
                  compl = compl_T, Tx_price = Tx_price_T, 
                  disc_events = disc_events, disc_costs = disc_costs,
                  lifetime = lifetime, years = years)
  
  # summary at the individual level
  output_ind <- merge(temp_C$df_ind, temp_T$df_ind, by = "id", suffixes = c("_C", "_T"))
  output_ind <- .add_inc(output_ind)
  
  # summary at the group level
  output_group <- merge(temp_C$df_group, temp_T$df_group, by = "id", suffixes = c("_C", "_T"))
  output_group <- .add_inc(output_group)
  
  # return
  return(list(output_ind = output_ind, output_group = output_group))
}

# wrapper for PSA
wrapper_PSA <- function(N_cores, 
                        dfBaselineMM_0, dfBaselineMM_T, 
                        states_and_endpoints, stageRand_2, stageRand_2_num, 
                        BVD = BVD, CVRand_num = CVRand_num,
                        ratesNVD, P_0, 
                        coeffs_PSA_VD,
                        coeffs_PSA_NFMAEorVD, coeffs_PSA_NFMAEorVD_gamma,
                        coeffs_PSA_NFMVEorVD, coeffs_PSA_NFMVEorVD_gamma,
                        coeffs_PSA_nonesrd, coeffs_PSA_dial2trans, coeffs_PSA_trans2dial,
                        coeffs_PSA_cost, coeffs_PSA_qol,
                        sims,
                        compl_C, compl_T,
                        Tx_price_C, Tx_price_T,
                        disc_events, disc_costs,
                        lifetime, years) {
  
  ### clean and prepare the data
  
  ids <- 1: nrow(dfBaselineMM_0)
  N_cores <- min(length(sims), N_cores)
  
  states0 <- states_and_endpoints$starting_states
  states_info <- states_and_endpoints$states_info
  endpts_CV <- states_and_endpoints$endpts_CV
  endpts_CKD <- states_and_endpoints$endpts_CKD
  endpts_F <- states_and_endpoints$endpts_F
  
  # non-vascular death
  ratesNVD_M <- ratesNVD$ratesNVD_M
  ratesNVD_M_ESRD <- ratesNVD$ratesNVD_M_ESRD
  ratesNVD_F <- ratesNVD$ratesNVD_F
  ratesNVD_F_ESRD <- ratesNVD$ratesNVD_F_ESRD
  
  #other
  cycleLengthInDays <- 365
  
  # survival distributions
  tp_survival_VD <-  "tp_survival_e"
  anc_VD <- NULL
  tp_survival_NFMAE <- "tp_survival_g"
  tp_survival_NFMVE <- "tp_survival_g"

  sfInit(cpu = N_cores, parallel = TRUE)
  cl <- sfGetCluster()
  clusterExport(cl, 
                c("H_coxph", "H_e", "H_g", "H_w", "c_lm", "gen_nextCycle", 
                  "get_state_combined", ".get_stateN",
                  "get_pNFMAE", "get_pNFMVEnotNFMAE", "get_pNVD_Ext",
                  "get_pcond",
                  "tp_logistic", "tp_multinomial",
                  "tp_survival_coxph", "tp_survival_e", "tp_survival_g", "tp_survival_w",
                  ".lamVec", ".lamM",
                  ".getX_CV_CKD", ".getX_QALY", 
                  "prelimRiskEquations", "prelimCost", "prelimQoL", 
                  ".cleanup_dfOutput", ".get_rates", "master_nonpar"))  
  registerDoSNOW(cl)
  
  retval <- foreach (n = sims, .packages = c("plyr")) %dopar% {
    
    # vascular death
    coeff_VD <- coeffs_PSA_VD[n, ]
    
    # NFMAE or VD
    coeff_NFMAE <- coeffs_PSA_NFMAEorVD[n, ]
    anc_NFMAE <- coeffs_PSA_NFMAEorVD_gamma[n]
    
    # NFMVE or VD
    coeff_NFMVE <- coeffs_PSA_NFMVEorVD[n, ]
    anc_NFMVE <- coeffs_PSA_NFMVEorVD_gamma[n]
    
    # transition from pre-ESRD
    coeff_from_nonESRD <- coeffs_PSA_nonesrd[[n]]
    
    # dialysis into transplant
    coeff_dial_into_tr <- coeffs_PSA_dial2trans[n, ]
    
    # transplant into dialysis
    coeff_tr_into_dial <- coeffs_PSA_trans2dial[n]
    
    # quality of life
    coeff_QoL <- coeffs_PSA_qol[n, ]
    
    # costs
    coeff_cost <- coeffs_PSA_cost[n, ]
    
    ### run analyses
  
    # treatment A
    temp_C <- master_nonpar(cycleLengthInDays = cycleLengthInDays,
                       ids = ids, 
                       dfBaselineMM_0 = dfBaselineMM_0, 
                       dfBaselineMM_T = dfBaselineMM_T, P_0 = P_0, 
                       stageRand_2 = stageRand_2, stageRand_2_num = stageRand_2_num, 
                       BVD = BVD, CVRand_num = CVRand_num,
                       states0 = states0, states_info = states_info, 
                       endpts_CV = endpts_CV, endpts_CKD = endpts_CKD, endpts_F = endpts_F,
                       ratesNVD_M = ratesNVD_M, ratesNVD_M_ESRD = ratesNVD_M_ESRD, 
                       ratesNVD_F = ratesNVD_M, ratesNVD_F_ESRD = ratesNVD_F_ESRD,
                       tp_survival_VD = tp_survival_VD, anc_VD = anc_VD, coeff_VD = coeff_VD,
                       tp_survival_NFMAE = tp_survival_NFMAE, anc_NFMAE = anc_NFMAE, coeff_NFMAE = coeff_NFMAE,
                       tp_survival_NFMVE = tp_survival_NFMVE, anc_NFMVE = anc_NFMVE, coeff_NFMVE = coeff_NFMVE,
                       coeff_from_nonESRD = coeff_from_nonESRD, coeff_dial_into_tr = coeff_dial_into_tr, 
                       coeff_tr_into_dial = coeff_tr_into_dial,  
                       coeff_QoL = coeff_QoL, coeff_cost = coeff_cost, 
                       compl = compl_C, Tx_price = Tx_price_C, 
                       disc_events = disc_events, disc_costs = disc_costs,
                       lifetime = lifetime, years = years)
  
    # treatment B
    temp_T <- master_nonpar(cycleLengthInDays = cycleLengthInDays,
                      ids = ids, 
                       dfBaselineMM_0 = dfBaselineMM_0, 
                       dfBaselineMM_T = dfBaselineMM_T, P_0 = P_0, 
                       stageRand_2 = stageRand_2, stageRand_2_num = stageRand_2_num, 
                       BVD = BVD, CVRand_num = CVRand_num,
                       states0 = states0, states_info = states_info, 
                       endpts_CV = endpts_CV, endpts_CKD = endpts_CKD, endpts_F = endpts_F,
                       ratesNVD_M = ratesNVD_M, ratesNVD_M_ESRD = ratesNVD_M_ESRD, 
                       ratesNVD_F = ratesNVD_M, ratesNVD_F_ESRD = ratesNVD_F_ESRD,
                       tp_survival_VD = tp_survival_VD, anc_VD = anc_VD, coeff_VD = coeff_VD,
                       tp_survival_NFMAE = tp_survival_NFMAE, anc_NFMAE = anc_NFMAE, coeff_NFMAE = coeff_NFMAE,
                       tp_survival_NFMVE = tp_survival_NFMVE, anc_NFMVE = anc_NFMVE, coeff_NFMVE = coeff_NFMVE,
                       coeff_from_nonESRD = coeff_from_nonESRD, coeff_dial_into_tr = coeff_dial_into_tr, 
                       coeff_tr_into_dial = coeff_tr_into_dial,  
                       coeff_QoL = coeff_QoL, coeff_cost = coeff_cost, 
                       compl = compl_T, Tx_price = Tx_price_T, 
                       disc_events = disc_events, disc_costs = disc_costs,
                       lifetime = lifetime, years = years)
    
    output_group <- merge(temp_C$df_group, temp_T$df_group, by = "id", suffixes = c("_C", "_T"))
    output_ind <- merge(temp_C$df_ind, temp_T$df_ind, by = "id", suffixes = c("_C", "_T"))
    output <- rbind(output_group, output_ind)
    
    # add simulation number
    output$sim <- n
    
    return(output)
  }
  
  sfStop()
  
  # rbind
  output <- as.data.frame(do.call("rbind", retval))
  
  # add incremental columns and ICERs
  output <- .add_inc(output)
  
  # CIs for events
  dfCI_events <- data.frame(setDT(output)[, 
                                   as.list(c(
                                     # Control group
                                     D_5_C = .PSA_CI_events(D_5_C),
                                     VD_5_C = .PSA_CI_events(VD_5_C),
                                     NFMVEorVD_first_5_C = .PSA_CI_events(NFMVEorVD_first_5_C),
                                     ESRD_first_5_C = .PSA_CI_events(ESRD_first_5_C),
                                     D_10_C = .PSA_CI_events(D_10_C),
                                     VD_10_C = .PSA_CI_events(VD_10_C),
                                     NFMVEorVD_first_10_C = .PSA_CI_events(NFMVEorVD_first_10_C),
                                     ESRD_first_10_C = .PSA_CI_events(ESRD_first_10_C),
                                     LY_C = .PSA_CI_events(LY_C),
                                     LY_disc_C = .PSA_CI_events(LY_disc_C),
                                     QALY_C = .PSA_CI_events(QALY_C),
                                     QALY_disc_C = .PSA_CI_events(QALY_disc_C),
                                     cost_hosp_C = .PSA_CI_events(cost_hosp_C),
                                     cost_hosp_disc_C = .PSA_CI_events(cost_hosp_disc_C),
                                     cost_tx_C = .PSA_CI_events(cost_tx_C),
                                     cost_tx_disc_C = .PSA_CI_events(cost_tx_disc_C),
                                     # Treatment group
                                     D_5_T = .PSA_CI_events(D_5_T),
                                     VD_5_T = .PSA_CI_events(VD_5_T),
                                     NFMVEorVD_first_5_T = .PSA_CI_events(NFMVEorVD_first_5_T),
                                     ESRD_first_5_T = .PSA_CI_events(ESRD_first_5_T),
                                     D_10_T = .PSA_CI_events(D_10_T),
                                     VD_10_T = .PSA_CI_events(VD_10_T),
                                     NFMVEorVD_first_10_T = .PSA_CI_events(NFMVEorVD_first_10_T),
                                     ESRD_first_10_T = .PSA_CI_events(ESRD_first_10_T),
                                     LY_T = .PSA_CI_events(LY_T),
                                     LY_disc_T = .PSA_CI_events(LY_disc_T),
                                     QALY_T = .PSA_CI_events(QALY_T),
                                     QALY_disc_T = .PSA_CI_events(QALY_disc_T),
                                     cost_hosp_T = .PSA_CI_events(cost_hosp_T),
                                     cost_hosp_disc_T = .PSA_CI_events(cost_hosp_disc_T),
                                     cost_tx_T = .PSA_CI_events(cost_tx_T),
                                     cost_tx_disc_T = .PSA_CI_events(cost_tx_disc_T),
                                     # incremental
                                     LY_inc = .PSA_CI_events(LY_inc),
                                     LY_inc_disc = .PSA_CI_events(LY_inc_disc),
                                     QALY_inc = .PSA_CI_events(QALY_inc),
                                     QALY_inc_disc = .PSA_CI_events(QALY_inc_disc),
                                     cost_hosp_inc = .PSA_CI_events(cost_hosp_inc),
                                     cost_hosp_inc_disc = .PSA_CI_events(cost_hosp_inc_disc),
                                     cost_tx_inc = .PSA_CI_events(cost_tx_inc),
                                     cost_tx_inc_disc = .PSA_CI_events(cost_tx_inc_disc),
                                     cost_total_inc = .PSA_CI_events(cost_total_inc),
                                     cost_total_inc_disc = .PSA_CI_events(cost_total_inc_disc))), .(id)])
  if (!(lifetime))
    dfCI_events <- merge(dfCI_events, data.frame(setDT(output)[, 
                                                 as.list(c(
                                                   # Control group
                                                   D_all_C = .PSA_CI_events(D_all_C),
                                                   VD_all_C = .PSA_CI_events(VD_all_C),
                                                   NFMVEorVD_first_all_C = .PSA_CI_events(NFMVEorVD_first_all_C),
                                                   ESRD_first_all_C = .PSA_CI_events(ESRD_first_all_C),
                                                   # Control group
                                                   D_all_T = .PSA_CI_events(D_all_T),
                                                   VD_all_T = .PSA_CI_events(VD_all_T),
                                                   NFMVEorVD_first_all_T = .PSA_CI_events(NFMVEorVD_first_all_T),
                                                   ESRD_first_all_T = .PSA_CI_events(ESRD_first_all_T))), .(id)]))
  # CIs for ICERs
  dfCI_ICERs <- data.frame(setDT(output)[, 
                           as.list(c(
                             cost_LY = .PSA_CI_ICERs(vec_IC = cost_total_inc, vec_IE =  LY_inc),
                             cost_QALY = .PSA_CI_ICERs(vec_IC = cost_total_inc, vec_IE = QALY_inc),
                             cost_LY_disc = .PSA_CI_ICERs(vec_IC = cost_total_inc_disc, vec_IE = LY_inc_disc),
                             cost_QALY_disc = .PSA_CI_ICERs(vec_IC = cost_total_inc_disc, vec_IE =  QALY_inc_disc))), .(id)])

  # combine
  dfCI <- merge(dfCI_events, dfCI_ICERs)
  
  # CEAC
  vec_R <- seq(from = 10000, to = 100000, by = 50)
  N <- length(sims)
  temp <- as.data.frame(subset(output, id == "all"))
  dfCEAC_undisc <- .PSA_CEAC(df = temp, 
                  vec_IC_lab = "cost_total_inc", vec_IE = "QALY_inc", 
                  vec_R = vec_R, N = N)
  dfCEAC_disc <- .PSA_CEAC(df = temp, 
                      vec_IC_lab = "cost_total_inc_disc", vec_IE = "QALY_inc_disc", 
                      vec_R = vec_R, N = N)
  dfCEAC <- merge(dfCEAC_undisc, dfCEAC_disc, by = c("id", "R"), suffixes = c("_undisc", "_disc"))
  
  return(list(output = output,
              dfCI = dfCI, dfCEAC = dfCEAC))

}