rm(list = ls())

# setup the directory and check that the required R packages have been installed

#setwd("L:/SHARP analyses/Shiny Model/") # Set this folder to be the one that contains server.R
setwd("/home/kidneymodel/ShinyApps/app")

# load functions and packages

source("code/master code.R")
source("code/shiny functions.R")

# load defaults values

set.seed(1234)
coeffs_default <- get(load("data/coeffs_default.Rdata"))
coeffs_PSA <- get(load("data/coeffs_PSA.Rdata"))
coeffs_PSA_VD <- coeffs_PSA$VD
coeffs_PSA_NFMAEorVD <- coeffs_PSA$NFMAEorVD
coeffs_PSA_NFMAEorVD_gamma <- coeffs_PSA$NFMAEorVD_gamma
coeffs_PSA_NFMVEorVD <- coeffs_PSA$NFMVEorVD
coeffs_PSA_NFMVEorVD_gamma <- coeffs_PSA$NFMVEorVD_gamma
coeffs_PSA_nonesrd <- coeffs_PSA$nonesrd
coeffs_PSA_dial2trans <- coeffs_PSA$dial2trans
coeffs_PSA_trans2dial <- coeffs_PSA$trans2dial
coeffs_PSA_cost <- coeffs_PSA$cost
coeffs_PSA_qol <- coeffs_PSA$qol

states_and_endpoints <- get(load("data/states_and_endpoints.Rdata"))

pat0 <- get(load("data/pat0.Rdata"))
data <- get(load("data/default_patient.Rdata"))
patT <- get(load("data/patT.Rdata"))
pat0_names <- names(pat0)
patT_names <- names(patT)

stages <- c("1-3b", "4", "5", "dialysis", "transplant")
names(stages) <- c("CKD 3B", "CKD 4", "CKD 5, not RRT", "Dialysis", "Transplant")

CV <- c("None",
        "Major atherosclerotic event in the last year",
        "Major atherosclerotic event 1-2 years ago",
        "Major atherosclerotic event >2 years ago",
        "No MAE, but haemorrhagic stroke in the last year",
        "No MAE, but haemorrhagic stroke 1-2 years ago",
        "No MAE, but haemorrhagic stroke >2 years ago",
        "Another cardiovascular event")

ratesNVD <- get(load("data/ratesNVD.Rdata"))

lifetime <- FALSE
years <- NULL

# max number of cores to be used in parallel programming

N_cores <- 8

################################################################################
### validation: setting up
################################################################################

### id tab

validate_id_screen <- function(x) {
  validate(
    need(x$age == TRUE & x$diab == TRUE & x$ckd_dur == TRUE & x$esrd_dur == TRUE, 
         "The model cannot be executed. Please check the following conditions:"),
    need(x$age == TRUE, "Age should be between 40 and 90"),
    need(x$diab == TRUE, "Patients with diabetic nephropathy must have a diagnosis of diabetes"),
    need(x$ckd_dur == TRUE, "CKD duration should be between 0 and the age value"),
    need(x$esrd_dur == TRUE, "RRT duration should be between 0 and the CKD duration value")
  )
}

req_colnames_id <- c("id", "age", "sex", "ethnicity",
                     "education", "adultDep", "smoker", "currentAlc",
                     "BMI_quant", "DBP_quant", "SBP_quant", 
                     "CholHDL_quant", "Albumin_quant", "Hemoglobin_quant", "Phosphate_quant",
                     "ACR_quant", "CVD", "DM", "CKDStage", "CKDDuration",
                     "renalDiagnosis", "RRTDuration", "TX")

req_format_id <- list(
  list(var = "id", numeric = TRUE),
  list(var = "age", numeric = TRUE),
  list(var = "sex", numeric = TRUE),
  list(var = "ethnicity", numeric = TRUE),
  list(var = "education", numeric = TRUE),
  #list(var = "childDep", numeric = TRUE),
  list(var = "adultDep", numeric = TRUE),
  list(var = "smoker", numeric = TRUE),
  list(var = "currentAlc", numeric = TRUE),
  list(var = "BMI_quant", numeric = TRUE), 
  list(var = "DBP_quant", numeric = TRUE),
  list(var = "SBP_quant", numeric = TRUE),
  list(var = "CholHDL_quant", numeric = TRUE),
  list(var = "Albumin_quant", numeric = TRUE),
  list(var = "Hemoglobin_quant", numeric = TRUE),
  list(var = "Phosphate_quant", numeric = TRUE),
  list(var = "ACR_quant", numeric = TRUE),
  list(var = "CVD", numeric = TRUE),
  list(var = "DM", numeric = TRUE),
  list(var = "CKDStage", numeric = TRUE),
  list(var = "CKDDuration", numeric = TRUE),
  list(var = "renalDiagnosis", numeric = TRUE),
  list(var = "RRTDuration", numeric = TRUE),
  list(var = "TX", numeric = TRUE)
)

req_values_id <- list(
  list(var = "id", othervars = NULL, 
       condition = "anyDuplicated(df$id) == 0", 
       error_text = "Patient id should be unique"),
  list(var = "age", othervars = NULL,
       condition = "df$age >= 40 && df$age <= 90",
       error_text = "age column can only take values between 40 and 90"),
  list(var = "sex", othervars = NULL,
       condition = "all(df$sex %in% c(0, 1))",
       error_text = "sex column can only take values 0, 1"),
  list(var = "ethnicity", othervars = NULL,
       condition = "all(df$ethnicity %in% c(0, 1, 2, 3, 4))",
       error_text = "ethnicity column can only take values 0, 1, 2, 3, 4"),
  list(var = "education", othervars = NULL,
       condition = "all(df$education %in% c(0, 1, 2))",
       error_text = "education column can only take values 0, 1, 2"),
  #list(var = "childDep", othervars = NULL,
  #     condition = "all(df$childDep %in% c(0, 1))",
  #     error_text = "childDep column can only take values 0, 1"),
  list(var = "adultDep", othervars = NULL,
       condition = "all(df$adultDep %in% c(0, 1))",
       error_text = "adultDep column can only take values 0, 1"),
  list(var = "smoker", othervars = NULL,
       condition = "all(df$smoker %in% c(0, 1, 2))",
       error_text = "smoker column can only take values 0, 1, 2"),
  list(var = "currentAlc", othervars = NULL,
       condition = "all(df$currentAlc %in% c(0, 1))",
       error_text = "currentAlc column can only take values 0, 1"),
  list(var = "BMI_quant", othervars = NULL,
       condition = "all(df$BMI_quant %in% c(0, 1, 2))",
       error_text = "BMI_quant column can only take values 0, 1, 2"),
  list(var = "DBP_quant", othervars = NULL,
       condition = "all(df$DBP_quant %in% c(0, 1, 2))",
       error_text = "DBP_quant column can only take values 0, 1, 2"),
  list(var = "SBP_quant", othervars = NULL,
       condition = "all(df$SBP_quant %in% c(0, 1, 2))",
       error_text = "SBP_quant column can only take values 0, 1, 2"),
  list(var = "CholHDL_quant", othervars = NULL,
       condition = "all(df$CholHDL_quant %in% c(0, 1, 2))",
       error_text = "CholHDL_quant column can only take values 0, 1, 2"),
  list(var = "Albumin_quant", othervars = NULL,
       condition = "all(df$Albumin_quant %in% c(0, 1, 2))",
       error_text = "Albumin_quant column can only take values 0, 1, 2"),
  list(var = "Hemoglobin_quant", othervars = NULL,
       condition = "all(df$Hemoglobin_quant %in% c(0, 1, 2))",
       error_text = "Hemoglobin_quant column can only take values 0, 1, 2"),
  list(var = "Phosphate_quant", othervars = NULL,
       condition = "all(df$Phosphate_quant %in% c(0, 1, 2))",
       error_text = "Phosphate_quant column can only take values 0, 1, 2"),
  list(var = "ACR_quant", othervars = c("CKDStage"),
       condition = "all(df$ACR_quant[df$CKDStage %in% c(0, 1, 2)] %in% c(0, 1, 2)) & 
       all(df$ACR_quant[df$CKDStage %in% c(3, 4)] %in% c(3))",
       error_text = "ACR_quant column can only take values 0, 1, 2 for pre-RRT participants and 3 for RRT participants"),
  list(var = "CVD", othervars = NULL,
       condition = "all(df$CVD %in% c(0, 1, 2, 3, 4, 5, 6, 7))",
       error_text = "CVD column can only take values 0, 1, 2, 3, 4, 5, 6, 7"),
  list(var = "DM", othervars = c("renalDiagnosis"),
       condition = "all(df$DM %in% c(0, 1)) & 
       all(df$DM[df$renalDiagnosis == 0] == 1)",
       error_text = "DM column can only take values 0, 1. Participants with diabetic nepropathy should be marked as having diabetes"),
  list(var = "CKDStage", othervars = NULL,
       condition = "all(df$CKDStage %in% c(0, 1, 2, 3, 4))",
       error_text = "CKD stage column can only take values 0, 1, 2, 3, 4"),
  list(var = "CKDDuration", othervars = "age",
       condition = "df$CKDDuration >= 0 && df$CKDDuration <= df$age",
       error_text = "CKDDuration column values should be between 0 and the participant's age"),
  list(var = "renalDiagnosis", othervars = NULL,
       condition = "all(df$renalDiagnosis %in% c(0, 1, 2))",
       error_text = "renalDiagnosis column values can only be 0, 1, 2"),
  list(var = "RRTDuration", othervars = c("CKDStage", "CKDDuration"),
       condition = "all(df$RRTDuration[df$CKDStage %in% c(0, 1, 2)] == 0) &
       all(df$RRTDuration[df$CKDStage %in% c(3, 4)] >= 0) &
       all(df$RRTDuration[df$CKDStage %in% c(3, 4)] <= df$CKDDuration[df$CKDStage %in% c(3, 4)])",
       error_text = "RRTDuration column values should be 0 for pre-RRT participants and between 0 and CKD duration for RRT participants"),
  list(var = "TX", othervars = "CKDStage",
       condition = "all(df$TX[df$CKDStage %in% c(0, 1, 2)] == 0) &
       all(df$TX[df$CKDStage %in% c(3, 4)] %in% c(0, 1))",
       error_text = "TX column can only take values 0, 1. Pre-RRT participants cannot have previous transplants.")
)

validate_file <- function(vars_missing, vars_format, vars_values) {
  validate(
    need(length(vars_missing) == 0 & length(vars_format) == 0 & length(vars_values) == 0,
         "The model cannot be executed. Please check the following conditions:"),
    need(length(vars_missing) == 0, paste("The following columns are missing:", 
                                          paste(vars_missing, sep = "", collapse = "; "))),
    need(length(vars_format) == 0, paste("The following columns are in the wrong format:", 
                                         paste(vars_format, sep = "", collapse = "; "))),
    need(length(vars_values) == 0, paste("The following columns contain disallowed values:", 
                                         paste(vars_values, sep = "", collapse = "; ")))
  ) 
}

### treatment tab

validate_trt <- function(x) {
  validate(
    need(x$RR == TRUE & x$compl1 == TRUE & x$price1 == TRUE, 
         "The model cannot be executed. Please check the following conditions:"),
    need(x$RR == TRUE, "Hazard ratios should be above 0"),
    need(x$compl1 == TRUE, "Compliance should be a percentage between 0 and 100"),
    need(x$price1 == TRUE, "Treatment cost should be non-negative")
  )
}

validate_trt_PSA <- function(x) {
  validate(
    need(x$RR == TRUE & x$RR_CI & x$compl1 == TRUE & x$price1 == TRUE, 
         "The model cannot be executed. Please check the following conditions:"),
    need(x$RR == TRUE, "Hazard ratios should be above 0"),
    need(x$RR_CI == TRUE, "Confidence intervals should contain the mean"),
    need(x$compl1 == TRUE, "Compliance should be a percentage between 0 and 100"),
    need(x$price1 == TRUE, "Treatment cost should be non-negative")
  )
}

### cost tab

validate_cost <- function(x) {
  validate(
    need(x$costs == TRUE, 
         "The model cannot be executed. Costs should be be non-negative")
  )
}

validate_cost_PSA <- function(x) {
  validate(
    need(x$costs == TRUE & x$SEs == TRUE, 
         "The model cannot be executed. Please check the following conditions:"),
    need(x$costs == TRUE, "Costs should be be non-negative"),
    need(x$SEs == TRUE, "Standard errors should be positive")
  )
}

### qol tab

validate_qol <- function(x) {
  validate(
    need(x$qol_min_check == TRUE & x$qol_min == TRUE & x$qol_max == TRUE, 
         "The model cannot be executed. Please check the following conditions:"),
    need(x$qol_min_check == TRUE, 
         "The pre-specified minimum QoL value should be between -2 and 0"),
    need(x$qol_min == TRUE,
         paste("Quality of life cannot be lower than ", x$x_lower, sep = "")),
    need(x$qol_max == TRUE, 
         "Quality of life cannot exceed 1")
  )
}

### nvd tab

req_colnames_nvd <- c("CKDStage", "ageBand", "sex", "p_NVD")

req_format_nvd <- list(
  list(var = "CKDStage", numeric = TRUE),
  list(var = "ageBand", numeric = FALSE),
  list(var = "sex", numeric = TRUE),
  list(var = "p_NVD", numeric = TRUE)
)

req_values_nvd <- list(
  list(var = "CKDStage", othervars = NULL, 
       condition = "setequal(df$CKDStage,  c(0, 1, 2, 3, 4))", 
       error_text = "CKD stage must contain each of 0, 1, 2, 3, 4, and nothing else"),
  list(var = "SEX", othervars = NULL, 
       condition = "setequal(df$sex,  c(0, 1))", 
       error_text = "Sex column must contain each of 0, 1, and nothing else"),
  list(var = "p_NVD", othervars = NULL, 
       condition = "all(df$p_NVD >= 0) & all(df$p_NVD <= 1)", 
       error_text = "Probabilities must take values between 0 and 1")
)

### nvd tab

req_colnames_dp <- c("age")

req_format_dp <- list(
  list(var = "age", numeric = TRUE)
)

req_values_dp <- list(
  list(var = "age", othervars = NULL,
       condition = "df$age >= 40 && df$age <= 90",
       error_text = "Age values should be between 40 and 90")
)

validate_dp_screen <- function(x) {
  if (x$fixed == TRUE)
    validate(
      need(x$cycles == TRUE, 
           paste("The model cannot be executed. The number of cycles should be between 0 and ", 
                 floor(100 - x$age), 
                 ". Predictions cannot be made beyond 100 years of age", sep = ""))
    ) else
      validate(need(x$stopping == TRUE, 
                    paste("The model cannot be executed. The stopping age should be between the patient's age and 100."))
      )
}

validate_dp_file <- function(x) {
  # check the patient input file is in the correct format
  validate(need(x$file == TRUE, 
                "The model cannot be executed. Check that the age column in the patient characteristics file is present, in a numeric format and contains values between 40 and 90"))
  # if all ok, proceed to actual validation
  if (x$fixed == TRUE)
    validate(
      need(x$cycles == TRUE, 
           paste("The model cannot be executed. The number of cycles should be between 0 and ", 
                 floor(100 - x$age), 
                 ". Predictions cannot be made beyond 100 years of age", sep = ""))
    ) else
      validate(need(x$stopping == TRUE, 
                    paste("The model cannot be executed. The stopping age should be between the patient's age and 100."))
      )
}

### server function

shinyServer(function(input, output, session) {
  
  ################################################################################
  ### Switch between tabs
  ################################################################################
  
  observeEvent(input$link_to_Glossary, {
    updateTabItems(session, "panels", "Glossary")
  })
  
  observeEvent(input$link_to_Results, {
    updateTabItems(session, "panels", "Results")
  })
  
  ################################################################################
  ### validation: extract values
  ################################################################################
  
  ### patient characteristics
  
  # screen values
  id_screen_value <- reactive({
    if (is.null(input$slider_age))
      return(list(
        age = TRUE, diab = TRUE, ckd_dur = TRUE, esrd_dur = TRUE
      )) else
        return(list(
          age = (input$slider_age >= 40 & input$slider_age <= 90),
          diab = (input$sel_rdiag != "Diabetic nephropathy" |
                    (input$sel_rdiag == "Diabetic nephropathy" & input$sel_diab == "Yes")),
          ckd_dur = (input$sel_ckddur >= 0 & input$sel_ckddur <= input$slider_age),
          esrd_dur = ((input$sel_ckd %in% c("Dialysis", "Transplant") & 
                         input$sel_esrddur >= 0 & input$sel_esrddur < input$sel_ckddur) | input$sel_ckd %in% c("CKD 3B", "CKD 4", "CKD 5, not RRT"))
        ))
  })
  
  output$validated_id_screen_value <- renderUI({
    validate_id_screen(id_screen_value())
  })
  
  # file
  id_file_value <- reactive({
    
    inFile <- input$file1
    
    
    if (is.null(inFile))
      return(NULL)
    
    alpha <- .validate_file_nonempty(inFile = inFile,
                                     req_colnames = req_colnames_id,
                                     req_format = req_format_id,
                                     req_values = req_values_id)
    return(list(vars_missing = alpha$vars_missing, 
                vars_format = alpha$vars_format, 
                vars_values = alpha$vars_values))
  })
  
  output$validated_id_file_value <- renderUI({
    validate_file(id_file_value()$vars_missing, 
                  id_file_value()$vars_format, 
                  id_file_value()$vars_values)
  })
  
  ### Treatment parameteres
  
  trt_value <- reactive({
    if (is.null(input$RR_VD))
      return(list(
        RR = TRUE, compl1 = TRUE, price1 = TRUE)) else
          return(list(
            RR = (input$RR_VD > 0 & input$RR_NFMAE > 0 & input$RR_NFMVE > 0),
            compl1 = (input$compl1 >= 0 & input$compl1 <= 100),
            price1 = (input$price1 >= 0)
          ))
  })
  
  trt_value_PSA <- reactive({
    if (is.null(input$RR_VD_2))
      return(list(
        RR = TRUE, RR_CI = TRUE, compl1 = TRUE, price1 = TRUE)) else
          return(list(
            RR = (input$RR_VD_2 > 0 & input$RR_NFMAE_2 > 0 & input$RR_NFMVE_2 > 0 &
                    input$RR_VD_l > 0 & input$RR_NFMAE_l > 0 & input$RR_NFMVE_l > 0 & 
                    input$RR_VD_r > 0 & input$RR_NFMAE_r > 0 & input$RR_NFMVE_r > 0),
            RR_CI = (input$RR_VD_l <= input$RR_VD_2 & input$RR_VD_2 <= input$RR_VD_r &
                       input$RR_NFMAE_l <= input$RR_NFMAE_2 & input$RR_NFMAE_2 <= input$RR_NFMAE_r &
                       input$RR_NFMVE_l <= input$RR_NFMVE_2 & input$RR_NFMVE_2 <= input$RR_NFMVE_r),
            compl1 = (input$compl1_2 >= 0 & input$compl1_2 <= 100),
            price1 = (input$price1_2 >= 0)
          ))
  })
  
  output$validated_trt_value <- renderUI({
    validate_trt(trt_value())
  })
  
  output$validated_trt_value_PSA <- renderUI({
    validate_trt_PSA(trt_value_PSA())
  })
  
  ### costs
  
  cost_value <- reactive({
    if (is.null(input$cost_ckd13b))
      return(list(
        costs = TRUE)) else
          return(list(
            costs = (input$cost_ckd13b >= 0 & input$cost_ckd4 >= 0 & input$cost_ckd5 >= 0 & 
                       input$cost_dial0 >= 0 & input$cost_dial1 >= 0 & input$cost_tx0 >= 0 & 
                       input$cost_tx1 >= 0 & input$cost_cvd5_nd >= 0 & input$cost_cvd5_d >= 0 & 
                       input$cost_cvd6 >= 0 & input$cost_cvd478 >= 0 & input$cost_vd >= 0 & 
                       input$cost_nvd >= 0 & input$cost_diab >= 0)
          ))
  })
  
  cost_value_PSA <- reactive({ 
    if (is.null(input$cost_ckd13b_2))
      return(list(
        costs = TRUE, SEs = TRUE)) else
          return(list(
            costs = (input$cost_ckd13b_2 >= 0 & input$cost_ckd4_2 >= 0 & input$cost_ckd5_2 >= 0 & 
                       input$cost_dial0_2 >= 0 & input$cost_dial1_2 >= 0 & input$cost_tx0_2 >= 0 & 
                       input$cost_tx1_2 >= 0 & input$cost_cvd5_nd_2 >= 0 & input$cost_cvd5_d_2 >= 0 & 
                       input$cost_cvd6_2 >= 0 & input$cost_cvd478_2 >= 0 & input$cost_vd_2 >= 0 & 
                       input$cost_nvd_2 >= 0 & input$cost_diab_2 >= 0),
            SEs = (input$cost_ckd13b_sd > 0 & input$cost_ckd4_sd > 0 & input$cost_ckd5_sd > 0 & 
                     input$cost_dial0_sd > 0 & input$cost_dial1_sd > 0 & input$cost_tx0_sd > 0 & 
                     input$cost_tx1_sd > 0 & input$cost_cvd5_nd_sd > 0 & input$cost_cvd5_d_sd > 0 & 
                     input$cost_cvd6_sd > 0 & input$cost_cvd478_sd > 0 & input$cost_vd_sd > 0 & 
                     input$cost_nvd_sd > 0 & input$cost_diab_sd > 0)
          ))
  })
  
  output$validated_cost_value <- renderUI({
    validate_cost(cost_value())
  })
  
  output$validated_cost_value_PSA <- renderUI({
    validate_cost_PSA(cost_value_PSA())
  })
  
  ### QoL
  
  qol_value <- reactive({
    if (is.null(input$qol_baseline))
      return(list(
        x = -1.54, 
        qol_min_check = TRUE, 
        qol_min = TRUE, 
        qol_max = TRUE
      )) else {
        x <- input$qol_min
        # TODO: save coordinates to avoid recalculation
        # minimal possible QoL
        qol_min <- input$qol_baseline + 
          (3 * input$qol_age) * (input$qol_age <= 0) + 
          (2 * input$qol_age) * (input$qol_age >= 0) +
          (input$qol_sexm) * (input$qol_sexm <= 0) +
          min((input$qol_eduGCSE) * (input$qol_eduGCSE <= 0), 
              (input$qol_eduBelowSec) * (input$qol_eduBelowSec <= 0)) +
          min((input$qol_smokerbefore) * (input$qol_smokerbefore <= 0), 
              (input$qol_smokercurrently) * (input$qol_smokercurrently <= 0)) +
          min((input$qol_BMI1) * (input$qol_BMI1 <= 0), 
              (input$qol_qol_BMI3) * (input$qol_BMI3 <= 0)) +
          (input$qol_transplant) * (input$qol_transplant <= 0) +
          (input$qol_diabeticNep) * (input$qol_diabeticNep <= 0) +
          (input$qol_dialysis) * (input$qol_dialysis <= 0) +
          min((input$qol_mve0) * (input$qol_mve0 <= 0), 
              (input$qol_mve1) * (input$qol_mve2 <= 0),
              (input$qol_mve2) * (input$qol_mve2 <= 0))
        # maximal possible QoL
        qol_max <- input$qol_baseline + 
          (3 * input$qol_age) * (input$qol_age >= 0) + 
          (2 * input$qol_age) * (input$qol_age <= 0) +
          (input$qol_sexm) * (input$qol_sexm >= 0) +
          max((input$qol_eduGCSE) * (input$qol_eduGCSE >= 0), 
              (input$qol_eduBelowSec) * (input$qol_eduBelowSec >= 0)) +
          max((input$qol_smokerbefore) * (input$qol_smokerbefore >= 0), 
              (input$qol_smokercurrently) * (input$qol_smokercurrently >= 0)) +
          max((input$qol_BMI1) * (input$qol_BMI1 >= 0), 
              (input$qol_BMI3) * (input$qol_BMI3 >= 0)) +
          (input$qol_transplant) * (input$qol_transplant >= 0) +
          (input$qol_diabeticNep) * (input$qol_diabeticNep >= 0) +
          (input$qol_dialysis) * (input$qol_dialysis >= 0) +
          max((input$qol_mve0) * (input$qol_mve0 >= 0), 
              (input$qol_mve1) * (input$qol_mve2 >= 0),
              (input$qol_mve2) * (input$qol_mve2 >= 0))
        
        return(list(
          x_lower = x,
          qol_min_check = (x >= -2 & x <= 0),
          qol_min = (qol_min >= x),
          qol_max = (qol_max <= 1)
        ))
      }
  })
  
  output$validated_qol_value <- renderUI({
    validate_qol(qol_value())
  })
  
  ### NVD
  
  nvd_value <- reactive({
    
    inFile <- input$file2
    
    if (is.null(inFile))
      return(NULL)
    
    # all values except for the ageBand
    alpha <- .validate_file_nonempty(inFile = inFile,
                                     req_colnames = req_colnames_nvd,
                                     req_format = req_format_nvd,
                                     req_values = req_values_nvd)
    
    #vars_not_validated <- c(alpha$vars_missing, alpha$vars_format_names, alpha$vars_values_names)
    #
    #ageBand_values <- TRUE
    #if (length(intersect(c("CKDStage", "ageBand",  "SEX"), alpha$vars_not_validated)) == 0) {
    #  x <- as.data.frame(t(sapply(df$y, function(x) unlist(strsplit(x, "-")))), stringsAsFactors = FALSE)
    #  x[, 1] <- as.numeric(x[, 1])
    #  x[, 2] <- as.numeric(x[, 2])
    #  colnames(x) <- c("T1", "T2")
    #  df <- cbind(df, x)
    #  check_Inf <- TRUE
    #  for (sex in c("F", "M"))
    #    for (ckdstage in c('3B', '4', '5', 'dialysis', 'transplant'))
    #      if (!(Inf %in% df$T2[df$SEX == sex & df$CKDStage == ckdstage]))
    #        check_Inf <- FALSE
    #  ageBand_values <- all(!is.na(x[, 1])) & all(!is.na(x[, 2])) & all(x[, 1] > 0) & all(x[, 2] > x[, 1]) & check_Inf = TRUE
    #}
    
    return(list(vars_missing = alpha$vars_missing, 
                vars_format = alpha$vars_format, 
                vars_values = alpha$vars_values))
  })
  
  output$validated_nvd_value <- renderUI({
    validate_file(nvd_value()$vars_missing, 
                  nvd_value()$vars_format, 
                  nvd_value()$vars_values)
  })
  
  ### decision parameters
  
  dp_screen_value <- reactive({
    if (is.null(input$dp_lifetime))
      return(list(
        fixed = TRUE, cycles = TRUE))
    if (is.null(input$slider_age)) {
      if (input$dp_lifetime == 'Simulation for a fixed number of years') {
        x <- input$dp_cycles
        return(list(
          fixed = TRUE,
          cycles = (x >= 0 & x + 65 <= 100),
          age = 65
        )) 
      } else {
        x <- input$dp_stopping_age
        return(list(
          fixed = FALSE,
          stopping_age = (x > 65 & x <= 100)
        ))
      }
    } else {
      if (input$dp_lifetime == 'Simulation for a fixed number of years') {
        x <- input$dp_cycles
        return(list(
          fixed = TRUE,
          cycles = (x >= 0 & x + input$slider_age <= 100),
          age = input$slider_age
        )) 
      } else {
        x <- input$dp_stopping_age
        return(list(
          fixed = FALSE,
          stopping_age = (x > input$slider_age & x <= 100)
        ))
      }
    }
  })
  
  dp_file_value <- reactive({
    
    inFile <- input$file1
    
    # check the file is in the correct format
    if (is.null(inFile))
      return(NULL)
    check_file <- FALSE
    df <- read.csv(inFile$datapath, stringsAsFactors=FALSE)
    if ("age" %in% colnames(df) & is.numeric(df$age) & all(df$age >= 40) & all(df$age <= 90))
      check_file <- TRUE
    if (check_file == FALSE)
      return(list(file = FALSE))
    
    # file in the correct format, proceed
    max_age <- max(df$age)
    
    if (is.null(input$dp_lifetime))
      return(list(
        file = TRUE, fixed = TRUE, cycles = max_age <= 70)) else {
          if (input$dp_lifetime == 'Simulation for a fixed number of years') {
            x <- input$dp_cycles
            return(list(
              file = TRUE,
              fixed = TRUE,
              cycles = (x >= 0 & x + max_age <= 100),
              age = max_age
            )) 
          } else {
            x <- input$dp_stopping_age
            return(list(
              file = TRUE,
              fixed = FALSE,
              stopping_age = (x > max_age & x <= 100)
            ))
          }
        }
  })
  
  output$validated_dp_value <- renderUI({
    if (input$cb_bl) 
      validate_dp_file(dp_file_value()) else
        validate_dp_screen(dp_screen_value())
  })
  
  ### validation results page
  
  output$validated_all <- renderUI({
    if (is.null(input$PSA) | input$PSA == 'No (deterministic analysis)') {
      if (input$cb_bl)  
        validate_file(id_file_value()$vars_missing, id_file_value()$vars_format, id_file_value()$vars_values) else 
          validate_id_screen(id_screen_value())
      validate_trt(trt_value())
      validate_cost(cost_value())
      validate_qol(qol_value())
      if (input$cb_nvd) 
        validate_file(nvd_value()$vars_missing, nvd_value()$vars_format, nvd_value()$vars_values)
      if (input$cb_bl)
        validate_dp_file(dp_file_value()) else
          validate_dp_screen(dp_screen_value())
    } else {
      if (input$cb_bl)  
        validate_file(id_file_value()$vars_missing, id_file_value()$vars_format, id_file_value()$vars_values) else 
          validate_id_screen(id_screen_value())
      validate_trt_PSA(trt_value_PSA())
      validate_cost_PSA(cost_value_PSA())
      validate_qol(qol_value())
      if (input$cb_nvd) 
        validate_file(nvd_value()$vars_missing, nvd_value()$vars_format, nvd_value()$vars_values)
      if (input$cb_bl)
        validate_dp_file(dp_file_value()) else
          validate_dp_screen(dp_screen_value())
    }
    actionButton("action1", "Run analyses")
  })
  
  ################################################################################
  ### Reset values to default
  ################################################################################
  
  output$reset_id <- renderUI({
    times <- input$reset_input_id
    list(
      br(),
      h4(em("Demographic and socio-economic characteristics")),
      fluidRow(column(4, 
                      numericInput("slider_age", label = "Age (years)", 
                                   min = 40, max = 90, value = 65)),
               column(4, 
                      selectInput("sel_sex", "Gender", 
                                  choices = c("Female", "Male"), selected = "Female")),
               column(4, selectInput("sel_ethn", "Ethnicity", 
                                     choices = c("White", "Asian, lives in China", "Asian, lives outside China",
                                                 "Black", "Other"), selected = "White"))),
      fluidRow(column(4, 
                      selectInput("sel_ed", "Highest educational attainment", 
                                  choices = c("Any post-secondary education", "Completed secondary education", "Below secondary education"),
                                  selected = "Any post-secondary education")),
                 #column(4, 
              #        selectInput("sel_childdep", "Child dependants", 
              #                    choices = c("No", "Yes"), selected = "No")),
               column(4, 
                      selectInput("sel_adultdep", "Adult dependants", 
                                  choices = c("No", "Yes"), selected = "No")),
              column(4, 
                     selectInput("sel_smok", "Smoking status", 
                                 choices = c("Never smoked", "Ex-smoker", "Current smoker"), 
                                 selected = "Never smoked"))),
      fluidRow(
               column(4, 
                      selectInput("sel_alc", "Alcohol drinker", 
                                  choices = c("No", "Yes"), selected = "No")),
               column(4, 
                      selectInput("sel_bmi", "Body mass index", 
                                  choices = c("<25 kg/m\u00B2", "25-29 kg/m\u00B2", "\u226530 kg/m\u00B2"), 
                                  selected = "25-29 kg/m\u00B2"))),
      h4(em("Clinical factors")),
      fluidRow(column(4, 
                      selectInput("sel_dbp", "Diastolic blood pressure", 
                                  choices = c("<75 mmHg", "75-84 mmHg","\u226585 mmHg"), 
                                  selected = "75-84 mmHg")),
               column(4, selectInput("sel_sbp", "Systolic blood pressure", 
                                     choices = c("<130 mmHg", "130-149 mmHg","\u2265150 mmHg"), 
                                     selected = "130-149 mmHg")),
               column(4, 
                      selectInput("sel_hdl", "HDL cholesterol", 
                                  choices = c("<0.9 mmol/L", "0.9-1.1 mmol/L","\u22651.2 mmol/L"), 
                                  selected = "0.9-1.1 mmol/L"))),
      fluidRow(column(4, 
                      selectInput("sel_alb", "Albumin", 
                                  choices = c("<3.9 g/dL","3.9-4.1 g/dL", "\u22654.2 g/dL"), 
                                  selected = "3.9-4.1 g/dL")),
               column(4, selectInput("sel_haem", "Haemoglobin", 
                                     choices = c("<11.6 g/dL","11.6-12.9 g/dL", "\u226513.0 g/dL"), 
                                     selected = "11.6-12.9 g/dL")),
               column(4, selectInput("sel_phos", "Phosphate", 
                                     choices = c("<1.2 mmol/L","1.2-1.4 mmol/L", "\u22651.5 mmol/L"), 
                                     selected = "1.2-1.4 mmol/L"))),
      conditionalPanel(
        condition = "input.sel_ckd != 'Dialysis' & input.sel_ckd != 'Transplant'",
        fluidRow(column(4, 
                        selectInput("sel_acr", "Urinary albumin:creatinine ratio", 
                                    choices = c("<30 mg/g","30-300 mg/g", ">300 mg/g"), 
                                    selected = "30-300 mg/g"))
        )),
      h4(em("Disease history")),
      fluidRow(column(4, selectInput("sel_vd", "Latest cardiovascular event", 
                                     choices = c("None", "Major atherosclerotic event in the last year",
                                                 "Major atherosclerotic event 1-2 years ago",
                                                 "Major atherosclerotic event >2 years ago",
                                                 "No MAE, but haemorrhagic stroke in the last year",
                                                 "No MAE, but haemorrhagic stroke 1-2 years ago",
                                                 "No MAE, but haemorrhagic stroke >2 years ago",
                                                 "Another cardiovascular event"), 
                                     selected = "None")),
               column(4, selectInput("sel_diab", "Diabetes", 
                                     choices = c("No", "Yes"), selected = "No"))
      ),
      fluidRow(column(4, selectInput("sel_ckd", "CKD stage", 
                                     choices = c("CKD 3B", "CKD 4", "CKD 5, not RRT", "Dialysis", "Transplant"), 
                                     selected = "CKD 3B")),
               column(4, numericInput("sel_ckddur", "CKD duration (years)", 
                                      value = 10, min = 0, max = 60)),
               column(4, selectInput("sel_rdiag", 
                                     "Renal diagnosis", 
                                     choices = c("Diabetic nephropathy", "Cystic kidney disease", "Other known or unknown cause"), 
                                     selected = "Other known or unknown cause"))),
      conditionalPanel(
        condition = "input.sel_ckd == 'Dialysis' || input.sel_ckd == 'Transplant'",
        fluidRow(column(4, numericInput("sel_esrddur", 
                                        "RRT duration (years)", value = 5, min = 0, max = 60)),
                 column(4, selectInput("sel_tx", "Previous (failed) transplant", 
                                       choices = c("No","Yes"), selected = "No"))
        )))
  })
  
  output$reset_trt <- renderUI({
    
    type <- input$anal_type
    
    if (is.null(type) | type == "Long-term projections") {
      RR <- 1
      compl <- 0
      price <- 0
    } else {
      RR <- 0.9
      compl <- 100
      price <- 1
    }
    
    times <- input$reset_input_trt
    list(
      h3("Treatment effects"),
      h4(em("Cardiovascular death")),
      numericInput("RR_VD", label = "Hazard ratio",  
                   value = RR, min = 0, step = 0.01),
      h4(em("Cardiovascular death or non-fatal major atherosclerotic event")),
      numericInput("RR_NFMAE", label = "Hazard ratio",  
                   value = RR, min = 0, step = 0.01),
      h4(em("Cardiovascular death or non-fatal major vascular event")),
      numericInput("RR_NFMVE", label = "Hazard ratio",  
                   value = RR, min = 0, step = 0.01),
      h3("Compliance (%)"),
      numericInput("compl1", 
                   label = NULL, 
                   value = compl, min = 0, max = 100),
      h3("Daily treatment cost (full use)"),
      numericInput("price1", NULL, 
                   value = price, min = 0, step=0.01)
    )
  })
  
  output$reset_trt_PSA <- renderUI({
    
    type <- input$anal_type
    
    if (is.null(type) | type == "Long-term projections") {
      RR <- 1
      RR_l <- 0.9
      RR_r <- 1.1
      compl <- 0
      price <- 0
    } else {
      RR <- 0.9
      RR_l <- 0.8
      RR_r <- 1
      compl <- 100
      price <- 1
    }
    
    times <- input$reset_input_trt_PSA
    list(
      h3("Treatment effects"),
      p("Treatment effects for the probabilistic sensitivity analyses 
        are sampled from log-normal distributions using the correlation matrix from the SHARP study. Enter the estimates for the hazard ratios 
        together with the 95% confidence interval (CI) on the exponential scale."),
      h4(em("Cardiovascular death")),
      fluidRow(column(4, 
                      numericInput("RR_VD_2", label = "Hazard ratio",
                                   value = RR, min = 0, max = 1.1, step = 0.01)),
               column(4, 
                      numericInput("RR_VD_l", label = "Lower 95% CI",
                                   value = RR_l, min = 0, max = 1.1, step = 0.01)),
               column(4, numericInput("RR_VD_r", label = "Upper 95% CI",
                                      value = RR_r, min = 0, max = 1.1, step = 0.01))),
      h4(em("Cardiovascular death or non-fatal major atherosclerotic event")),
      fluidRow(column(4, 
                      numericInput("RR_NFMAE_2", label = "Hazard ratio",
                                   value = RR, min = 0, max = 1.1, step = 0.01)),
               column(4, 
                      numericInput("RR_NFMAE_l", label = "Lower 95% CI",
                                   value = RR_l, min = 0, max = 1.1, step = 0.01)),
               column(4, numericInput("RR_NFMAE_r", label = "Upper 95% CI",
                                      value = RR_r, min = 0, max = 1.1, step = 0.01))),
      h4(em("Cardiovascular death or non-fatal major vascular event")),
      fluidRow(column(4, 
                      numericInput("RR_NFMVE_2", label = "Hazard ratio",
                                   value = RR, min = 0, max = 1.1, step = 0.01)),
               column(4, 
                      numericInput("RR_NFMVE_l", label = "Lower 95% CI",
                                   value = RR_l, min = 0, max = 1.1, step = 0.01)),
               column(4, numericInput("RR_NFMVE_r", label = "Upper 95% CI",
                                      value = RR_r, min = 0, max = 1.1, step = 0.01))),
      h4("Compliance (%)"),
      numericInput("compl1_2", 
                   label = NULL, 
                   value = compl, min = 0, max = 100),
      h3("Daily treatment cost (full use)"),
      numericInput("price1_2", NULL, 
                   value = price, min = 0, step=0.01)
      )
  })
  
  output$reset_cost <- renderUI({
    times <- input$reset_input_cost
    list(
      h3("Annual cost of CKD"),
      h4(em("CKD stage 3B")),
      numericInput("cost_ckd13b", "mean estimate", 
                   value = 427, min = 0, step = 10),
      h4(em("CKD stage 4")),
      numericInput("cost_ckd4", "mean estimate", 
                   value = 417, min = 0, step = 10),
      h4(em("CKD stage 5")),
      numericInput("cost_ckd5", "mean estimate", 
                   value = 556, min = 0, step = 10),
      h4(em("On dialysis, for year of dialysis initiation")),
      numericInput("cost_dial0", "mean estimate", 
                   value = 20112, min = 0, step = 10),
      h4(em("On dialysis, not for year of dialysis initiation")),
      numericInput("cost_dial1", "mean estimate", 
                   value = 24709, min = 0, step = 10),
      h4(em("With kidney transplant, for year of kidney transplant")),
      numericInput("cost_tx0", "mean estimate", 
                   value = 26061, min = 0, step = 10),
      h4(em("With kidney transplant, not for year of kidney transplant")),
      numericInput("cost_tx1", "mean estimate", 
                   value = 1216, min = 0, step = 10),
      h3("Additional annual costs"),
      h4(em("Diabetes")),
      numericInput("cost_diab", "mean estimate", 
                   value = 181, min = 0, step = 10),
      h4(em("Major non-fatal vascular event this year, not on dialysis")),
      numericInput("cost_cvd5_nd", "mean estimate", 
                   value = 4607, min = 0, step = 10),
      h4(em("Major non-fatal vascular event this year, on dialysis")),
      numericInput("cost_cvd5_d", "mean estimate", 
                   value = 6497, min = 0, step = 10),
      h4(em("Major non-fatal vascular event in the previous year")),
      numericInput("cost_cvd6", "mean estimate", 
                   value = 782, min = 0, step = 10),
      h4(em("Other cardiovascular event")),
      numericInput("cost_cvd478", "mean estimate", 
                   value = 182, min = 0, step = 10),
      h4(em("Vascular death")),
      numericInput("cost_vd", "mean estimate", 
                   value = 1204, min = 0, step = 10),
      h4(em("Non-vascular death")),
      numericInput("cost_nvd", "mean estimate", 
                   value = 1474, min = 0, step = 10))
    #h3("Other costs"),
    #h4(em("Reduction of dialysis costs in the year of death")),
    #numericInput("cost_death_dial", "mean estimate", 
    #             value = 9324, min = 0, step = 10))    
  })
  
  output$reset_cost_PSA <- renderUI({
    times <- input$reset_input_cost_PSA
    list(
      br(),
      p("The default costs for the probabilistic sensitivity analyses are derived 
        from the SHARP data using the bootstrap method. 
        To provide alternative costs, enter the means
        and the standard errors below, and the costs will be sampled from gamma distributions.
        The displayed values are based on SHARP data and UK 2014 prices [1]."),
      h3("Annual cost of CKD"),
      h4(em("CKD stage 3B")),
      fluidRow(column(6, 
                      numericInput("cost_ckd13b_2", label = "mean estimate ",
                                   value = 427, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_ckd13b_sd", label = "standard error",
                                   value = 32, min = 1, step = 10))),
      h4(em("CKD stage 4")),
      fluidRow(column(6, 
                      numericInput("cost_ckd4_2", label = "mean estimate",
                                   value = 417, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_ckd4_sd", label = "standard error",
                                   value = 27, min = 1, step = 10))),
      h4(em("CKD stage 5")),
      fluidRow(column(6, 
                      numericInput("cost_ckd5_2", label = "mean estimate",
                                   value = 556, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_ckd5_sd", label = "standard error",
                                   value = 41, min = 1, step = 10))),
      h4(em("On dialysis, for year of dialysis initiation")),
      fluidRow(column(6, 
                      numericInput("cost_dial0_2", label = "mean estimate",
                                   value = 20112, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_dial0_sd", label = "standard error",
                                   value = 198, min = 1, step = 10))),
      h4(em("On dialysis, not for year of dialysis initiation")),
      fluidRow(column(6, 
                      numericInput("cost_dial1_2", label = "mean estimate",
                                   value = 24709, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_dial1_sd", label = "standard error",
                                   value = 51, min = 1, step = 10))),
      h4(em("With kidney transplant, for year of kidney transplant")),
      fluidRow(column(6, 
                      numericInput("cost_tx0_2", label = "mean estimate",
                                   value = 26061, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_tx0_sd", label = "standard error",
                                   value = 311, min = 1, step = 10))),
      h4(em("With kidney transplant, not for year of kidney transplant")),
      fluidRow(column(6, 
                      numericInput("cost_tx1_2", label = "mean estimate",
                                   value = 1216, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_tx1_sd", label = "standard error",
                                   value = 92, min = 1, step = 10))),
      h3("Additional annual costs"),
      h4(em("Diabetes")),
      fluidRow(column(6, 
                      numericInput("cost_diab_2", label = "mean estimate",
                                   value = 181, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_diab_sd", label = "standard error",
                                   value = 63, min = 1, step = 10))),
      h4(em("Major non-fatal vascular event this year, not on dialysis")),
      fluidRow(column(6, 
                      numericInput("cost_cvd5_nd_2", label = "mean estimate",
                                   value = 4607, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_cvd5_nd_sd", label = "standard error",
                                   value = 287, min = 1, step = 10))),
      h4(em("Major non-fatal vascular event this year, on dialysis")),
      fluidRow(column(6, 
                      numericInput("cost_cvd5_d_2", label = "mean estimate",
                                   value = 6497, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_cvd5_d_sd", label = "standard error",
                                   value = 284, min = 1, step = 10))),
      h4(em("Major non-fatal vascular event in the previous year")),
      fluidRow(column(6, 
                      numericInput("cost_cvd6_2", label = "mean estimate",
                                   value = 782, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_cvd6_sd", label = "standard error",
                                   value = 209, min = 1, step = 10))),
      h4(em("Other prior vascular disease")),
      fluidRow(column(6, 
                      numericInput("cost_cvd478_2", label = "mean estimate",
                                   value = 182, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_cvd478_sd", label = "standard error",
                                   value = 286, min = 1, step = 10))),
      h4(em("Vascular death")),
      fluidRow(column(6, 
                      numericInput("cost_vd_2", label = "mean estimate",
                                   value = 1204, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_vd_sd", label = "standard error",
                                   value = 361, min = 1, step = 10))),
      h4(em("Non-vascular death")),
      fluidRow(column(6, 
                      numericInput("cost_nvd_2", label = "mean estimate",
                                   value = 1474, min = 0, step = 10)),
               column(6, 
                      numericInput("cost_nvd_sd", label = "standard error",
                                   value = 201, min = 1, step = 10))),
      #h3("Other costs"),
      #h4(em("Reduction of dialysis costs in the year of death")),
      #fluidRow(column(6, 
      #                numericInput("cost_death_dial_2", label = "mean estimate",
      #                             value = 9324, min = 0, step = 10)),
      #         column(6, 
      #                numericInput("cost_death_dial_sd", label = "standard error",
      #                             value = 340, min = 1, step = 10))),
      br(),
      br(),
      p("[1] Kent S, Schlackow I, Lozano-K\u00FChne J, Reith C, Emberson J, Haynes R, 
        Gray A, Cass A, Baigent C, Landray MJ, Herrington W, Mihaylova B, 
        on behalf of the SHARP Collaborative Group.,", 
        em("What is the impact of chronic kidney disease stage and cardiovascular 
           disease on the annual cost of hospital care in moderate-to-severe kidney disease?"), 
        "BMC Nephrology 2015; 16:65." ))})
  
  output$reset_QoL <- renderUI({
    times <- input$reset_input_QoL
    list(
      h3("Baseline QoL"),
      fluidRow(column(6, 
                      numericInput("qol_baseline", NULL, 
                                   value = 0.860, min = 0, step = 0.01))),
      h3("Additional effects"),
      h4(em("Demographic and socio-economic characteristics")),
      fluidRow(column(6, 
                      numericInput("qol_age", "Age (per 10 years)", 
                                   value = -0.048, min = 0, step = 0.01)),
               column(6, 
                      numericInput("qol_sexm", "Male", 
                                   value = 0.059, min = 0, step = 0.01))),
      fluidRow(column(6, 
                      numericInput("qol_eduGCSE", "Completed secondary education", 
                                   value = -0.017, min = 0, step = 0.01)),
               column(6, 
                      numericInput("qol_eduBelowSec", "Below secondary education", 
                                   value = -0.036, min = 0, step =0.01))),
      fluidRow(column(6, 
                      numericInput("qol_smokerbefore", "Ex-smoker", 
                                   value = -0.009, min = 0, step = 0.01)),
               column(6, 
                      numericInput("qol_smokercurrently", "Current smoker", 
                                   value = -0.037, min = 0, step = 0.01))),
      fluidRow(column(6, 
                      numericInput("qol_BMI1", "BMI <25 kg/m\u00B2", 
                                   value = 0.011, min = 0, step = 0.01)),
               column(6, 
                      numericInput("qol_BMI3", "BMI \u226530 kg/m\u00B2", 
                                   value = -0.043, min = 0, step = 0.01))),
      h4(em("Disease history")),
      fluidRow(column(6, 
                      numericInput("qol_mve0", "Major vascular event in the last year", 
                                   value = -0.173, min = 0, step = 0.01)),
               column(6, 
                      numericInput("qol_mve1", "Major vascular event more than a year ago, but not in the last year", 
                                   value = -0.103, min = 0, step = 0.01))),
      fluidRow(column(6, 
                      numericInput("qol_mve2", "No major vascular events, but another cardiovascular event", 
                                   value = -0.071, min = 0, step = 0.01)),
               column(6, 
                      numericInput("qol_dialysis", "Dialysis", 
                                   value = -0.056, min = 0, step = 0.01))),
      fluidRow(column(6, 
                      numericInput("qol_diabeticNep", "Diabetic nephropathy", 
                                   value = -0.059, min = 0, step = 0.01)),
               column(6, 
                      numericInput("qol_transplant", "Previous (failed) kidney transplant", 
                                   value = -0.071, min = 0, step = 0.01))),
      h4("Validation check: minimum admissible utility", style = "color:rgb(0, 102, 213)"),
      fluidRow(column(6, 
                      numericInput("qol_min", NULL, 
                                   value = -0.594, min = -2, max = 0, step = 0.01))))
  }) 
  
  output$reset_dp <- renderUI({
    times <- input$reset_input_dp
    list(
      sliderInput("dr_cost", label = "Discount rate: costs", 
                  min = 0, max = 10, value = 3.5, step = .5),
      br(),
      sliderInput("dr_ev", label = "Discount rate: health outcomes", 
                  min = 0, max = 10, value = 3.5, step = .5),
      br(),
      selectInput("dp_lifetime", "Duration of model execution",  choices = c("Lifetime simulation", "Simulation for a fixed number of years"), selected = "Simulation for a fixed number of years"),
      br(),
      conditionalPanel(          
        condition = "input.dp_lifetime == 'Simulation for a fixed number of years'",    
        numericInput("dp_cycles", label = "Number of years of analysis", value = 30, min = 1)),
      conditionalPanel(          
        condition = "input.dp_lifetime == 'Lifetime simulation'",    
        numericInput("dp_stopping_age", label = "Stopping age", value = 95, min = 41))
    )
  })
  
  newdata <- eventReactive(input$action1, {
    
    ### this is dependency on the button action1, ie calculations are only executed if the button is clicked
    
    ################################################################################
    ### Type of analysis
    ################################################################################

    LE <- (is.null(input$anal_type) | input$anal_type == "Long-term projections")
    PSA <- (input$PSA == "Yes (probabilistic analysis)")
    
    ################################################################################
    ### Patient characteristics
    ################################################################################
    
    if (is.null(input$slider_age)) {
      # the demographics screen was not clicked; use default values
      stageRand_2 <- stages[data$CKDStage + 1]
      CVRand <- CV[data$CVD + 1]
    } else {
      if (input$cb_bl == TRUE) {
        
        ### if tickbox ticked, read off from the provided file  
        
        inFile <- input$file1
        data <- read.csv(inFile$datapath, check.names = F)
        
      } else {
        
        ### the demographics screen was clicked and the information is read off from the screen
        
        data <- data.frame(id = 1, 
                           age = input$slider_age)
        
        data$sex <- as.numeric(input$sel_sex == "Male")
        
        ethnicity <- input$sel_ethn
        data$ethnicity <- if (ethnicity == "White") 0 else
          if (ethnicity == "Asian, lives in China") 1 else
            if (ethnicity == "Asian, lives outside China") 2 else
              if (ethnicity == "Black") 3 else
                4
        
        education <- input$sel_ed
        data$education <- if (education == "Any post-secondary education") 0 else
          if (education == "Completed secondary education") 1 else
            2
        
        #data$childDep <- as.numeric(input$sel_childdep == "Yes")
        
        data$adultDep <- as.numeric(input$sel_adultdep == "Yes")
        
        smoker <- input$sel_smok     
        data$smoker <- if (smoker == "Never smoked") 0 else
          if (smoker == "Ex-smoker") 1 else
            2
        
        data$currentAlc <- as.numeric(input$sel_alc == "Yes")
        
        BMI_quant <- input$sel_bmi
        data$BMI_quant <- if (BMI_quant == "25-29 kg/m\u00B2") 0 else
          if (BMI_quant == "<25 kg/m\u00B2") 1 else
            2
        
        DBP_quant <- input$sel_dbp
        data$DBP_quant <- if (DBP_quant == "75-84 mmHg") 0 else
          if (DBP_quant == "<75 mmHg") 1 else
            2
        
        SBP_quant <- input$sel_sbp
        data$SBP_quant <- if (SBP_quant == "130-149 mmHg") 0 else
          if (SBP_quant == "<130 mmHg") 1 else
            2
        
        CholHDL_quant <- input$sel_hdl
        data$CholHDL_quant <- if (CholHDL_quant == "0.9-1.1 mmol/L") 0 else
          if (CholHDL_quant == "<0.9 mmol/L") 1 else
            2
        
        Albumin_quant <- input$sel_alb
        data$Albumin_quant <- if (Albumin_quant == "3.9-4.1 g/dL") 0 else
          if (Albumin_quant == "<3.9 g/dL") 1 else
            2
        
        Hemoglobin_quant <- input$sel_haem
        data$Hemoglobin_quant <- if (Hemoglobin_quant == "11.6-12.9 g/dL") 0 else
          if (Hemoglobin_quant == "<11.6 g/dL") 1 else
            2
        
        Phosphate_quant <- input$sel_phos
        data$Phosphate_quant <- if (Phosphate_quant == "1.2-1.4 mmol/L") 0 else
          if (Phosphate_quant == "<1.2 mmol/L") 1 else
            2
        
        ACR_quant <- input$sel_acr
        data$ACR_quant <- if (input$sel_ckd %in% c("Dialysis", "Transplant")) 3 else 
          if (ACR_quant == "30-300 mg/g") 0 else
            if (ACR_quant == "<30 mg/g") 1 else
              2
        
        CVD <- input$sel_vd
        data$CVD <- if (CVD == "None") 0 else
          if (CVD == "Major atherosclerotic event in the last year") 1 else
            if (CVD == "Major atherosclerotic event 1-2 years ago") 2 else
              if (CVD == "Major atherosclerotic event >2 years ago") 3 else
                if (CVD == "No MAE, but haemorrhagic stroke in the last year") 4 else
                  if (CVD == "No MAE, but haemorrhagic stroke 1-2 years ago") 5 else
                    if (CVD == "No MAE, but haemorrhagic stroke >2 years ago") 6 else
                      7
        
        data$DM <- as.numeric(input$sel_diab == "Yes")
        
        CKDStage <- input$sel_ckd
        data$CKDStage <- if (CKDStage == "CKD 3B") 0 else
          if (CKDStage == "CKD 4") 1 else
            if (CKDStage == "CKD 5, not RRT") 2 else
              if (CKDStage == "Dialysis") 3 else
                if (CKDStage == "Transplant") 4
        
        data$CKDDuration <- input$sel_ckddur
        
        renalDiagnosis <- input$sel_rdiag
        data$renalDiagnosis <- if (renalDiagnosis == "Other known or unknown cause") 2 else
          if (renalDiagnosis == "Diabetic nephropathy") 0 else
            1
        
        data$RRTDuration <- if (CKDStage %in% c("CKD 3B", "CKD 4", "CKD 5, not RRT")) 0 else
          input$sel_esrddur
        
        data$TX <- if (CKDStage %in% c("Dialysis", "Transplant"))
          as.numeric(input$sel_tx == "Yes") else
            0
        
      }
      
      ### put the data into the format suitable for the code
      
      # time-updated variables
      patT <- as.matrix(data.frame(id = data$id, 
                                   age_T = data$age, age_T2 = data$age + 1,
                                   CKDDuration_T = data$CKDDuration, 
                                   ESRDDuration_T = data$RRTDuration))
      
      # not time-updated variables
      
      pat0 <- data
      
      colnames(pat0)[which(colnames(pat0) == "sex")] <- "SEX"
      
      pat0$SEX <- factor(pat0$SEX, 
                         levels = c(0, 1),
                         labels = c("F", "M"))
      pat0$ethnicity <- factor(pat0$ethnicity,
                               levels = c(0, 1, 2, 3, 4),
                               labels = c("White", "Asian: China", "Asian: other", "Black", "Other"))
      pat0$education <- factor(pat0$education, 
                               levels = c(0, 1, 2),
                               labels = c("A-levels and above", "GCSE/vocational", "below secondary"))
      #pat0$childDep <- factor(pat0$childDep, 
      #                        levels = c(0, 1),
      #                        labels = c("0", ">0"))
      pat0$adultDep <- factor(pat0$adultDep, 
                              levels = c(0, 1),
                              labels = c("1", ">1"))
      pat0$smoker <- factor(pat0$smoker, 
                            levels = c(0, 1, 2),
                            labels = c("never", "before", "currently"))
      pat0$currentAlc <- factor(pat0$currentAlc, 
                                levels = c(0, 1),
                                labels = c(0, 1))
      pat0$BMI_quant <- factor(pat0$BMI_quant, 
                               levels = c(0, 1, 2),
                               labels = c("T2", "T1", "T3"))
      pat0$DBP_quant <- factor(pat0$DBP_quant, 
                               levels = c(0, 1, 2), 
                               labels = c("T2", "T1", "T3"))
      pat0$SBP_quant <- factor(pat0$SBP_quant, 
                               levels = c(0, 1, 2),
                               labels = c("T2", "T1", "T3"))
      pat0$CholHDL_quant <- factor(pat0$CholHDL_quant, 
                                   levels = c(0, 1, 2),
                                   labels = c("T2", "T1", "T3"))
      pat0$Albumin_quant <- factor(pat0$Albumin_quant, 
                                   levels = c(0, 1, 2),
                                   labels = c("T2", "T1", "T3"))
      pat0$Hemoglobin_quant <- factor(pat0$Hemoglobin_quant, 
                                      levels = c(0, 1, 2),
                                      labels = c("T2", "T1", "T3"))
      pat0$Phosphate_quant <- factor(pat0$Phosphate_quant, 
                                     levels = c(0, 1, 2),
                                     labels = c("T2", "T1", "T3"))
      pat0$ACR_quant <- factor(pat0$ACR_quant, 
                               levels = c(0, 1, 2, 3),
                               labels = c("T2", "T1", "T3", "ESRD"))
      pat0$DM <- factor(pat0$DM, 
                        levels = c(0, 1),
                        labels = c(0, 1))
      pat0$TX <- factor(pat0$TX, 
                        levels = c(0, 1),
                        labels = c(0, 1))
      
      pat0$renalDiagnosis <- factor(pat0$renalDiagnosis, 
                                    levels = c(2, 0, 1),
                                    labels = c("Other (known or unknown cause)", "Diabetic nephropathy", "Cystic kidney disease"))
      
      vars <- c(
        "SEX",
        "ethnicity", "smoker", "currentAlc", "adultDep", "education", "BMI_quant", 
        "DM", "TX", "DBP_quant", "SBP_quant", "Albumin_quant", "Hemoglobin_quant", "Phosphate_quant", 
        "ACR_quant", "CholHDL_quant", "renalDiagnosis")
      
      pat0 <- model.matrix(
        as.formula(paste("~", paste(vars, collapse = "+"), sep = "")), data = pat0)
      
      pat0 <- cbind(pat0, id = data$id)
      
      stageRand_2 <- stages[data$CKDStage + 1]
      CVRand <- CV[data$CVD + 1]
      
    }
    
    stageRand_2 <- as.character(stageRand_2)
    stageRand_2_num <- mapply(.get_N_CKD_0, stageRand_2 = as.character(stageRand_2), ESRDDuration_T = patT[, "ESRDDuration_T"])
    CVRand_num <- mapply(.get_N_CV_0, CVRand = CVRand)
    BVD <- (CVRand == "Another cardiovascular event")
    
    ################################################################################
    ### Treatment paramenters 
    ################################################################################
    
    compl_C <- 0
    Tx_price_C <- 0
    
    if (input$PSA == "No (deterministic analysis)") {
      if (is.null(input$RR_VD)) {
        
        # default values
        
        if (LE) {
        # compliance
        compl_T <- 0
        # intervention price
        Tx_price_T <- 0
        coeffs_default$VD["complAll"] <- log(1)
        coeffs_default$NFMAE_or_VD2[[1]]["complAll"] <- log(1)
        coeffs_default$NFMVE_or_VD2[[1]]["complAll"] <- log(1)
        } else {
          # compliance
          compl_T <- 1
          # intervention price
          Tx_price_T <- 1
          coeffs_default$VD["complAll"] <- log(0.9)
          coeffs_default$NFMAE_or_VD2[[1]]["complAll"] <- log(0.9)
          coeffs_default$NFMVE_or_VD2[[1]]["complAll"] <- log(0.9)
        }
      } else {
        # compliance
        compl_T <- input$compl1/100
        # intervention price
        Tx_price_T <- input$price1
        # det Tx effects
        coeffs_default$VD["complAll"] <- log(input$RR_VD)
        coeffs_default$NFMAE_or_VD2[[1]]["complAll"] <- log(input$RR_NFMAE)
        coeffs_default$NFMVE_or_VD2[[1]]["complAll"] <- log(input$RR_NFMVE)
      }
    } else if (input$PSA == "Yes (probabilistic analysis)") {
      if (is.null(input$RR_VD_2)) {
        
        # default values
        if (LE) {
        # compliance
        compl_T <- 0
        # intervention price
        Tx_price_T <- 0
        # values for deterministic analysis
        coeffs_default$VD["complAll"] <- log(1)
        coeffs_default$NFMAE_or_VD2[[1]]["complAll"] <- log(1)
        coeffs_default$NFMVE_or_VD2[[1]]["complAll"] <- log(1)
        # values for the PSA
        vd <- log(rlnorm(1000, meanlog = log(1), sdlog = (log(1.1) - log(0.9)) / 3.92))
        mae <- log(rlnorm(1000, meanlog = log(1), sdlog = (log(1.1) - log(0.9)) / 3.92))
        mve <- log(rlnorm(1000, meanlog = log(1), sdlog = (log(1.1) - log(0.9)) / 3.92))
        # create variables with given correlation
        coeffs_PSA_VD[, "complAll"] <- vd
        coeffs_PSA_NFMAEorVD[, "complAll"] <- 0.622585243 * vd + 0.782551989 * mae
        coeffs_PSA_NFMVEorVD[, "complAll"] <- 0.616308462 * vd + 0.774208435 * mae + 0.144101282 * mve
        } else {
          # compliance
          compl_T <- 1
          # intervention price
          Tx_price_T <- 1
          # values for deterministic analysis
          coeffs_default$VD["complAll"] <- log(0.9)
          coeffs_default$NFMAE_or_VD2[[1]]["complAll"] <- log(0.9)
          coeffs_default$NFMVE_or_VD2[[1]]["complAll"] <- log(0.9)
          # values for the PSA
          vd <- log(rlnorm(1000, meanlog = log(0.9), sdlog = (log(1) - log(0.8)) / 3.92))
          mae <- log(rlnorm(1000, meanlog = log(0.9), sdlog = (log(1) - log(0.8)) / 3.92))
          mve <- log(rlnorm(1000, meanlog = log(0.9), sdlog = (log(1) - log(0.8)) / 3.92))
          # create variables with given correlation
          coeffs_PSA_VD[, "complAll"] <- vd
          coeffs_PSA_NFMAEorVD[, "complAll"] <- 0.622585243 * vd + 0.782551989 * mae
          coeffs_PSA_NFMVEorVD[, "complAll"] <- 0.616308462 * vd + 0.774208435 * mae + 0.144101282 * mve
        }
      } else {
        # compliance
        compl_T <- input$compl1_2/100
        # intervention price
        Tx_price_T <- input$price1_2
        # values for the deterministic analysis
        coeffs_default$VD["complAll"] <- log(input$RR_VD_2)
        coeffs_default$NFMAE_or_VD2[[1]]["complAll"] <- log(input$RR_NFMAE_2)
        coeffs_default$NFMVE_or_VD2[[1]]["complAll"] <- log(input$RR_NFMVE_2)
        # values for the PSA
        # sample Tx effect
        vd <- log(rlnorm(1000, meanlog = log(input$RR_VD_2), sdlog = (log(input$RR_VD_r) - log(input$RR_VD_l)) / 3.92))
        mae <- log(rlnorm(1000, meanlog = log(input$RR_NFMAE_2), sdlog = (log(input$RR_NFMAE_r) - log(input$RR_NFMAE_l)) / 3.92))
        mve <- log(rlnorm(1000, meanlog = log(input$RR_NFMVE_2), sdlog = (log(input$RR_NFMVE_r) - log(input$RR_NFMVE_l)) / 3.92))
        # create variables with given correlation
        coeffs_PSA_VD[, "complAll"] <- vd
        coeffs_PSA_NFMAEorVD[, "complAll"] <- 0.622585243 * vd + 0.782551989 * mae
        coeffs_PSA_NFMVEorVD[, "complAll"] <- 0.616308462 * vd + 0.774208435 * mae + 0.144101282 * mve
      }
    } 
    
    
    ################################################################################
    ### Costs
    ################################################################################
    
    if ((is.null(input$PSA) | input$PSA == "No (deterministic analysis)") & !(is.null(input$cost_ckd13b))) {
      # det costs
      cost_coef <- coeffs_default$cost2
      cost_coef["(Intercept)"] <- input$cost_ckd13b
      cost_coef["rs2"] <- input$cost_ckd4 - input$cost_ckd13b
      cost_coef["rs3"] <- input$cost_ckd5 - input$cost_ckd13b
      cost_coef["rs4"] <- input$cost_tx0 - input$cost_ckd13b
      cost_coef["rs5"] <- input$cost_tx1 - input$cost_ckd13b
      cost_coef["rs6"] <- input$cost_dial0 - input$cost_ckd13b
      cost_coef["rs7"] <- input$cost_dial1 - input$cost_ckd13b
      cost_coef["cvd2"] <- input$cost_vd
      cost_coef["cvd3"] <- input$cost_nvd
      cost_coef["cvd5b"] <- input$cost_cvd5_nd
      cost_coef["cvd6b"] <- input$cost_cvd6
      cost_coef["cvd478b"] <- input$cost_cvd478
      cost_coef["dial.cvd5bTRUE"] <- input$cost_cvd5_d - input$cost_cvd5_nd
      #cost_coef["death_dialTRUE"] <- -input$cost_death_dial
      cost_coef["death_dialTRUE"] <- -(cost_coef["rs6"] + cost_coef["rs7"]) * 0.42 / 2

      cost_coef["DM1"] <- input$cost_diab
      coeffs_default$cost2 <- cost_coef
    } else 
      if (input$PSA == "Yes (probabilistic analysis)" & !(is.null(input$cost_ckd13b_2))) {
        # det costs
        cost_coef <- coeffs_default$cost2
        cost_coef["(Intercept)"] <- input$cost_ckd13b_2
        cost_coef["rs2"] <- input$cost_ckd4_2 - input$cost_ckd13b_2
        cost_coef["rs3"] <- input$cost_ckd5_2 - input$cost_ckd13b_2
        cost_coef["rs4"] <- input$cost_tx0_2 - input$cost_ckd13b_2
        cost_coef["rs5"] <- input$cost_tx1_2 - input$cost_ckd13b_2
        cost_coef["rs6"] <- input$cost_dial0_2 - input$cost_ckd13b_2
        cost_coef["rs7"] <- input$cost_dial1_2 - input$cost_ckd13b_2
        cost_coef["cvd2"] <- input$cost_vd_2
        cost_coef["cvd3"] <- input$cost_nvd_2
        cost_coef["cvd5b"] <- input$cost_cvd5_nd_2
        cost_coef["cvd6b"] <- input$cost_cvd6_2
        cost_coef["cvd478b"] <- input$cost_cvd478_2
        cost_coef["dial.cvd5bTRUE"] <- input$cost_cvd5_d_2 - input$cost_cvd5_nd_2
        cost_coef["death_dialTRUE"] <- -(cost_coef["rs6"] + cost_coef["rs7"]) * 0.42 / 2
        cost_coef["DM1"] <- input$cost_diab_2
        coeffs_default$cost2 <- cost_coef
        
        # PSA costs
        
        mu <- input$cost_ckd13b_2
        sigma <- input$cost_ckd13b_sd
        t <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu)
        coeffs_PSA_cost[, "(Intercept)"] <- t
        
        mu <- input$cost_ckd4_2
        sigma <- input$cost_ckd4_sd
        coeffs_PSA_cost[, "rs2"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu) - t
        
        mu <- input$cost_ckd5_2
        sigma <- input$cost_ckd5_sd
        coeffs_PSA_cost[, "rs3"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu) - t
        
        mu <- input$cost_tx0_2
        sigma <- input$cost_tx0_sd
        coeffs_PSA_cost[, "rs4"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu) - t
        
        mu <- input$cost_tx1_2
        sigma <- input$cost_tx1_sd
        coeffs_PSA_cost[, "rs5"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu) - t
        
        mu <- input$cost_dial0_2
        sigma <- input$cost_dial0_sd
        coeffs_PSA_cost[, "rs6"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu) - t
        
        mu <- input$cost_dial1_2
        sigma <- input$cost_dial1_sd
        coeffs_PSA_cost[, "rs7"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu) - t
        
        mu <- input$cost_vd_2
        sigma <- input$cost_vd_sd
        coeffs_PSA_cost[, "cvd2"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu)
        
        mu <- input$cost_nvd_2
        sigma <- input$cost_nvd_sd
        coeffs_PSA_cost[, "cvd3"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu)
        
        mu <- input$cost_cvd5_nd_2
        sigma <- input$cost_cvd5_nd_sd
        coeffs_PSA_cost[, "cvd5b"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu)    
        
        mu <- input$cost_cvd6_2
        sigma <- input$cost_cvd6_sd
        coeffs_PSA_cost[, "cvd6b"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu)
        
        mu <- input$cost_cvd478_2
        sigma <- input$cost_cvd478_sd
        coeffs_PSA_cost[, "cvd478b"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu)
        
        mu <- input$cost_cvd5_d_2
        sigma <- input$cost_cvd5_d_sd
        coeffs_PSA_cost[, "dial.cvd5bTRUE"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu) - 
          coeffs_PSA_cost[, "cvd5b"]
        
        coeffs_PSA_cost[, "death_dialTRUE"] <- -(coeffs_PSA_cost[, "rs6"] + coeffs_PSA_cost[, "rs7"]) * 0.42 / 2
        
        mu <- input$cost_diab_2
        sigma <- input$cost_diab_sd
        coeffs_PSA_cost[, "DM1"] <- rgamma(1000, shape = (mu / sigma)^2, scale = (sigma)^2 / mu)    
        
      }
    
    ################################################################################
    ## Quality of life
    ################################################################################
    
    if (!is.null(input$qol_baseline)) {   
      # screen clicked, overwrite default values
      qol_coef <- coeffs_default$qol2 # extract QoL coeffs
      qol_coef["(Intercept)"] <- input$qol_baseline - 6 * isolate(input$qol_age)
      qol_coef["SEXM"] <- input$qol_sexm
      qol_coef["smokerbefore"] <- input$qol_smokerbefore
      qol_coef["smokercurrently"] <- input$qol_smokercurrently
      qol_coef["educationGCSE/vocational"] <- input$qol_eduGCSE
      qol_coef["educationbelow secondary"] <- input$qol_eduBelowSec
      qol_coef["BMI_quantT1"] <- input$qol_BMI1
      qol_coef["BMI_quantT3"] <- input$qol_BMI3
      qol_coef["TX1"] <- input$qol_transplant
      qol_coef["renalDiagnosisDiabetic nephropathy"] <- input$qol_diabeticNep
      qol_coef["dialysis1"] <- input$qol_dialysis
      qol_coef["DH_CV_5no NFMVE, baseline VD"] <- input$qol_mve2
      qol_coef["DH_CV_5NFMVE <1 year ago"] <- input$qol_mve0
      qol_coef["DH_CV_5NFMVE >1 year ago"] <- input$qol_mve1
      qol_coef["age_T2"] <- input$qol_age/10
      # replace Qol coeffs
      coeffs_default$qol2 <- qol_coef   
    }
    
    ################################################################################
    ### Non-vascular death probabilities
    ################################################################################
    
    if (input$cb_nvd == TRUE) {
      
      inFile2 <- input$file2
      ratesNVD <- read.csv(inFile2$datapath, check.names = F)
      
      # create T0 & T1 columns
      T <- strsplit(as.character(ratesNVD$ageBand), "-")
      ratesNVD$T0 <- as.numeric(lapply(T, "[[", 1))
      ratesNVD$T1 <- as.numeric(lapply(T, "[[", 2))
      
      # make CKD stages numeric to subsequently save as a matrix
      ratesNVD$N_CKD <- NA
      ratesNVD$N_CKD[ratesNVD$CKDStage == 0] <- 1
      ratesNVD$N_CKD[ratesNVD$CKDStage == 1] <- 2
      ratesNVD$N_CKD[ratesNVD$CKDStage == 2] <- 3
      ratesNVD$N_CKD[ratesNVD$CKDStage == 3] <- 8
      ratesNVD$N_CKD[ratesNVD$CKDStage == 4] <- 4
      
      colnames(ratesNVD)[which(colnames(ratesNVD) == "sex")] <- "SEX"
      ratesNVD <- subset(ratesNVD, select = c("N_CKD", "T0", "T1", "SEX", "p_NVD")) 
      
      # split by sex & ESRD status
      ratesNVD_M <- as.matrix(subset(ratesNVD, SEX == 1, select = c("N_CKD", "T0", "T1", "p_NVD")))
      ratesNVD_M_ESRD <- as.matrix(subset(ratesNVD, SEX == 1 & N_CKD %in% c(4, 8), select = c("N_CKD", "T0", "T1", "p_NVD")))
      ratesNVD_F <- as.matrix(subset(ratesNVD, SEX == 0, select = c("N_CKD", "T0", "T1", "p_NVD")))
      ratesNVD_F_ESRD <- as.matrix(subset(ratesNVD, SEX == 0 & N_CKD %in% c(4, 8), select = c("N_CKD", "T0", "T1", "p_NVD")))
      
      # save
      ratesNVD <- list(ratesNVD = ratesNVD, ratesNVD_M = ratesNVD_M, ratesNVD_M_ESRD = ratesNVD_M_ESRD,
                       ratesNVD_F = ratesNVD_F, ratesNVD_F_ESRD = ratesNVD_F_ESRD)
      
    }
    
    
    ################################################################################
    ### Decision paramenters 
    ################################################################################
    
    if (is.null(input$dr_ev)) {
      # screen was not clicked so load default values
      disc_events <- 0.035
      disc_costs <- 0.035
      years <- 30
    } else{
      # screen clicked and can read off the values
      disc_events <- input$dr_ev/100
      disc_costs <- input$dr_cost/100
      if (input$dp_lifetime == "Simulation for a fixed number of years") 
        years <- floor(input$dp_cycles) else {
          lifetime <- TRUE
          years <- input$dp_stopping_age
        }
    }
    
    ################################################################################
    ### Matrix with initial probabilities
    ################################################################################
    
    P_0 <- data.frame(id = pat0[, "id"])
    P_0$ESRDDuration_T <- patT[, "ESRDDuration_T"]
    P_0$cycle <- 0
    P_0$stageRand_2 <- stageRand_2
    
    # add probability columns
    P_0[, states_and_endpoints$endpts_CV] <- 0
    P_0[, states_and_endpoints$endpts_CKD] <- 0
    
    # initialise states
    labs <- unique(sapply(states_and_endpoints$states_info, "[[", 1))
    P_0[, labs] <- 0
    for (i in 1 : nrow(P_0))
      P_0[i, paste(CVRand_num[i], stageRand_2_num[i], sep = "_")] <- 1
    
    # remove redundant columns
    P_0 <- P_0[, !(colnames(P_0) %in% c("stageRand_2", "ESRDDuration_T"))]
    
    # add the QALY columns
    P_0$QALY_NF <- 0
    P_0$QALY_F <- 0
    
    # add the cost column
    P_0$cost_hosp <- 0
    
    ### save as matrix
    P_0 <- as.matrix(P_0)
    
    # default values
    
    output_LE <- NULL
    output_LE_PSA <- NULL
    output_CE <- NULL
    output_CE_PSA <- NULL
    
    single_pat <- nrow(data) == 1
    
    ### TODO: align column names so that don't need to change
    # eg ESRD -> RRT; NFMVE -> MVE
    
    withProgress(message = 'Please wait...', value = 0, {
      
      if (LE) {
        ################################################################################
        ### long-term projections
        ################################################################################
        
        ### run deterministic analysis first
        
        alpha <- wrapper_LE(N_cores = N_cores, 
                            dfBaselineMM_0 = pat0, dfBaselineMM_T = patT, 
                            states_and_endpoints = states_and_endpoints, 
                            stageRand_2 = stageRand_2, stageRand_2_num = stageRand_2_num, 
                            BVD = BVD, CVRand_num = CVRand_num,
                            ratesNVD = ratesNVD, P_0 = P_0, 
                            coeffs = coeffs_default,
                            compl = compl_T, Tx_price = Tx_price_T,
                            disc_events = disc_events, disc_costs = disc_costs,
                            lifetime = lifetime, years = years)
        
        ## patient-level summary - available as download only
        results_ind <- alpha$output_ind
        results_download_ind <- merge(data, results_ind, sort = FALSE)
        
        # summary at group level - available to view and as download
        if (single_pat) {
          results_group <- results_ind
          results_download_group <- results_download_ind
        } else {
          results_group <- alpha$output_group
          results_download_group <- results_group
        }
        
        results_view <- .output_LE(df = results_download_group, lifetime = lifetime)
        
        if (PSA == 0 | is.null(input$PSA)) {
          
          ### change column names
          
          colnames(results_download_ind) <- gsub("ESRD", "RRT", colnames(results_download_ind), fixed = TRUE)
          colnames(results_download_group) <- gsub("ESRD", "RRT", colnames(results_download_group), fixed = TRUE)
          
          colnames(results_download_ind) <- gsub("NFMVE", "MVE", colnames(results_download_ind), fixed = TRUE)
          colnames(results_download_group) <- gsub("NFMVE", "MVE", colnames(results_download_group), fixed = TRUE)
          
          ### output for the deterministic analysis ###
          
          output_LE <- list(results_download_ind = results_download_ind,
                            results_download_group = results_download_group, 
                            results_view = results_view)
          
        } else {
          
          ### analyses and output for the PSA ###
          
          Nsim <- input$Nsamp
          sims <- sort(sample(1 : 1000, Nsim, replace = F))
          
          alpha_PSA <- wrapper_LE_PSA(N_cores = N_cores, 
                                      dfBaselineMM_0 = pat0, dfBaselineMM_T = patT, 
                                      states_and_endpoints = states_and_endpoints, 
                                      stageRand_2 = stageRand_2, stageRand_2_num = stageRand_2_num, 
                                      BVD = BVD, CVRand_num = CVRand_num,
                                      ratesNVD = ratesNVD, P_0 = P_0, 
                                      sims = sims, 
                                      coeffs_PSA_VD = coeffs_PSA_VD,
                                      coeffs_PSA_NFMAEorVD = coeffs_PSA_NFMAEorVD,
                                      coeffs_PSA_NFMAEorVD_gamma = coeffs_PSA_NFMAEorVD_gamma,
                                      coeffs_PSA_NFMVEorVD = coeffs_PSA_NFMVEorVD,
                                      coeffs_PSA_NFMVEorVD_gamma = coeffs_PSA_NFMVEorVD_gamma,
                                      coeffs_PSA_nonesrd = coeffs_PSA_nonesrd,
                                      coeffs_PSA_dial2trans = coeffs_PSA_dial2trans,
                                      coeffs_PSA_trans2dial = coeffs_PSA_trans2dial,
                                      coeffs_PSA_cost = coeffs_PSA_cost,
                                      coeffs_PSA_qol = coeffs_PSA_qol,
                                      compl = compl_T, Tx_price = Tx_price_T,
                                      disc_events = disc_events, disc_costs = disc_costs,
                                      lifetime = lifetime, years = years)
          
          time_suffixes <- if (lifetime) c("_5", "_10") else c("_5", "_10", "_all")
          event_suffixes <- as.vector(t(outer(time_suffixes, c("", "_l", "_u"), paste, sep = "")))
          events <- as.vector(t(outer(c("NFMVEorVD_first", "ESRD_first", "VD", "D"), event_suffixes, paste, sep = "")))
          LYs <- as.vector(t(outer(c("LY", "QALY", "cost_hosp", "cost_tx"), c("", "_l", "_u"), paste, sep = "")))
          cols <- c("id", events, LYs)
          
          # summary at individual level
          results_PSA_ind <- subset(alpha_PSA$dfCI, id != "all")
          colnames(results_PSA_ind) <- gsub(".", "_", colnames(results_PSA_ind), fixed = TRUE)
          results_PSA_ind <- merge(results_ind, results_PSA_ind, sort = FALSE)
          results_PSA_ind <- results_PSA_ind[, cols]
          results_PSA_download_ind <- merge(data, results_PSA_ind, sort = FALSE)
          
          # summary at group level
          if (single_pat)
            results_PSA_download_group <- results_PSA_download_ind else {
              results_PSA_group <- subset(alpha_PSA$dfCI, id == "all")
              colnames(results_PSA_group) <- gsub(".", "_", colnames(results_PSA_group), fixed = TRUE)
              results_PSA_group <- merge(results_group, results_PSA_group, sort = FALSE)
              results_PSA_download_group <- results_PSA_group[, cols]
            }
          
          results_PSA_view <- .output_LE_PSA(df = results_PSA_download_group, lifetime = lifetime)
          
          ### change column names
          
          colnames(results_PSA_download_ind) <- gsub("ESRD", "RRT", colnames(results_PSA_download_ind), fixed = TRUE)
          colnames(results_PSA_download_group) <- gsub("ESRD", "RRT", colnames(results_PSA_download_group), fixed = TRUE)
          
          colnames(results_PSA_download_ind) <- gsub("NFMVE", "MVE", colnames(results_PSA_download_ind), fixed = TRUE)
          colnames(results_PSA_download_group) <- gsub("NFMVE", "MVE", colnames(results_PSA_download_group), fixed = TRUE)
          
          ### return
          
          output_LE_PSA <- list(Nsim = input$Nsamp,
                                results_download_ind = results_PSA_download_ind,
                                results_download_group = results_PSA_download_group,
                                results_view = results_PSA_view)
        }
      } else {
        
        ################################################################################
        ### cost-effectiveness analysis
        ################################################################################
        
        ### run deterministic analysis first
        
        alpha <- wrapper(N_cores = N_cores, 
                         dfBaselineMM_0 = pat0, dfBaselineMM_T = patT, 
                         states_and_endpoints = states_and_endpoints, 
                         stageRand_2 = stageRand_2, stageRand_2_num = stageRand_2_num, 
                         BVD = BVD, CVRand_num = CVRand_num,
                         ratesNVD = ratesNVD, P_0 = P_0, 
                         coeffs = coeffs_default,
                         compl_C = compl_C, compl_T = compl_T,
                         Tx_price_C = Tx_price_C, Tx_price_T = Tx_price_T,
                         disc_events = disc_events, disc_costs = disc_costs,
                         lifetime = lifetime, years = years)
        
        ## patient-level summary 
        results_ind <- alpha$output_ind
        results_download_ind <- merge(data, results_ind, sort = FALSE)
        
        # summary at group level 
        if (single_pat) {
          results_group <- results_ind
          results_download_group <- results_download_ind
        } else {
          results_group <- alpha$output_group
          results_download_group <- results_group
        }
        
        results_view <- .output(df = results_download_group, lifetime = lifetime)
        
        if (PSA == 0 | is.null(input$PSA)) {
          
          ### change column names
          
          colnames(results_download_ind) <- gsub("ESRD", "RRT", colnames(results_download_ind), fixed = TRUE)
          colnames(results_download_group) <- gsub("ESRD", "RRT", colnames(results_download_group), fixed = TRUE)
          
          colnames(results_download_ind) <- gsub("NFMVE", "MVE", colnames(results_download_ind), fixed = TRUE)
          colnames(results_download_group) <- gsub("NFMVE", "MVE", colnames(results_download_group), fixed = TRUE)
          
          ### output for the deterministic analysis ###
          
          output_CE <- list(results_download_ind = results_download_ind,
                            results_download_group = results_download_group, 
                            results_view = results_view)
        } else {
          
          ### analyses and output for the PSA ###
          
          Nsim <- input$Nsamp
          sims <- sort(sample(1 : 1000, Nsim, replace = F))
          
          alpha_PSA <- wrapper_PSA(N_cores = N_cores, 
                                   dfBaselineMM_0 = pat0, dfBaselineMM_T = patT, 
                                   states_and_endpoints = states_and_endpoints, 
                                   stageRand_2 = stageRand_2, stageRand_2_num = stageRand_2_num, 
                                   BVD = BVD, CVRand_num = CVRand_num,
                                   ratesNVD = ratesNVD, P_0 = P_0, 
                                   sims = sims, 
                                   coeffs_PSA_VD = coeffs_PSA_VD,
                                   coeffs_PSA_NFMAEorVD = coeffs_PSA_NFMAEorVD,
                                   coeffs_PSA_NFMAEorVD_gamma = coeffs_PSA_NFMAEorVD_gamma,
                                   coeffs_PSA_NFMVEorVD = coeffs_PSA_NFMVEorVD,
                                   coeffs_PSA_NFMVEorVD_gamma = coeffs_PSA_NFMVEorVD_gamma,
                                   coeffs_PSA_nonesrd = coeffs_PSA_nonesrd,
                                   coeffs_PSA_dial2trans = coeffs_PSA_dial2trans,
                                   coeffs_PSA_trans2dial = coeffs_PSA_trans2dial,
                                   coeffs_PSA_cost = coeffs_PSA_cost,
                                   coeffs_PSA_qol = coeffs_PSA_qol,
                                   compl_C = compl_C, compl_T = compl_T,
                                   Tx_price_C = Tx_price_C, Tx_price_T = Tx_price_T,
                                   disc_events = disc_events, disc_costs = disc_costs,
                                   lifetime = lifetime, years = years)
          time_suffixes <- if (lifetime) c("_5", "_10") else c("_5", "_10", "_all")
          events <- as.vector(t(outer(c("NFMVEorVD_first", "ESRD_first", "VD", "D"), time_suffixes, paste, sep = "")))
          LYs <- as.vector(t(outer(c("LY", "QALY", "cost_hosp", "cost_tx"), 
                                   c("", "_disc"), paste, sep = "")))
          
          # control group
          events_C <- paste(events, "C", sep = "_")
          LY_C <- paste(LYs, "C", sep = "_")
          cols_C <- as.vector(t(outer(c(events_C, LY_C), c("", "_l", "_u"), paste, sep = "")))
          
          # treatment group
          events_T <- paste(events, "T", sep = "_")
          LY_T <- paste(LYs, "T", sep = "_")
          cols_T <- as.vector(t(outer(c(events_T, LY_T), c("", "_l", "_u"), paste, sep = "")))
          
          # incremental columns
          cols_inc <- as.vector(t(outer(paste(c("LY", "QALY", "cost_hosp", "cost_tx", "cost_total"), "inc", sep = "_"), 
                                        c("", "_disc"), paste, sep = "")))
          cols_inc <- as.vector(t(outer(cols_inc, c("", "_l", "_u"), paste, sep = "")))
          
          # icers
          cols_ICER <- c(c("cost_LY", "cost_QALY"),
                         paste(c("cost_LY", "cost_QALY"), "disc", sep = "_"))
          cols_ICER <- as.vector(t(outer(cols_ICER, c("", "_l", "_u"), paste, sep = "")))
          
          # combine
          cols <- c("id", cols_C, cols_T, cols_inc, cols_ICER)
          
          # summary at individual level
          results_PSA_ind <- subset(alpha_PSA$dfCI, id != "all")
          colnames(results_PSA_ind) <- gsub(".", "_", colnames(results_PSA_ind), fixed = TRUE)
          results_PSA_ind <- merge(results_ind, results_PSA_ind, sort = FALSE)
          results_PSA_ind <- results_PSA_ind[, cols]
          results_PSA_download_ind <- merge(data, results_PSA_ind, sort = FALSE)
          
          if (single_pat)
            results_PSA_download_group <- results_PSA_download_ind else {
              results_PSA_group <- subset(alpha_PSA$dfCI, id == "all")
              colnames(results_PSA_group) <- gsub(".", "_", colnames(results_PSA_group), fixed = TRUE)
              results_PSA_group <- merge(results_group, results_PSA_group, sort = FALSE)
              results_PSA_download_group <- results_PSA_group[, cols]
            }
          
          results_PSA_view <- .output_PSA(df = results_PSA_download_group, lifetime = lifetime)
          
          # CEAC
          dfCEAC <- alpha_PSA$dfCEAC
          p_CEAC_undisc <- .output_p_CEAC(df = dfCEAC, CE_lab = "CE_undisc")
          p_CEAC_disc <- .output_p_CEAC(df = dfCEAC, CE_lab = "CE_disc")
          
          ### change column names
          
          colnames(results_PSA_download_ind) <- gsub("ESRD", "RRT", colnames(results_PSA_download_ind), fixed = TRUE)
          colnames(results_PSA_download_group) <- gsub("ESRD", "RRT", colnames(results_PSA_download_group), fixed = TRUE)
          
          colnames(results_PSA_download_ind) <- gsub("NFMVE", "MVE", colnames(results_PSA_download_ind), fixed = TRUE)
          colnames(results_PSA_download_group) <- gsub("NFMVE", "MVE", colnames(results_PSA_download_group), fixed = TRUE)
          
          # return output
          output_CE_PSA <- list(Nsim = Nsim,
                                results_download_ind = results_PSA_download_ind,
                                results_download_group = results_PSA_download_group,
                                results_view = results_PSA_view,
                                p_CEAC_undisc = p_CEAC_undisc,
                                p_CEAC_disc = p_CEAC_disc)
        }
      }
      
    })
    
    return(list(single_pat = single_pat,
                output_LE = output_LE, output_LE_PSA = output_LE_PSA,
                output_CE = output_CE, output_CE_PSA = output_CE_PSA))
  })
  
  ################################################################################
  ### Outputs: common
  ################################################################################
  
  # TODO: one command for all downloads
  # possibly one command for long-term projections (ie PSA and not PSA in one file) & one for CE
  
  output$text_output_detail <- renderUI({
    if (is.null(newdata()))
      return(NULL)
    if (input$anal_type == "Long-term projections" & input$PSA == "No (deterministic analysis)")
      if (is.null(newdata()$output_LE))
        return(NULL)
    if (input$anal_type == "Long-term projections" & input$PSA == "Yes (probabilistic analysis)")
      if (is.null(newdata()$output_LE_PSA))
        return(NULL)
    if (input$anal_type == "Cost-effectiveness analysis" & input$PSA == "No (deterministic analysis)")
      if (is.null(newdata()$output_CE))
        return(NULL)
    if (input$anal_type == "Cost-effectiveness analysis" & input$PSA == "Yes (probabilistic analysis)")
      if (is.null(newdata()$output_CE_PSA))
        return(NULL)
    if (newdata()$single_pat)
      return(list(
        br(),
        p(em("Detailed results are available in the downloadable summary file.")),
        br(),
        p(em("Probability of a major vascular event or vascular death is only available for participants entering model
        without a history of a major atherosclerotic event or haemorrhagic stroke."))))
    return(list(
      br(),
      p(em("The results are presented at the group level. 
        Detailed and patient-level results are available in the downloadable summary files.")),
      br(),
      p(em("Probability of a major vascular event or vascular death is only available for participants entering model
        without a history of a major atherosclerotic event or haemorrhagic stroke."))))
  })
  
  ################################################################################
  ### Outputs: long-term projections
  ################################################################################
  
  ### deterministic analysis ###
  
  output$table_LE <- renderTable(
    {newdata()$output_LE$results_view[["Tx"]]}, 
    caption = 'Long-term projections (cumulative probabilities per 1,000 participants)',
    caption.placement = 'top',
    digits = 0,
    align = 'r',
    include.rownames = FALSE
  )
  
  output$link_to_download_LE <- renderUI({
    times <- input$action1
    if (is.null(newdata()$output_LE))
      return(NULL) 
    if (newdata()$single_pat)
      return(list(downloadLink("download_LE", "Download summary")))
    return(list(downloadLink("download_LE", "Download summary"),
                br(),
                downloadLink("download_ind_LE", "Download patient-level summary")))
  })
  
  #download group summary (NB only works in an external browser)
  output$download_LE <- downloadHandler(
    filename = paste0("download_", Sys.Date(),".csv"),
    content = function(file) {
      write.csv(newdata()$output_LE$results_download_group, file, row.names = F)}
  )
  
  #download patient-level summary (NB only works in an external browser)
  output$download_ind_LE <- downloadHandler(
    filename = paste0("download_", Sys.Date(),".csv"),
    content = function(file) {
      write.csv(newdata()$output_LE$results_download_ind, file, row.names = F)}
  )
  
  ### probabilistic analysis ###
  
  output$table_LE_PSA <- renderTable(
    {newdata()$output_LE_PSA$results_view[["Tx"]]}, 
    caption = 'Long-term projections (cumulative probabilities per 1,000 participants)',
    caption.placement = 'top',
    digits = 0,
    align = 'r',
    include.rownames = FALSE)
  
  output$link_to_download_LE_PSA <- renderUI({
    times <- input$action1
    if (is.null(newdata()$output_LE_PSA))
      return(NULL) 
    if (nrow(newdata()$output_LE_PSA$results_download_ind) == 1)
      return(list(textOutput("text_PSA_Nsim"),
                  br(),
                  downloadLink("download_LE_PSA", "Download summary")))
    if (nrow(newdata()$output_LE_PSA$results_download_ind) > 1)
      return(list(textOutput("text_PSA_Nsim"),
                  downloadLink("download_LE_PSA", "Download summary"),
                  br(),
                  downloadLink("download_ind_LE_PSA", "Download patient-level summary")))
  })
  
  #download group summary (NB only works in an external browser)
  output$download_LE_PSA <- downloadHandler(
    filename = paste0("download_", Sys.Date(),".csv"),
    content = function(file) {
      write.csv(newdata()$output_LE_PSA$results_download_group, file, row.names = F)}
  )
  
  #download patient-level summary (NB only works in an external browser)
  output$download_ind_LE_PSA <- downloadHandler(
    filename = paste0("download_", Sys.Date(),".csv"),
    content = function(file) {
      write.csv(newdata()$output_LE_PSA$results_download_ind, file, row.names = F)}
  )
  
  ################################################################################
  ### Outputs: Cost-effectiveness analysis
  ################################################################################
  
  output$display_disc <- renderUI({
    times <- input$action1
    if (input$PSA == "No (deterministic analysis)" & is.null(newdata()$output_CE) |
        input$PSA == "Yes (probabilistic analysis)" & is.null(newdata()$output_CE_PSA))
      return(NULL) 
    checkboxInput("disc", 
                  "Discount cost-effectiveness results")
  })
  
  ### deterministic analysis ###
  
  output$table_Tx_C <- renderTable(
    {newdata()$output_CE$results_view[["Tx_C"]]}, 
    caption = 'Long-term projections in the control group (cumulative probabilities per 1,000 participants)',
    caption.placement = 'top',
    digits = 0,
    align = 'r',
    include.rownames = FALSE)
  
  output$table_Tx_T <- renderTable(
    {newdata()$output_CE$results_view[["Tx_T"]]}, 
    caption = 'Long-term projections in the treatment group (cumulative probabilities per 1,000 participants)',
    caption.placement = 'top',
    digits = 0,
    align = 'r',
    include.rownames = FALSE)
  
  output$table_CE_undisc <- renderTable(
    {newdata()$output_CE$results_view[["CE_undisc"]]}, 
    caption = 'Incremental cost-effectiveness over the simulation duration (results per 1,000 participants)',
    caption.placement = 'top',
    digits = 0,
    align = 'r',
    include.rownames = FALSE)
  
  output$table_CE_disc <- renderTable(
    {newdata()$output_CE$results_view[["CE_disc"]]}, 
    caption = 'Incremental cost-effectiveness over the simulation duration (results per 1,000 participants)',
    caption.placement = 'top',
    digits = 0,
    align = 'r',
    include.rownames = FALSE)
  
  output$link_to_download <- renderUI({
    times <- input$action1
    if (is.null(newdata()$output_CE))
      return(NULL) 
    if (newdata()$single_pat) 
      return(list(downloadLink("download", "Download summary")))
    return(list(downloadLink("download", "Download summary"),
                br(),
                downloadLink("download_ind", "Download patient-level summary")))
  })
  
  #download group summary (NB only works in an external browser)
  output$download <- downloadHandler(
    filename = paste0("download_", Sys.Date(),".csv"),
    content = function(file) {
      write.csv(newdata()$output_CE$results_download_group, file, row.names = F)}
  )
  
  #download patient-level summary (NB only works in an external browser)
  output$download_ind <- downloadHandler(
    filename = paste0("download_", Sys.Date(),".csv"),
    content = function(file) {
      write.csv(newdata()$output_CE$results_download_ind, file, row.names = F)}
  )
  
  ### probabilistic analysis ###
  
  output$table_Tx_C_PSA <- renderTable(
    {newdata()$output_CE_PSA$results_view[["Tx_C"]]}, 
    caption = 'Long-term projections in the control group (cumulative probabilities per 1,000 participants)',
    caption.placement = 'top',
    digits = 0,
    align = 'r',
    include.rownames = FALSE)
  
  output$table_Tx_T_PSA <- renderTable(
    {newdata()$output_CE_PSA$results_view[["Tx_T"]]}, 
    caption = 'Long-term projections in the treatment group (cumulative probabilities per 1,000 participants)',
    caption.placement = 'top',
    digits= 0,
    align = 'r',
    include.rownames = FALSE)
  
  output$table_CE_undisc_PSA <- renderTable(
    {newdata()$output_CE_PSA$results_view[["CE_undisc"]]}, 
    caption = 'Incremental cost-effectiveness over the simulation duration (results per 1,000 participants)',
    caption.placement = 'top',
    digits = 0,
    align = 'r',
    include.rownames = FALSE)
  
  output$table_CE_disc_PSA <- renderTable(
    {newdata()$output_CE_PSA$results_view[["CE_disc"]]}, 
    caption = 'Incremental cost-effectiveness over the simulation duration (results per 1,000 participants)',
    caption.placement = 'top',
    digits= 0,
    align = 'r',
    include.rownames = FALSE)
  
  output$p_CEAC_undisc <- renderPlot(
    print(newdata()$output_CE_PSA$p_CEAC_undisc))
  
  output$p_CEAC_disc <- renderPlot(
    print(newdata()$output_CE_PSA$p_CEAC_disc))
  
  output$link_to_download_PSA <- renderUI({
    times <- input$action1
    if (is.null(newdata()$output_CE_PSA))
      return(NULL) 
    if (newdata()$single_pat)
      return(list(textOutput("text_PSA_Nsim"),
                  br(),
                  downloadLink("download_PSA", "Download summary")))
    return(list(textOutput("text_PSA_Nsim"),
                downloadLink("download_PSA", "Download summary"),
                br(),
                downloadLink("download_PSA_ind", "Download patient-level summary")))
  })
  
  #download group summary (NB only works in an external browser)
  output$download_PSA <- downloadHandler(
    filename = paste0("download_", Sys.Date(),".csv"),
    content = function(file) {
      write.csv(newdata()$output_CE_PSA$results_download_group, file, row.names = F)}
  )
  
  #download patient-level summary (NB only works in an external browser)
  output$download_PSA_ind <- downloadHandler(
    filename = paste0("download_", Sys.Date(),".csv"),
    content = function(file) {
      write.csv(newdata()$output_CE_PSA$results_download_ind, file, row.names = F)}
  )
  
  })