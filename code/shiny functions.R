###################################################################
### internal functions for cleaning the input files                                     		      
###################################################################

# validating an input while, which should be non-empty for convenience
.validate_file_nonempty <- function(inFile, 
                                    req_colnames, req_format, req_values) {
  
  df <- read.csv(inFile$datapath, stringsAsFactors=FALSE)
  colnames_file <- colnames(df)
  
  # missing variables
  vars_missing <- c()
  for (name in req_colnames)
    if (!(name %in% colnames_file))
      vars_missing <- c(vars_missing, name)
  
  # variables in the wrong format
  vars_format <- list()
  vars_format_names <- c()
  for (name in req_format) {
    var_temp <- name$var
    if (!(var_temp %in% vars_missing)) {
      if (name$numeric == TRUE) {
        if (!is.numeric(df[, var_temp])) {
          vars_format <- c(vars_format, paste(var_temp, " (needs to be numeric)")) 
          vars_format_names <- c(vars_format_names, var_temp)
        }
      } else
        if (!is.character(df[, var_temp])) {
          vars_format <- c(vars_format, paste(var_temp, " (needs to be character)")) 
          vars_format_names <- c(vars_format_names, var_temp)
        }
    }
  }
  
  # variables with disallowed values
  vars_values <- list() 
  vars_values_names <- list()
  for (name in req_values) {
    var_temp <- name$var
    if (!(var_temp %in% vars_missing))
      if (!(var_temp %in% vars_format_names))
        if (length(intersect(name$othervars, vars_missing)) == 0
            & length(intersect(name$othervars, vars_format_names)) == 0)
          if (!eval(parse(text = name$condition))) {
            vars_values <- c(vars_values, paste(var_temp, " (", name$error_text, ")", sep = ""))
            vars_values_names <- c(vars_values_names, var_temp)
          } 
  }
  
  return(list(vars_missing = vars_missing, 
              vars_format_names = vars_format_names, vars_format = vars_format, 
              vars_values = vars_values, vars_values_names = vars_values_names))
}

###################################################################
### internal functions for producing suitable output: long-term projections
###################################################################

# format tables for the printed output (deterministic analysis)
.output_LE <- function(df, lifetime) {
  
  ### outcomes
  
  colnames_Tx <- c("", "MVE or VD", "RRT", "Vascular deaths", "All deaths")
  
  Tx <- data.frame(
    c("At 5 years", "At 10 years"),
    c(1000 * df$NFMVEorVD_first_5, 1000 * df$NFMVEorVD_first_10),
    c(1000 * df$ESRD_first_5, 1000 * df$ESRD_first_10),
    c(1000 * df$VD_5, 1000 * df$VD_10),
    c(1000 * df$D_5, 1000 * df$D_10))
  colnames(Tx) <- colnames_Tx    
  
  if (!(lifetime)) {
    Tx_2 <- data.frame(
      c("Over simulation duration"),
      c(1000 * df$NFMVEorVD_first_all),
      c(1000 * df$ESRD_first_all),
      c(1000 * df$VD_all),
      c(1000 * df$D_all))
    colnames(Tx_2) <- colnames_Tx
    Tx <- rbind(Tx, Tx_2)
  }
  
  return(list(Tx = Tx))
  
}

# format tables for the printed output (probabilistic analysis)
.output_LE_PSA <- function(df, lifetime) {
  
  vars <- list(
    list(var = "NFMVEorVD_first_5", digits = 0, negate = FALSE),
    list(var = "ESRD_first_5", digits = 0, negate = FALSE),
    list(var = "VD_5", digits = 0, negate = FALSE),
    list(var = "D_5", digits = 0, negate = FALSE),
    list(var = "NFMVEorVD_first_10", digits = 0, negate = FALSE),
    list(var = "ESRD_first_10", digits = 0, negate = FALSE),
    list(var = "VD_10", digits = 0, negate = FALSE),
    list(var = "D_10", digits = 0, negate = FALSE))
  
  if (!(lifetime))
    vars <- append(vars, list(
      list(var = "NFMVEorVD_first_all", digits = 0, negate = FALSE),
      list(var = "ESRD_first_all", digits = 0, negate = FALSE),
      list(var = "VD_all", digits = 0, negate = FALSE),
      list(var = "D_all", digits = 0, negate = FALSE)
    ))
  
  for (var_n in vars) {
    
    var <- var_n$var
    dg <- var_n$digits
    var_l <- paste(var, "l", sep = "_")
    var_u <- paste(var, "u", sep = "_")
    temp <- df[, c(var, var_l, var_u)]
    varOld <- paste(var, "Est", sep = "_")
    df[, varOld] <- df[, var]
    if (var_n$negate) {
      #df[, var] <- comma_format(accuracy = dg, nsmall = 0)(-df[, var] * 1000)
      #df[, var_l] <- comma_format(accuracy = dg, nsmall = 0)(-df[, var_l] * 1000)
      #df[, var_u] <- comma_format(accuracy = dg, nsmall = 0)(-df[, var_u] * 1000)
      df[, var] <- round(-df[, var] * 1000)
      df[, var_l] <- round(-df[, var_l] * 1000)
      df[, var_u] <- round(-df[, var_u] * 1000)
      df[, var] <- paste(df[, var], " (", df[, var_u], ", ", df[, var_l], ")", sep = "")
    } else {
      #df[, var] <- comma_format(accuracy = dg, nsmall = 0)(df[, var] * 1000)
      #df[, var_l] <- comma_format(accuracy = dg, nsmall = 0)(df[, var_l] * 1000)
      #df[, var_u] <- comma_format(accuracy = dg, nsmall = 0)(df[, var_u] * 1000)
      df[, var] <- round(df[, var] * 1000)
      df[, var_l] <- round(df[, var_l] * 1000)
      df[, var_u] <- round(df[, var_u] * 1000)
      df[, var] <- paste(df[, var], " (", df[, var_l], ", ", df[, var_u], ")", sep = "")
    }
  }
  
  ### outcomes
  
  colnames_Tx <- c("", "MVE or VD", "RRT", "Vascular deaths", "All deaths")
  
  Tx <- data.frame(
    c("At 5 years", "At 10 years"),
    c(df$NFMVEorVD_first_5, df$NFMVEorVD_first_10),
    c(df$ESRD_first_5, df$ESRD_first_10),
    c(df$VD_5, df$VD_10),
    c(df$D_5, df$D_10))
  colnames(Tx) <- colnames_Tx    
  
  if (!(lifetime)) {
    Tx_2 <- data.frame(
      c("Over simulation duration"),
      c(df$NFMVEorVD_first_all),
      c(df$ESRD_first_all),
      c(df$VD_all),
      c(df$D_all))
    colnames(Tx_2) <- colnames_Tx
    
    Tx <- rbind(Tx, Tx_2)
  }
  
  return(list(Tx = Tx))
  
}


###################################################################
### internal functions for producing suitable output: CE analysis
###################################################################


# format tables for the printed output (deterministic analysis)
.output <- function(df, lifetime) {
  
  ### outcomes
  
  colnames_Tx <- c("", "MVE or VD", "RRT", "Vascular deaths", "All deaths")
  
  # Placebo group (Treatment A)
  Tx_C <- data.frame(
    c("At 5 years", "At 10 years"),
    c(1000 * df$NFMVEorVD_first_5_C, 1000 * df$NFMVEorVD_first_10_C),
    c(1000 * df$ESRD_first_5_C, 1000 * df$ESRD_first_10_C),
    c(1000 * df$VD_5_C, 1000 * df$VD_10_C),
    c(1000 * df$D_5_C, 1000 * df$D_10_C))
  colnames(Tx_C) <- colnames_Tx
  
  # Treatment group (Treatment B)
  Tx_T <- data.frame(
    c("At 5 years", "At 10 years"),
    c(1000 * df$NFMVEorVD_first_5_T, 1000 * df$NFMVEorVD_first_10_T),
    c(1000 * df$ESRD_first_5_T, 1000 * df$ESRD_first_10_T),
    c(1000 * df$VD_5_T, 1000 * df$VD_10_T),
    c(1000 * df$D_5_T, 1000 * df$D_10_T))
  colnames(Tx_T) <- colnames_Tx
  
  # add "all" column if required
  if (!(lifetime)) {
    
    Tx_C_2 <- data.frame(
      c("Over simulation duration"),
      c(1000 * df$NFMVEorVD_first_all_C),
      c(1000 * df$ESRD_first_all_C),
      c(1000 * df$VD_all_C),
      c(1000 * df$D_all_C))
    colnames(Tx_C_2) <- colnames_Tx
    Tx_C <- rbind(Tx_C, Tx_C_2)
    
    Tx_T_2 <- data.frame(
      c("Over simulation duration"),
      c(1000 * df$NFMVEorVD_first_all_T),
      c(1000 * df$ESRD_first_all_T),
      c(1000 * df$VD_all_T),
      c(1000 * df$D_all_T))
    colnames(Tx_T_2) <- colnames_Tx
    Tx_T <- rbind(Tx_T, Tx_T_2)
  }
  
  ### CE results
  
  colnames_CE <- c("LYs gained", "QALYs gained", "Incremental hospital costs",
                   "Treatment costs",
                   "Cost per LY gained", "Cost per QALY gained")
  
  # undiscounted
  CE_undisc <- data.frame(
    1000 * df$LY_inc, 1000 * df$QALY_inc,
    comma_format(accuracy = 0)(1000 * (df$cost_hosp_T - df$cost_hosp_C)),
    comma_format(accuracy = 0)(1000 * (df$cost_tx_T - df$cost_tx_C)),
    comma_format(accuracy = 0)(df$cost_LY), 
    comma_format(accuracy = 0)(df$cost_QALY))
  colnames(CE_undisc) <- colnames_CE
  
  # discounted
  CE_disc <- data.frame(
    1000 * df$LY_inc_disc, 1000 * df$QALY_inc_disc,
    comma_format(accuracy = 0)(1000 * (df$cost_hosp_disc_T - df$cost_hosp_disc_C)),
    comma_format(accuracy = 0)(1000 * (df$cost_tx_disc_T - df$cost_tx_disc_C)),
    comma_format(accuracy = 0)(df$cost_LY_disc), 
    comma_format(accuracy = 0)(df$cost_QALY_disc))
  colnames(CE_disc) <- colnames_CE
  
  return(list(Tx_C = Tx_C, Tx_T = Tx_T, 
              CE_undisc = CE_undisc, CE_disc = CE_disc))
}

# format tables for the printed output (probabilistic analysis)
.output_PSA <- function(df, lifetime) {
  
  vars <- list(
    list(var = "NFMVEorVD_first_5_C", scale = 1000, digits = 0, negate = FALSE),
    list(var = "ESRD_first_5_C", scale = 1000, digits = 0, negate = FALSE),
    list(var = "VD_5_C", scale = 1000, digits = 0, negate = FALSE),
    list(var = "D_5_C", scale = 1000, digits = 0, negate = FALSE),
    list(var = "NFMVEorVD_first_10_C", scale = 1000, digits = 0, negate = FALSE),
    list(var = "ESRD_first_10_C", scale = 1000, digits = 0, negate = FALSE),
    list(var = "VD_10_C", scale = 1000, digits = 0, negate = FALSE),
    list(var = "D_10_C", scale = 1000, digits = 0, negate = FALSE),
    list(var = "NFMVEorVD_first_5_T", scale = 1000, digits = 0, negate = FALSE),
    list(var = "ESRD_first_5_T", scale = 1000, digits = 0, negate = FALSE),
    list(var = "VD_5_T", scale = 1000, digits = 0, negate = FALSE),
    list(var = "D_5_T", scale = 1000, digits = 0, negate = FALSE),
    list(var = "NFMVEorVD_first_10_T", scale = 1000, digits = 0, negate = FALSE),
    list(var = "ESRD_first_10_T", scale = 1000, digits = 0, negate = FALSE),
    list(var = "VD_10_T", scale = 1000, digits = 0, negate = FALSE),
    list(var = "D_10_T", scale = 1000, digits = 0, negate = FALSE),
    list(var = "LY_inc", scale = 1000, digits = 0, negate = FALSE),
    list(var = "QALY_inc", scale = 1000, digits = 0, negate = FALSE),
    list(var = "LY_inc_disc", scale = 1000, digits = 0, negate = FALSE),
    list(var = "QALY_inc_disc", scale = 1000, digits = 0, negate = FALSE),
    list(var = "cost_hosp_inc", scale = 1000, digits = 0, negate = FALSE),
    list(var = "cost_hosp_inc_disc", scale = 1000, digits = 0, negate = FALSE),
    list(var = "cost_tx_inc", scale = 1000, digits = 0, negate = FALSE),
    list(var = "cost_tx_inc_disc", scale = 1000, digits = 0, negate = FALSE),
    list(var = "cost_LY", scale = 1, digits = 0, negate = FALSE),
    list(var = "cost_QALY", scale = 1, digits = 0, negate = FALSE),
    list(var = "cost_LY_disc", scale = 1, digits = 0, negate = FALSE),
    list(var = "cost_QALY_disc", scale = 1, digits = 0, negate = FALSE))
  
  if (!(lifetime))
    vars <- append(vars, list(
      list(var = "NFMVEorVD_first_all_C", scale = 1000, digits = 0, negate = FALSE),
      list(var = "ESRD_first_all_C", scale = 1000, digits = 0, negate = FALSE),
      list(var = "VD_all_C", scale = 1000, digits = 0, negate = FALSE),
      list(var = "D_all_C", scale = 1000, digits = 0, negate = FALSE),
      list(var = "NFMVEorVD_first_all_T", scale = 1000, digits = 0, negate = FALSE),
      list(var = "ESRD_first_all_T", scale = 1000, digits = 0, negate = FALSE),
      list(var = "VD_all_T", scale = 1000, digits = 0, negate = FALSE),
      list(var = "D_all_T", scale = 1000, digits = 0, negate = FALSE)))
  
  for (var_n in vars) {
    var <- var_n$var
    scale <- var_n$scale
    
    dg <- var_n$digits
    var_l <- paste(var, "l", sep = "_")
    var_u <- paste(var, "u", sep = "_")
    temp <- df[, c(var, var_l, var_u)]
    varOld <- paste(var, "Est", sep = "_")
    df[, varOld] <- df[, var]
    if (var_n$negate) {
      #df[, var] <- comma_format(accuracy = dg, nsmall = 0)(-df[, var] * scale)
      #df[, var_l] <- comma_format(accuracy = dg, nsmall = 0)(-df[, var_l] * scale)
      #df[, var_u] <- comma_format(accuracy = dg, nsmall = 0)(-df[, var_u] * scale)
      df[, var] <- round(-df[, var] * scale)
      df[, var_l] <- round(-df[, var_l] * scale)
      df[, var_u] <- round(-df[, var_u] * scale)
      df[, var] <- paste(df[, var], " (", df[, var_u], ", ", df[, var_l], ")", sep = "")
    } else {
      #df[, var] <- comma_format(accuracy = dg, nsmall = 0)(df[, var] * scale)
      #df[, var_l] <- comma_format(accuracy = dg, nsmall = 0)(df[, var_l] * scale)
      #df[, var_u] <- comma_format(accuracy = dg, nsmall = 0)(df[, var_u] * scale)
      df[, var] <- round(df[, var] * scale)
      df[, var_l] <- round(df[, var_l] * scale)
      df[, var_u] <- round(df[, var_u] * scale)
      df[, var] <- paste(df[, var], " (", df[, var_l], ", ", df[, var_u], ")", sep = "")
    }
  }
  
  ### outcomes
  
  colnames_Tx <- c("", "MVE or VD", "RRT", "Vascular deaths", "All deaths")
  
  # Placebo group (Treatment A)
  Tx_C <- data.frame(
    c("At 5 years", "At 10 years"),
    c(df$NFMVEorVD_first_5_C, df$NFMVEorVD_first_10_C),
    c(df$ESRD_first_5_C, df$ESRD_first_10_C),
    c(df$VD_5_C, df$VD_10_C),
    c(df$D_5_C, df$D_10_C))
  colnames(Tx_C) <- colnames_Tx  
  
  # Treatment group (Treatment B)
  Tx_T <- data.frame(
    c("At 5 years", "At 10 years"),
    c(df$NFMVEorVD_first_5_T, df$NFMVEorVD_first_10_T),
    c(df$ESRD_first_5_T, df$ESRD_first_10_T),
    c(df$VD_5_T, df$VD_10_T),
    c(df$D_5_T, df$D_10_T))
  colnames(Tx_T) <- colnames_Tx    
  
  # add "all" column if required
  if (!(lifetime)) {
    
    Tx_C_2 <- data.frame(
      c("Over simulation duration"),
      c(df$NFMVEorVD_first_all_C),
      c(df$ESRD_first_all_C),
      c(df$VD_all_C),
      c(df$D_all_C))
    colnames(Tx_C_2) <- colnames_Tx
    Tx_C <- rbind(Tx_C, Tx_C_2)
    
    Tx_T_2 <- data.frame(
      c("Over simulation duration"),
      c(df$NFMVEorVD_first_all_T),
      c(df$ESRD_first_all_T),
      c(df$VD_all_T),
      c(df$D_all_T))
    colnames(Tx_T_2) <- colnames_Tx
    Tx_T <- rbind(Tx_T, Tx_T_2)
  }
  
  ### CE results
  
  colnames_CE <- c("LYs gained", "QALYs gained", "Incremental hospital costs",
                   "Treatment costs",
                   "Cost per LY gained", "Cost per QALY gained")
  
  # undiscounted
  CE_undisc <- data.frame(
    df$LY_inc, df$QALY_inc,
    df$cost_hosp_inc, df$cost_tx_inc, 
    df$cost_LY, df$cost_QALY)
  colnames(CE_undisc) <- colnames_CE
  
  # discounted
  CE_disc <- data.frame(
    df$LY_inc_disc, df$QALY_inc_disc,
    df$cost_hosp_inc_disc, df$cost_tx_inc_disc,
    df$cost_LY_disc, df$cost_QALY_disc)
  colnames(CE_disc) <- colnames_CE
  
#print(Tx_C)  
#print(Tx_T)
#print(CE_undisc)
#print(CE_disc)
  
  return(list(Tx_C = Tx_C, Tx_T = Tx_T, 
              CE_undisc = CE_undisc, CE_disc = CE_disc))
  
}

.output_p_CEAC <- function(df, CE_lab){
  
  p <- ggplot(df, aes_string(x = "R", y = CE_lab))
  #p <- p + geom_smooth(se = FALSE)
  p <- p + geom_line()
  p <- p + facet_wrap( ~ id, ncol = 2)
  p <- p + scale_x_continuous(name = "Value of cost-effectiveness threshold", 
                              breaks = seq(from = 10000, to = 100000, by = 10000),
                              labels = comma)
  p <- p + scale_y_continuous(name = "Probability cost-effective", 
                              breaks = seq(from = 0.10, to = 1, by = 0.1), 
                              labels = percent)
  p <- p + coord_cartesian(ylim = c(0, 1.05))
  p <- p + theme_bw() + theme(
    axis.text.x = element_text(size = 13, angle = 30),
    axis.text.y = element_text(size = 13),
    strip.text = element_text(size = 15),
    strip.background = element_blank(),
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.7))
  
  return(p)
  
}