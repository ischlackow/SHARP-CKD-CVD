#This file provides code to form the user interface - it defines the inputs into 
#the model sent to 'server' and calls the resulting output. 

#to run shiny both 'ui.R' and 'server.R' must be open (and must have these precise names)
#there will be a button in the top right called 'Run App'. Alternatively you can 
#press CTRL+SHIFT+ENTER. 

library(shiny)
library(shinydashboard)

shinyUI(fluidPage(
  
  # format the validation error message
  tags$head(
    tags$style(HTML("
                    .shiny-output-error-validation {
                    color: red;
                    font-size: 17px;
                    font-style: bold;
                    }
                    "))
    ),
  
  #progress bar
  tags$head(
    tags$style(HTML("
                    .shiny-progress-container {
                    top: 20%;
                    font-size: 17px;
                    text-align: left;
                    font-style: bold;
                    }"))),
  
  #progress bar
  tags$head(
    tags$style(HTML("
                    .shiny-progress .progress-text {
                    left: 34%;
                    background-color: white;
                    width = 500px; 
                    }"))),
  
  tags$head(
    tags$style(HTML("
                    .shiny-progress .progress {
                    height: 0px;
                    }"))),

  # keywords and description
  tags$head(
    tags$meta(name = "keywords", content = "HERC, SHARP, kidney, model, CKD model, chronic kidney disease model, lifetime, long-term, longterm, Iryna Schlackow, Borislava Mihaylova")),

  tags$head(
    tags$meta(name = "description", content = "The SHARP CKD-CVD outcomes model simulates long-term cardiovascular event rates, kidney disease progression, (quality-of-life adjusted) survival and healthcare costs associated with individual patient profiles and treatments. It can be applied to patient populations with moderate-to-severe chronic kidney disease who are over 40 years of age, and can be used with individual patients as well as groups of patients. The model was developed by Iryna Schlackow and Borislava Mihaylova, based at Health Economics Research Centre (HERC), Nuffield Department of Population Health, University of Oxford")),


  
#  titlePanel("SHARP CKD-CVD outcomes model (beta version)"),
  titlePanel("SHARP CKD-CVD outcomes model"),
  sidebarLayout(
    navlistPanel(id = "panels", "Introduction", 
                 tabPanel("Model overview",
  
                          ################################################################################
                          ### Model overview
                          ################################################################################
                          
			h3(span("Please note that the model is no longer maintaned. The model code has been made available ", 
			HTML("<a href='http://tiny.cc/sharp-outcomes-model'> online </a>"), style = "color:red")),
                          h3("Introduction"),
                          p("The SHARP CKD-CVD outcomes model simulates long-term cardiovascular event rates, 
                            kidney disease progression, (quality-of-life adjusted) survival and healthcare costs 
                            associated with individual patient profiles and treatments. 
                            It can be applied to patient populations with moderate-to-severe chronic kidney disease
                            who are over 40 years of age, and can be used with individual patients as well as groups of patients."),
                          p("The model reports long-term projections as well as cost-effectiveness results 
                            comparing against the 'no treatment' strategy. 
                            The evaluated health outcomes and costs are reported separately for each treatment arm. 
                            The user can vary parameters to assess sensitivity of the results."),
                          p("To perform the analysis, specify the required parameters
                            using the 'Model parameter' tabs and click on the 'Run analyses' button
                            on the ", actionLink("link_to_Results", "Results"), " tab. 
                            Please refer to the ", HTML("<a href='userguide.pdf'> User guide</a>"), 
			    "and the ", HTML("<a href='http://heart.bmj.com/content/early/2017/08/03/heartjnl-2016-310970'> published manuscript</a>"), 
                            "for further information."),
                          p("The ", actionLink("link_to_Glossary", "Glossary"),  
                            " tab contains a list of commonly used definitions."),
                          h3("Citation"),
                          p("When referring to this program in publications, please cite the following references:"),
			    p("Schlackow I, Kent S, Herrington W, Emberson J, Haynes R, 
                            Reith C, Wanner C, Fellstr\u00F6m B, Gray A, Landray MJ, Baigent C, 
                            Mihaylova B, on behalf of the SHARP Collaborative Group.",  
                            em("A policy model of cardiovascular disease in moderate-to-advanced chronic kidney disease."), 
			    "Heart. Published Online First: 05 August 2017. doi: 10.1136/heartjnl-2016-310970"),
                          p("Schlackow I, Mihaylova B.", em("The SHARP outcomes CKD-CVD outcomes model."), 
                            "2016; available at http://dismod.ndph.ox.ac.uk/kidneymodel/app/"),
                          h3("Contact"),
                          p("For queries, bug reports and suggestions, please email ", 
                            HTML("<a href='mailto:kidneymodel@ndph.ox.ac.uk?'>
                                 kidneymodel@ndph.ox.ac.uk</a>")),
                          h3("Acknowledgements"),
                          p("We thank Oliver Verran and Seamus Kent for their contribution 
                            to the development of the first version of the model and providing further feedback. 
                            We are also grateful to the IT team of the Oxford University's 
                            Nuffield Department of Population Health for their support in installing and running the software."),
                          h3("Disclaimer"),
                          p("The web interface for the SHARP CKD-CVD outcomes model is freely available for use."),
                          p("The University of Oxford is a charitable foundation devoted to education and research, 
                            and in order to protect its assets for the benefit of those objects, 
                            the University must make it clear that no condition is made or to be implied 
                            nor is any warranty given or to be implied as to the quality or accuracy of the Tool, 
                            or that it will be suitable for any particular purpose or for use under any specific conditions. 
                            The University and its staff accept no responsibility for the use which you make of the Tool. 
                            The University's liability for anything arising out of or in connection with the Tool supplied 
                            will not extend to loss of business or profit, or to any indirect or consequential damages or losses.")
                          ),
                 tabPanel("Glossary",
                          
                          ################################################################################
                          ### Glossary
                          ################################################################################
                          
                          p("This page contains a list of specialist terms that are used throughout the interface.
                            A copy of the Glossary is also available", 
                            HTML("<a href='glossary.pdf'> here</a>")),
                          
                          p(span("Body-mass index (BMI):", style = "color:rgb(0, 102, 213)"),
                            "the weight of an individual divided by their height squared, measured in 
                            kg/m\u00B2."),
                          
                          p(span("Chronic kidney disease (CKD):", style = "color:rgb(0, 102, 213)"),
                            "a long-term condition characterised by an impaired kidney function.
                            The diagnosis is usually based on estimating or measuring patient's 
                            glomerular filtration rate (GFR) at least twice 90 days apart, with CKD defined
                            as (e)GFR <90 ml/min/1.73m\u00B2."),
                          
                          p(span("CKD stage 3B:", style = "color:rgb(0, 102, 213)"),
                            "mild-to-moderate chronic kidney disease, defined as eGFR 30-45 ml/min/1.73m\u00B2."),
                          
                          p(span("CKD stage 4:", style = "color:rgb(0, 102, 213)"),
                            "Moderate chronic kidney disease, defined as eGFR 15-29 ml/min/1.73m\u00B2."),
                          
                          p(span("CKD stage 5:", style = "color:rgb(0, 102, 213)"),
                            "Advanced chronic kidney disease, defined as eGFR <15 ml/min/1.73m\u00B2;
                            not on renal replacement therapy."),
                          
                          p(span("Compliance:", style = "color:rgb(0, 102, 213)"),
                            "the degree to which a patient takes their medication, 
                            expressed as a percentage of the time for which the patient is compliant."),    
                          
                          p(span("Cost-effectiveness acceptability curve (CEAC):", style = "color:rgb(0, 102, 213)"),
                            "a summary of the uncertainty around a cost-effectiveness estimate. 
                            It is derived from probabilistic analysis and presents the probability of the intervention 
                            being cost-effective across a range of threshold values of cost-effectiveness 
                            (also known as maximum willingness to pay for a unit of benefit, 
                            decision-maker's willingness to pay)."),
                          
                          p(span("Cost-effectiveness analysis:", style = "color:rgb(0, 102, 213)"),
                            "an economic analysis that calculates the additional/incremental costs required to realize 
                            a unit of additional benefits when comparing two interventions."),
                          
                          p(span("Deterministic analysis:", style = "color:rgb(0, 102, 213)"),
                            "results derived using the mean estimates of contributing parameters and reporting only 
                            mean estimates of results without allowing for parameter uncertainty."),
                          
                          p(span("Incremental cost-effectiveness ratio (ICER):", style = "color:rgb(0, 102, 213)"),
                            "a statistic produced by the cost-effectiveness analysis, equal to the 
                            ratio of the cost difference between two interventions 
                            and their effect difference."),
                          
                          p(span("Life-year (LY):", style = "color:rgb(0, 102, 213)"),
                            "a normal (calendar) year."),
                          
                          p(span("Renal replacement therapy (RRT):", style = "color:rgb(0, 102, 213)"),
                            "therapy used in severe chronic kidney disease, defined as undergoing long-term dialysis or 
                            being in receipt of a kidney transplant."),
                          
                          p(span("Major atherosclerotic event (MAE):", style = "color:rgb(0, 102, 213)"),
                            "non-fatal myocardial infarction or coronary death, 
                            non-haemorrhagic stroke, or arterial revascularisation procedure 
                            excluding dialysis access procedures."),
                          
                          p(span("Major vascular event (MVE):", style = "color:rgb(0, 102, 213)"),
                            "non-fatal myocardial infarction or any cardiac death,
                            any stroke, or any arterial revascularisation procedure 
                            excluding dialysis access procedures."),
                          
                          p(span("Probabilistic analysis:", style = "color:rgb(0, 102, 213)"),
                            "analysis that takes into account uncertainty in contributing parameters"),
                          
                          p(span("Quality-adjusted life-year (QALY):", style = "color:rgb(0, 102, 213)"),
                            "a measure of health, which combine survival (ie, life-years) 
                            and health-related quality of life. For example, 1 QALY is equivalent to 
                            one year in full health."),
                          
                          p(span("The Study of Heart and Renal Protection (SHARP):", 
                                 style = "color:rgb(0, 102, 213)"),
                            "a 9,270-large multinational randomised controlled trial, which compared 
                            the use of simvastatin plus ezetimibe with placebo in participants 
                            with moderate-to-severe chronic kidney disease but no major coronary disease at recruitment."),
                          
                          p(span("Treatment effect:", style = "color:rgb(0, 102, 213)"),
                            "treatment effects in the model are presented with hazard ratios 
                            typically estimated in proportional hazards survival models."),
                          
                          p(span("Vascular death (VD):", style = "color:rgb(0, 102, 213)"),
                            "death from coronary heart disease or other cardiac disease, 
                            or from any type of stroke or other vascular causes.")),
                 
                 ################################################################################
                 ### Example files
                 ################################################################################
                 
                 tabPanel("File specifications",
                          p("The following example files are provided to help with the model use, see ", 
                            HTML("<a href='userguide.pdf'> User guide</a>"), "for detailed file descriptions."),
                          br(),
                          p("Input file with patient characteristics:", 
                            HTML("<a href = 'default_patient.csv'> default values (one patient)</a> "), "and ", 
                            HTML("<a href = 'example_input.csv'> several patients</a> ")),
                          p("Non-vascular death probabilities:", 
                            HTML("<a href = '2014_UK_CKD_NVD.csv'> 2014 UK non-vascular death probabilities </a> ")),
                          p("Output analysis files: long-term projections (deterministic)", 
                            HTML("<a href = 'example_output_LP_sum.csv'> summary </a> "), "and ",
                            HTML("<a href = 'example_output_LP_ind.csv'> patient-level </a> ")),
                          p("Output analysis files: long-term projections (probabilistic)", 
                            HTML("<a href = 'example_output_LP_PSA_sum.csv'> summary</a> "), "and ",
                            HTML("<a href = 'example_output_LP_PSA_ind.csv'> patient-level </a> ")),
                          p("Output analysis files: cost-effectiveness analysis (deterministic)", 
                            HTML("<a href = 'example_output_CE_sum.csv'> summary</a> "), "and ",
                            HTML("<a href = 'example_output_CE_ind.csv'> patient-level </a> ")),
                          p("Output analysis files: cost-effectiveness analysis (probabilistic)", 
                            HTML("<a href = 'example_output_CE_PSA_sum.csv'> summary</a> "), "and ", 
                            HTML("<a href = 'example_output_CE_PSA_ind.csv'> patient-level </a> "))
                 ),
                 
                 
                 "Model parameters",
                 tabPanel("Type of analysis",
                          
                          ################################################################################
                          ### Type of analysis
                          ################################################################################
                          
                          selectInput("anal_type", "Type of analysis", 
                                      choices = c("Long-term projections",
                                                  "Cost-effectiveness analysis"), 
                                      selected = "Long-term projections"),
                          br(), 
                          selectInput("PSA", "Include uncertainty?", 
                                      choices = c("No (deterministic analysis)", 
                                                  "Yes (probabilistic analysis)"), 
                                      selected = "Deterministic analysis"),
                          conditionalPanel(condition = "input.PSA == 'Yes (probabilistic analysis)'", 
                                           sliderInput("Nsamp", label = "Number of samples", 
                                                       min = 100, max = 1000, value = 100, step = 100),
                                           
                                           br(),
                                           p("The probabilistic sensitivity analysis is currently implemented 
                                             for treatment effects, disease rirks and hospital care costs.
                                             The default coefficients from the risk equations are derived from the SHARP data 
                                             using the bootstrap method."))
                                           ),
                 tabPanel("Patient characteristics",
                          
                          
                          ################################################################################
                          ### Patient characteristics
                          ################################################################################
                          
                          p("Select characteristics for a single patient 
                            or import a text file with these characteristics for one or more patients."),
                          checkboxInput("cb_bl", 
                                        "Import a file with patient characteristics"),
                          conditionalPanel(
                            condition = "input.cb_bl == true", 
                            fileInput('file1', ""),
                            uiOutput("validated_id_file_value")),
                          conditionalPanel(
                            condition = "input.cb_bl != true ",
                            uiOutput("validated_id_screen_value"),
                            br(),
                            actionButton("reset_input_id", "Reset inputs"), 
                            br(),
                            uiOutput("reset_id")
                          )),
                 tabPanel("Treatment parameters",
                          
                          ################################################################################
                          ### Treatment parameters
                          ################################################################################
                          
                          p("Hazard ratios should correspond to full compliance with 
                            treatment for each of the outcomes below. 
                            The rates should be on the exponential scale."), 
                          conditionalPanel(
                            condition = "input.PSA == 'No (deterministic analysis)'",
                            actionButton("reset_input_trt", "Reset inputs"),
                            uiOutput("validated_trt_value"),
                            uiOutput("reset_trt")
                          ),
                          conditionalPanel(
                            condition = "input.PSA == 'Yes (probabilistic analysis)'",
                            actionButton("reset_input_trt_PSA", "Reset inputs"),
                            uiOutput("validated_trt_value_PSA"),
                            uiOutput("reset_trt_PSA")
                          )),
                 
                 tabPanel("Annual healthcare costs", 
                          
                          ################################################################################
                          ### Annual healthcare costs
                          ################################################################################
                          
                          p("The default values are based on SHARP data and UK 2014 prices."),
                          conditionalPanel(
                            condition = "input.PSA == 'No (deterministic analysis)'",
                            actionButton("reset_input_cost", "Reset inputs"),
                            uiOutput("validated_cost_value"),
                            uiOutput('reset_cost')
                          ),
                          conditionalPanel(
                            condition = "input.PSA == 'Yes (probabilistic analysis)'",
                            actionButton("reset_input_cost_PSA", "Reset inputs"),
                            uiOutput("validated_cost_value_PSA"),
                            uiOutput('reset_cost_PSA')
                          )
                 ),
                 
                 tabPanel("Health-related quality of life",
                          
                          ################################################################################
                          ### Health-related quality of life
                          ################################################################################
                          
                          p("The default values are UK quality of life (QoL) utilities estimates derived
                            from the SHARP data."),
                          p("Baseline QoL is the quality of life utility of a 60 year old female, non-smoker, 
                            with above secondary education, 
                            with BMI 25-30 kg/m\u00B2, 
                            pre-RRT CKD and without diabetic nephropathy 
                            or vascular disease."),
                          br(),
                          actionButton("reset_input_QoL", "Reset inputs"),
                          uiOutput("validated_qol_value"),
                          br(),
                          uiOutput('reset_QoL')
                          ),
                 
                 tabPanel("Non-vascular death probabilities",
                          
                          ################################################################################
                          ### Non-vascular death probabilities
                          ################################################################################
                          
                          p("The default age and sex-specific non-vascular death probabilities 
                            were derived from the 2014 UK population data."), 
                          checkboxInput("cb_nvd", 
                                        "Import a file with non-vascular death probabilities"), 
                          br(),
                          conditionalPanel(
                            condition = "input.cb_nvd == true",
                            fileInput('file2', "Choose CSV file"),
                            uiOutput("validated_nvd_value")
                          )
                          ),
                 
                 tabPanel("Decision parameters", 
                          
                          ################################################################################
                          ### Decision parameters
                          ################################################################################
                          
                          p("The default setting for the discount rates for both costs and health outcomes 
                            is 3.5% (National Institute for Health and Care Excellence, 2013)"), 
                          br(),
                          actionButton("reset_input_dp", 
                                       "Reset inputs"),
                          uiOutput("validated_dp_value"),
                          uiOutput("reset_dp")
                          ),
                 
                 "Analyses",
                 
                 tabPanel("Results",
                          
                          ################################################################################
                          ### Model overview
                          ################################################################################
                          
                          # action button  
                          uiOutput("validated_all"),
                          
                          uiOutput("text_output_detail"),
                          
                          ########## Long-term projections ##########
                          
                          conditionalPanel(
                            condition = "input.anal_type == 'Long-term projections'",
                            
                            ### deterministic analysis
                            
                            conditionalPanel(
                              condition = "input.PSA == 'No (deterministic analysis)'", 
                              tableOutput("table_LE"),
                              br(),
                              uiOutput("link_to_download_LE")
                            ),
                            
                            ### probabilistic analysis
                            
                            conditionalPanel(
                              condition = "input.PSA == 'Yes (probabilistic analysis)'",
                              tableOutput("table_LE_PSA"),
                              br(),
                              uiOutput("link_to_download_LE_PSA")
                            )
                          ),
                          
                          ########## Cost-effectiveness analysis ##########
                          
                          conditionalPanel(
                            condition = "input.anal_type == 'Cost-effectiveness analysis'",
                            
                            uiOutput("display_disc"),
                            
                            ### deterministic analysis
                            
                            conditionalPanel(
                              condition = "input.PSA == 'No (deterministic analysis)'", 
                              tableOutput("table_Tx_C"), 
                              tableOutput("table_Tx_T"),
                              conditionalPanel(
                                condition = "input.disc == false",
                                tableOutput("table_CE_undisc")
                              ),
                              conditionalPanel(
                                condition = "input.disc == true",
                                tableOutput("table_CE_disc")
                              ), 
                              br(),
                              uiOutput("link_to_download")
                            ),
                            
                            ### probabilistic analysis
                            
                            conditionalPanel(
                              condition = "input.PSA == 'Yes (probabilistic analysis)'",
                              tableOutput("table_Tx_C_PSA"), 
                              tableOutput("table_Tx_T_PSA"),
                              conditionalPanel(
                                condition = "input.disc == false",
                                tableOutput("table_CE_undisc_PSA"),
                                br(),
                                plotOutput("p_CEAC_undisc")
                              ),
                              conditionalPanel(
                                condition = "input.disc == true",
                                tableOutput("table_CE_disc_PSA"),
                                br(),
                                plotOutput("p_CEAC_disc")
                              ), 
                              br(),
                              uiOutput("link_to_download_PSA")
                            )
                          ))
                 ),
    mainPanel())
                                           ))