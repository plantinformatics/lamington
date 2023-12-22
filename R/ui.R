# UI Code
navbarPage("Lamington",
           theme = shinythemes::shinytheme("spacelab"),
           
           tabPanel("Convert VCF File",
                    sidebarLayout(
                      sidebarPanel(width =2,
                                   radioButtons("file_location", "Select file from",
                                                choices = c(Desktop = "Client",
                                                            Server = "Server"),
                                                selected = "Client"),           
                                   
                                   textInput("GDS_filename", "GDS file name"),
                                   
                                   conditionalPanel(
                                     condition = "input.file_location == 'Client'",
                                     fileInput("file1", "Choose VCF File",
                                               multiple = FALSE,
                                               accept = c(".vcf.gz",
                                                          ".vcf",
                                                          ".csv"))
                                     ),
                                   conditionalPanel(
                                     condition = "input.file_location == 'Server'",
                                     selectInput("file2", "Select a File:", selected = NULL, choices = list.files(path ="../Data/",pattern = "\\.(vcf|vcf.gz)$")),  
                                     actionButton("refresh_list_server", "Refresh File List")
                                   ),
                                   conditionalPanel(
                                     condition = "input.GDS_filename != ''",
                                     tags$br(),
                                     actionButton("convert_button", "Convert and Display")
                                     
                                   )
                      ),
                      mainPanel(
                        verbatimTextOutput("output_text", placeholder = FALSE)%>% shinycssloaders::withSpinner(2,custom.css= TRUE)
                      )
                    )
           ),
           tabPanel("Select GDS file",
                    sidebarLayout(
                      sidebarPanel(width =2,
                                   selectInput("file_list1", "Select a GDS File:", selected = NULL, choices = list.files(path ="../Data/",pattern = "\\.gds$")),  
                                   tags$div('class'="btn-group-sm",
                                            actionButton("refresh_list", "Refresh List"),
                                            actionButton("check_file", "Check File")
                                            
                                   ),
                                   tags$br(),
                                   checkboxInput("pca_datafiles","Use the data downloaded from the specified GDS file",value = FALSE),
                                   conditionalPanel(
                                     condition = "input.pca_datafiles == '1'",
                                     tags$hr(),
                                     fileInput("data_file", "Upload a CSV file for PCA",
                                               multiple = FALSE,
                                               accept = c("text/csv",
                                                          "text/comma-separated-values,text/plain",
                                                          ".csv")),
                                     fileInput("data_file2", "Upload a RDS file with Snpsetid",
                                               multiple = FALSE,
                                               accept = c(".rds"))
                                   )
                      ),
                      mainPanel(
                        verbatimTextOutput("file_summary"),
                        tableOutput("tab2_df")%>% shinycssloaders::withSpinner(2,custom.css= TRUE)
                      )
                    )
           ),
           tabPanel("Histrogram",
                    sidebarLayout(
                      sidebarPanel(width =2,
                                   sliderInput("n", "Slide to change missing rate:",
                                               min = 0, max = 1, value = 0.3, step = 0.01),
                                   br(),
                                   sliderInput("maf", "Slide to change MAF:",
                                               min = 0, max = 0.5, value = 0.2, step = 0.001),
                                   br(),
                                   tags$label("GDS file:"),
                                   verbatimTextOutput("filename_histro",placeholder = TRUE),
                                   tags$label("Total Samples:"),
                                   verbatimTextOutput("GDS_sample",placeholder = TRUE),
                                   tags$label("Total SNPS:"),
                                   verbatimTextOutput("GDS_snps",placeholder = TRUE)
                                   
                      ),
                      mainPanel(
                        plotOutput(outputId = "histro")%>% shinycssloaders::withSpinner(2,custom.css= TRUE),
                        tags$hr(),
                        plotOutput(outputId = "histro2") %>% shinycssloaders::withSpinner(2,custom.css= TRUE)
                      )
                    )
           ),
           tabPanel("Genotype Matrix",
                    sidebarLayout(
                      sidebarPanel(width =2,
                                   conditionalPanel(
                                     condition = "input.pca_datafiles == '0'",
                                     radioButtons("pruning", "Select Data:",
                                                  choices = c(Complete = "complete",
                                                              "Filter Data" = "sample"),
                                                  selected = "complete"),
                                     conditionalPanel(
                                       condition = "input.pruning == 'sample'",
                                       numericInput("sample_size","Enter Sample Size",value= 10,min=1),
                                       numericInput("thres_size","Enter Threshold Value for LD - Pruning",value= 0.0001,min=0.001, max = 1,step = 0.001),
                                       numericInput("missing_rate","Missing Rate",value= 0.5,min=0.001, max = 1,step = 0.001),
                                       numericInput("maf_rate","MAF",value= 0.5,min=0.001, max = 0.5,step = 0.001),
                                       #selectInput("autosome", "Use autosomal SNPs", selected = "TRUE", choices = c("TRUE"="TRUE","FALSE"="FALSE"))
                                       checkboxInput("autosome", "Autosome.only", value = FALSE),
                                     ),
                                     conditionalPanel(
                                       condition = "input.pruning == 'complete'",
                                       p("Using all Samples and SNP Data"),
                                       tags$label("Convert GDS File")      
                                     ),
                                     actionButton("get_geno", "Get Geno")
                                   ),
                                   ######################### Regenerate ############################
                                   tags$hr(),
                                   tags$label("Regenerate the Genotype Matrix"), 
                                   actionButton("regeno", "Regenerate Geno"),
                                   tags$br(),
                                   p("Regenerate with updated results from PCA Graph or using previously downloaded results"),
                                   tags$br(),
                                   tags$label("Deleted Samples"),
                                   verbatimTextOutput("del_samples",placeholder = TRUE)
                                   
                                   ##################################################################
                      ),
                      mainPanel(
                        verbatimTextOutput("gds_summary") %>% shinycssloaders::withSpinner(2,custom.css= TRUE)
                      )
                    )
           ),
           tabPanel("PCA",
                    sidebarLayout(
                      sidebarPanel(width =2,
                                   tags$label("Principal Component Analysis (PCA)"),
                                   #selectInput("autosome", "Use autosomal SNPs", selected = "TRUE", choices = c("TRUE"="TRUE","FALSE"="FALSE")),  
                                   verbatimTextOutput("autosome_out"),
                                   numericInput("PCA_Thread","Enter CPU Cores",value= 4,min=1),
                                   actionButton("get_PCA", "Start PCA")
                      ),
                      mainPanel(
                        verbatimTextOutput("pca_summary") %>% shinycssloaders::withSpinner(2,custom.css= TRUE)
                      )
                    )
           ),
           tabPanel("CoreHunter",
                    sidebarLayout(
                      sidebarPanel(width =2,
                                   tags$label("Select Core Hunter Cores"), 
                                   checkboxGroupInput(inputId = "id_check_2", label = NULL, selected = 3,
                                                      choices = c("Core10" = 10, "Core20" = 20, "Core50" = 50, "Core100" = 100,  "Core200" = 200)),
                                   textInput("core_values", "Enter Cores",value = "", placeholder ="10,15,16"),
                                   checkboxGroupInput("select_cores", "Selected Cores",
                                                      choices = c("Always Selected", "Never Selected"),
                                                      selected = character(0)),
                                   conditionalPanel(
                                     condition = "input.select_cores.includes('Always Selected')",
                                     textInput("always_cores", "Always Selected Cores",placeholder = 'sample1 sample2 sample3')
                                   ),
                                   
                                   conditionalPanel(
                                     condition = "input.select_cores.includes('Never Selected')",
                                     textInput("never_cores", "Never Selected Cores",placeholder = 'sample1 sample2 sample3')
                                   ),
                                   checkboxInput("core_adv", "Add Objectives (Advance)", FALSE),
                                   conditionalPanel(
                                     condition = "input.core_adv== '1'",
                                     selectInput(inputId = "obj_type", label = "Type", choices = c("EN (Default)" = "EN", "AN", "EE", "SH", "HE", "CV"), selected = "EN (Default)"),
                                     selectInput(inputId = "obj_measure", label = "Measure", choices = c("MR (Default)"= "MR", "CE"), selected = "MR (Default)"),
                                     numericInput("obj_weight", "Weight:", 1, min = 0.1, max = 1,step= 0.1)
                                     
                                   ),
                                   tags$p("This may take a long time"), 
                                   actionButton("get_core", "Run Corehunter"),
                                   tags$p("Note: Before running CoreHunter, increase memory in Settings if needed otherwise default memory will be selected")
                      ),
                      mainPanel(
                        verbatimTextOutput("hunter_summary") %>% shinycssloaders::withSpinner(2,custom.css= TRUE)
                      )
                    )
           ),
           tabPanel("Add POP Data",
                    sidebarLayout(
                      sidebarPanel(width =2,
                                   # Input: Select a file ----
                                   tags$label("Add Population Data CSV"),
                                   fileInput("csv_file", NULL,
                                             multiple = FALSE,
                                             accept = c("csv",
                                                        ".csv")),
                                   checkboxInput("pop_adv", "Options", FALSE),
                                   conditionalPanel(
                                     condition = "input.pop_adv == '1'",
                                     tags$hr(),
                                     checkboxInput("header", "Header", TRUE),
                                     radioButtons("sep", "Separator",
                                                  choices = c(Comma = ",",
                                                              Semicolon = ";",
                                                              Tab = "\t"),
                                                  selected = ","),
                                     radioButtons("quote", "Quote",
                                                  choices = c(None = "",
                                                              "Double Quote" = '"',
                                                              "Single Quote" = "'"),
                                                  selected = '"'),
                                     tags$hr(),
                                     radioButtons("disp", "Display",
                                                  choices = c(Head = "head",
                                                              All = "all"),
                                                  selected = "head")
                                   ),
                                   conditionalPanel(
                                     condition = "input.sample_id !== ''",
                                     tags$hr(),
                                     selectInput(
                                       "data_primarykey", "Select primary key present in data to join ",choices = NULL, selected = NULL, multiple = FALSE
                                     ),
                                     tags$hr(),
                                     selectInput(
                                       "sample_id", "Primary Key or Sample ID",choices = NULL, selected = NULL, multiple = FALSE
                                     ),
                                     selectInput(
                                       "population_key", "Select population data",choices = NULL, selected = NULL, multiple = TRUE
                                     ),
                                     # Only show this panel if Custom is selected
                                     conditionalPanel(
                                       condition = "input.sample_id !== '' && input.population_key !== ''&& input.data_primarykey !== ''",
                                       actionButton("add_toDF", "Add columns to data")
                                     )
                                   ),
                                   tags$hr(),
                                   tags$label("Show PCA DataFrame"),
                                   tags$br(),
                                   actionButton("show_df", "PCA Data")
                      ),
                      mainPanel(
                        tableOutput("contents")%>% shinycssloaders::withSpinner(2,custom.css= TRUE),
                        tags$hr(),
                        tableOutput("final_df")%>% shinycssloaders::withSpinner(2,custom.css= TRUE)
                      )
                    )
           ),
           tabPanel("PCA Plot",
                    sidebarLayout(
                      sidebarPanel(width =2,
                                   selectInput(inputId = "x", label = "X Axis", choices = c("EV1","EV2","EV3","EV4"), selected = "EV1"),
                                   selectInput(inputId = "y", label = "Y Axis", choices = c("EV1","EV2","EV3","EV4"), selected = "EV2"),
                                   selectInput(inputId = "title", label = "Type", choices = NULL),
                                   actionButton("refresh_selection", "Refresh Selection"),
                                   tags$br(),
                                   tags$br(),
                                   checkboxInput("ellipse", "Confidence ellipses", value = FALSE),
                                   sliderInput("opacity", "Points opacity :", min = 0, max = 1, value = 1, step = 0.05),
                                   sliderInput("point_size", "Points size :", min = 0, max = 100, value = 50, step = 1),
                                   tags$div('class'= "btn-group-sm",
                                            actionButton("scatterD3-reset-zoom", HTML("<span class='glyphicon glyphicon-search' aria-hidden='true'></span> Reset Zoom"),class="mb-2"),
                                            actionButton("scatterD3-lasso-toggle", HTML("<span class='glyphicon glyphicon-screenshot' aria-hidden='true'></span> Toggle Lasso"), "data-toggle" = "button")
                                   ),
                                   tags$hr(),
                                   actionButton("delete_rows", "Deleted selected rows",class= "btn-sm"),
                                   tags$hr(),
                                   h6("Download"),
                                   tags$div('class'= "btn-group-sm",
                                            downloadButton("downloadData", "PCA"),
                                            downloadButton("downloadData2", "SNPSet ID")
                                   )
                      ),
                      mainPanel(
                        fluidRow(
                          shinycssloaders::withSpinner(scatterD3Output("scatterPlot", height = "700px"),2,custom.css= TRUE),
                          wellPanel(
                            tags$label("Selected Samples:"),
                            verbatimTextOutput("sel_samples",placeholder = TRUE), 
                            tags$label("Deleted Samples:"),
                            verbatimTextOutput("brush_info",placeholder = TRUE)
                          ),
                          dataTableOutput("table")
                        )
                      )
                    )
           ),
           tabPanel("Final Plot",
                        tabsetPanel(
                          tabPanel(title = "Final Plot Download",
                                   esquisse_ui(id = "esquisse",header = FALSE ),
                                   tags$hr(),
                                   actionButton("Refresh_p","Refresh Plot",class = "btn btn-default btn-lg btn-block")
                          ),
                          tabPanel( title = "Output",
                                    tags$b("Code:"),
                                    verbatimTextOutput("code"),
                                    tags$b("Filters:"),
                                    verbatimTextOutput("filters"),
                                    tags$b("Data:"),
                                    verbatimTextOutput("data")
                         ))
           ),
           tabPanel("Settings",
                    sidebarLayout(
                      sidebarPanel(width =2,
                                   numericInput("Corehunter_Me","Set CoreHunter Memory - Default is 0.5G",value= 10,min=1,step = 1),
                                   tags$div('class'="btn-group-vertical btn-group-sm",
                                            actionButton("check_memory", "Check heap memory"),
                                            br(),
                                            actionButton("set_mem", "Set Heap memory"))
                      ),
                      mainPanel(
                        verbatimTextOutput("java_mem")
                      )
                    )
           )
           
)