# UI Code
navbarPage(
  "Lamington",
  id = "mainnavbar",
  windowTitle = "Lamington|PCA and Core Diversity Set Explorer",
  theme = shinythemes::shinytheme("spacelab"),
  footer = uiOutput("testname"),
  ####========Tab for login and user rego===========================####
  tabPanel("User", value = "usertab", sidebarLayout(
    sidebarPanel(
      width = 2,
      wellPanel(
        id = "userContent",
        textInput("username", "Username"),
        passwordInput("password", "Password", ),
        conditionalPanel(
          condition = "output.loginoption=='login'",
          actionButton(
            "login",
            HTML("<span class='glyphicon glyphicon-log-in'></span> Login")
          ),
          tags$hr(),
          actionLink(
            "signup",
            HTML("or <span class='glyphicon glyphicon-user'></span> Sign Up here")
          )
        ),
        uiOutput("signupoptions"),
        ##Include option to listen to return/enter key event in password textfield
        tags$script(
          '$(document).on("keyup", "#password", function(e) {
            if(e.keyCode == 13) {
              Shiny.onInputChange("passwordEntered", Math.random());
            }
          });
        '
        )
      )
    ), mainPanel(textOutput("userlog"))
  )),
  ####========Tab for adding population data=======================####
  tabPanel("Add POP Data", id = "popdata_tab", sidebarLayout(
    sidebarPanel(
      width = 2,
      # Input: Select a file ----
      tags$label("Add Population Data"),
      fileInput(
        "csv_file",
        NULL,
        placeholder = 'select a csv file',
        multiple = FALSE,
        accept = c("csv", ".csv")
      ),
      checkboxInput("pop_adv", "File options", FALSE),
      tags$hr(),
      #-----------Add options for table ---------------
      conditionalPanel(
        condition = "input.pop_adv == '1'",
        tags$hr(),
        checkboxInput("header", "Header", TRUE),
        radioButtons(
          "sep",
          "Separator",
          choices = c(
            Comma = ",",
            Semicolon = ";",
            Tab = "\t"
          ),
          selected = ","
        ),
        radioButtons(
          "quote",
          "Quote",
          choices = c(
            None = "",
            "Double Quote" = '"',
            "Single Quote" = "'"
          ),
          selected = '"'
        )
      ),
      #-------Option to add more columns to the data frame-------
      conditionalPanel(
        condition = "output.csvtest",
        tags$h4("Update population data"),
        textInput("column_name", "New column name:", value = ""),
        actionButton(
          "add_toDF",
          HTML(
            "<span class='glyphicon glyphicon-edit'></span> Add new column to data"
          )
        ),
        tags$hr()
      ),
      
      conditionalPanel(
        condition = "input.add_toDF > 0 & input.column_name!=''",
        radioButtons(
          "updateopt",
          "Option to update table",
          choices = c(
            "Select rows in table" = 'ids_selected',
            "Paste list of IDs to update" = "ids_pasted"
          ),
          selected = 'ids_selected'
        ),
      ),
      conditionalPanel(
        condition = "input.updateopt == 'ids_pasted'",
        # Text input for new column values
        textAreaInput("pastedIds", "Enter IDs as single lines:", rows = 5)
      ),
      uiOutput("selsampleID"),
      uiOutput("popColnameText"),
      uiOutput("popColnameBtn"),
      tags$hr()
    ),
    mainPanel(
      DTOutput("pop_table") %>% shinycssloaders::withSpinner(2, custom.css = TRUE),
      # Button to clear selected rows (conditionally displayed)
      uiOutput("clearButton"),
      tags$hr(),
      #textOutput("csvtest"),
      #textOutput("pop_colname")
      #tableOutput("contents") %>% shinycssloaders::withSpinner(2, custom.css = TRUE)
    )
  )),
  ##-------Tab panel for uploading VCFs and converting to GDS format---
  tabPanel("Convert VCF File", sidebarLayout(
    sidebarPanel(
      width = 2,
      radioButtons(
        "vcf_location",
        "Select file from",
        choices = c(Desktop = "Client", Server = "Server"),
        selected = "Client"
      ),
      
      textInput("GDS_filename", "GDS file name"),
      
      conditionalPanel(
        condition = "input.vcf_location == 'Client'",
        fileInput(
          "vcf_desktop",
          "Choose VCF File",
          multiple = FALSE,
          accept = c(".vcf.gz", ".vcf", ".csv")
        )
      ),
      conditionalPanel(
        condition = "input.vcf_location == 'Server'",
        selectInput(
          "vcf_server",
          "Select a File:",
          selected = NULL,
          choices = list.files(path = "../VCFs/", pattern = "\\.(vcf|vcf.gz)$")
        ),
        actionButton(
          "refresh_list_server",
          HTML(
            "<span class='glyphicon glyphicon-refresh'></span> Refresh List"
          )
        )
      ),
      conditionalPanel(
        condition = "input.GDS_filename != ''",
        tags$br(),
        actionButton("convert_button", "Convert and Display")
        
      )
    ),
    mainPanel(
      verbatimTextOutput("output_text", placeholder = FALSE) %>% shinycssloaders::withSpinner(2, custom.css = TRUE)
    )
  )),
  ##-------Tab panel for selecting GDS file----------------------------
  tabPanel("Select GDS file", sidebarLayout(
    sidebarPanel(
      width = 2,
      selectInput(
        "gds_filelist",
        "Select a GDS File:",
        selected = NULL,
        choices = NULL
      ),
      tags$div(
        'class' = "btn-group-sm",
        actionButton(
          "refresh_gds_list",
          HTML(
            "<span class='glyphicon glyphicon-refresh'></span> Refresh List"
          )
        ),
        actionButton(
          "check_gdsfile",
          HTML("<span class='glyphicon glyphicon-check'></span> Check File")
        )
        
      ),
      uiOutput("deleteGDSButton")
    ),
    mainPanel(
      verbatimTextOutput("file_summary") %>% shinycssloaders::withSpinner(2, custom.css = TRUE),
      tableOutput("tab2_df") %>% shinycssloaders::withSpinner(2, custom.css = TRUE)
    )
  )),
  ##--------------Tab panel for plotting histogram---------------------
  tabPanel("Histogram", sidebarLayout(
    sidebarPanel(
      width = 2,
      sliderInput(
        "n",
        "Slide to change missing rate:",
        min = 0,
        max = 1,
        value = 0.3,
        step = 0.1
      ),
      br(),
      sliderInput(
        "maf",
        "Slide to change MAF:",
        min = 0,
        max = 0.5,
        value = 0.2,
        step = 0.01
      ),
      br(),
      tags$label("GDS file:"),
      verbatimTextOutput("filename_histo", placeholder = TRUE),
      tags$label("Total Samples:"),
      verbatimTextOutput("GDS_sample", placeholder = TRUE),
      tags$label("Total SNPS:"),
      verbatimTextOutput("GDS_snps", placeholder = TRUE),
      tags$hr(),
      uiOutput("pophistoColname"),
      uiOutput("pophistogroups"),
      uiOutput("pophistosample"),
      uiOutput("pophistoupdateBtn")
    ),
    mainPanel(
      #fluidRow(
      plotOutput(outputId = "histo") %>% shinycssloaders::withSpinner(2, custom.css = TRUE),
      tags$hr(),
      plotOutput(outputId = "histo2") %>% shinycssloaders::withSpinner(2, custom.css = TRUE)
      #)
    )
  )),
  ##-------------Tab panel for Genotype matrix---------------------------------------------
  tabPanel("Genotype Matrix", sidebarLayout(
    sidebarPanel(
      width = 2,
      radioButtons(
        "geno_filter",
        "Select Data:",
        choices = c(
          "No Filter" = "complete",
          "Filtering options" = "sample",
          "Regenerate from existing data" = "regenerate"
        ),
        selected = "complete"
      ),
      conditionalPanel(
        condition = "input.geno_filter == 'sample'",
        numericInput("sample_size", "Sample Size", value = 10, min = 1),
        numericInput(
          "ld_cutoff",
          "LD - Pruning",
          value = 0.9,
          min = 0,
          max = 1,
          step = 0.1
        ),
        numericInput(
          "missing_rate",
          "Missing Rate",
          value = 0.5,
          min = 0.1,
          max = 1,
          step = 0.1
        ),
        numericInput(
          "maf_rate",
          "MAF",
          value = 0.5,
          min = 0.001,
          max = 0.5,
          step = 0.01
        ),
        checkboxInput("filtsamples", "Filter samples", value = FALSE),
        uiOutput("txtsamplefilt"),
        conditionalPanel(
          condition = "output.filtsampselected",
          uiOutput("radioSampleIdfilt"),
          uiOutput("selectsampleIdcol"),
          uiOutput("selectsamplegroups"),
          uiOutput("selectgrouptoretain")
        ),
        tags$hr(),
        tags$script(
          '$(document).on("click", "#samplestofiltlist", function(e) {
              Shiny.onInputChange("samplelistEntered", Math.random());
          });
        '
        ),
        tags$script('
          $(document).on("change", "input[type=radio][name=selectpopfiltOption]", function() {
            Shiny.onInputChange("selectpopfiltOptionEntered", this.value);
          });
        ')
        
      ),
      #selectInput("autosome", "Use autosomal SNPs", selected = "TRUE", choices = c("TRUE"="TRUE","FALSE"="FALSE"))
      checkboxInput("autosome", "Autosome only", value = FALSE),
      conditionalPanel(
        condition = "input.geno_filter == 'complete'",
        p("Using all Samples and SNP Data"),
        tags$label("Convert GDS File")
      ),
      ######################### Regenerate ############################
      conditionalPanel(
        condition = "input.geno_filter == 'regenerate'",
        p(
          "Regenerate with updated results from PCA Graph or using previously saved results"
        ),
        conditionalPanel(
          condition = "output.pca_status",
          tags$br(),
          tags$label("Show PCA DataFrame"),
          tags$br(),
          actionButton("show_df", "PCA Data"),
          tags$br(),
          tags$label("Deleted Samples"),
          verbatimTextOutput("del_samples", placeholder = TRUE)
        ),
        tags$hr(),
        checkboxInput(
          "pca_datafiles",
          "Use data from saved/existing PCA and GDS files",
          value = FALSE
        ),
        #---Options to add extra information to
        conditionalPanel(
          condition = "input.pca_datafiles == '1'",
          tags$hr(),
          fileInput(
            "data_file",
            "Upload CSV file with PCA info",
            multiple = FALSE,
            accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
          ),
          fileInput(
            "data_file2",
            "Upload RDS file with Snpsetid",
            multiple = FALSE,
            accept = c(".rds")
          )
        ),
        tags$hr(),
        #tags$label("Regenerate the Genotype Matrix"),
        #actionButton("regeno", "Regenerate Geno"),
        #tags$br(),
        
      ),
      actionButton(
        "get_geno",
        HTML("<span class='glyphicon glyphicon-import'></span> Get Geno")
      ),
    ),
    mainPanel(
      verbatimTextOutput("gds_summary") %>% shinycssloaders::withSpinner(2, custom.css = TRUE),
      tags$br(),
      DTOutput("final_df") %>% shinycssloaders::withSpinner(2, custom.css = TRUE),
      tags$br(),
      #verbatimTextOutput("pca_status",placeholder = T)
    )
  )),
  ############################ PCA Tab ###############################################################
  tabPanel("PCA", sidebarLayout(
    sidebarPanel(
      width = 2,
      tags$label("Principal Component Analysis (PCA)"),
      #selectInput("autosome", "Use autosomal SNPs", selected = "TRUE", choices = c("TRUE"="TRUE","FALSE"="FALSE")),
      uiOutput("show_pcaoptions"),
      uiOutput("show_loadcalcpca"),
      #-----Option to map pop data to PCA data
      conditionalPanel(
        condition = "output.csvtest && input.get_PCA>0",
        tags$hr(),
        tags$h4("Map population data to genotype data"),
        selectInput(
          "data_primarykey",
          "Primary key in data to join ",
          choices = NULL,
          selected = NULL,
          multiple = FALSE
        ),
        tags$hr(),
        selectInput(
          "sample_id",
          "Primary Key or Sample ID",
          choices = NULL,
          selected = NULL,
          multiple = FALSE
        ),
        selectInput(
          "population_key",
          "Select population data",
          choices = NULL,
          selected = NULL,
          multiple = TRUE
        ),
        tags$hr(),
        actionButton(
          "add_toPCA",
          HTML(
            "<span class='glyphicon glyphicon-saved'></span> Update PCA data"
          )
        ),
        tags$hr()
        
      )
      
    ),
    mainPanel(
      verbatimTextOutput("pca_summary") %>% shinycssloaders::withSpinner(2, custom.css = TRUE),
      DTOutput("pca_dt") %>% shinycssloaders::withSpinner(2, custom.css = TRUE),
      
    )
  )),
  ########################## Core Hunter Tab ###########################################################
  tabPanel("CoreHunter", sidebarLayout(
    sidebarPanel(
      width = 2,
      tags$label("Select Core Hunter Cores"),
      checkboxGroupInput(
        inputId = "id_check_2",
        label = NULL,
        selected = 3,
        choices = c(
          "Core10" = 10,
          "Core20" = 20,
          "Core50" = 50,
          "Core100" = 100,
          "Core200" = 200
        )
      ),
      textInput(
        "core_values",
        "Enter Cores",
        value = "",
        placeholder = "10,15,16"
      ),
      checkboxGroupInput(
        "select_cores",
        "Selected Cores",
        choices = c("Always Selected", "Never Selected"),
        selected = character(0)
      ),
      conditionalPanel(
        condition = "input.select_cores.includes('Always Selected')",
        textInput("always_cores", "Always Selected Cores", placeholder = 'sample1 sample2 sample3')
      ),
      
      conditionalPanel(
        condition = "input.select_cores.includes('Never Selected')",
        textInput("never_cores", "Never Selected Cores", placeholder = 'sample1 sample2 sample3')
      ),
      checkboxInput("core_adv", "Add Objectives (Advance)", FALSE),
      conditionalPanel(
        condition = "input.core_adv== '1'",
        selectInput(
          inputId = "obj_type",
          label = "Type",
          choices = c("EN (Default)" = "EN", "AN", "EE", "SH", "HE", "CV"),
          selected = "EN (Default)"
        ),
        selectInput(
          inputId = "obj_measure",
          label = "Measure",
          choices = c("MR (Default)" = "MR", "CE"),
          selected = "MR (Default)"
        ),
        numericInput(
          "obj_weight",
          "Weight:",
          1,
          min = 0.1,
          max = 1,
          step = 0.1
        )
        
      ),
      tags$p("This may take a long time"),
      actionButton(
        "get_core",
        HTML(
          "<span class='glyphicon glyphicon-play-circle'></span> Run Corehunter"
        )
      ),
      tags$p(
        "Note: Before running CoreHunter, increase memory in Settings if needed otherwise default memory will be selected"
      )
    ),
    mainPanel(
      verbatimTextOutput("hunter_summary") %>% shinycssloaders::withSpinner(2, custom.css = TRUE)
    )
  )),
  
  ####========Tab for generating PCA Plot=======================####
  tabPanel("PCA Plot", sidebarLayout(
    sidebarPanel(
      width = 2,
      selectInput(
        inputId = "x",
        label = "X Axis",
        choices = c("EV1", "EV2", "EV3", "EV4"),
        selected = "EV1"
      ),
      selectInput(
        inputId = "y",
        label = "Y Axis",
        choices = c("EV1", "EV2", "EV3", "EV4"),
        selected = "EV2"
      ),
      selectInput(
        inputId = "title",
        label = "Type",
        choices = NULL
      ),
      actionButton(
        "refresh_selection",
        HTML(
          "<span class='glyphicon glyphicon-refresh'></span> Refresh Selection"
        ),
        class = "mb-2"
      ),
      tags$br(),
      tags$br(),
      checkboxInput("ellipse", "Confidence ellipses", value = FALSE),
      sliderInput(
        "opacity",
        "Points opacity :",
        min = 0,
        max = 1,
        value = 1,
        step = 0.05
      ),
      sliderInput(
        "point_size",
        "Points size :",
        min = 0,
        max = 100,
        value = 50,
        step = 1
      ),
      tags$div(
        'class' = "btn-group-sm",
        actionButton(
          "scatterD3-reset-zoom",
          HTML(
            "<span class='glyphicon glyphicon-zoom-out' aria-hidden='true'></span> Reset Zoom"
          ),
          class = "mb-2"
        ),
        actionButton(
          "scatterD3-lasso-toggle",
          HTML(
            "<span class='glyphicon glyphicon-screenshot' aria-hidden='true'></span> Toggle Lasso"
          ),
          "data-toggle" = "button"
        ),
        tags$hr(),
        actionButton(
          "reset_selected_rows",
          HTML(
            "<span class='glyphicon glyphicon-repeat'></span> Reset selection"
          ),
          class = "btn-sm"
        )
      ),
      tags$hr(),
      actionButton(
        "delete_rows",
        HTML(
          "<span class='glyphicon glyphicon-remove'></span> Delete selected rows"
        )
        ,
        class = "btn-sm"
      ),
      
      tags$hr(),
      h6("Download"),
      tags$div(
        'class' = "btn-group-sm",
        downloadButton("downloadData", "PCA"),
        downloadButton("downloadData2", "SNPSet ID")
      )
    ),
    mainPanel(
      fluidRow(
        shinycssloaders::withSpinner(
          scatterD3Output("scatterPlot", height = "700px"),
          2,
          custom.css = TRUE
        ),
        wellPanel(
          tags$label("Selected Samples:"),
          verbatimTextOutput("sel_samples", placeholder = TRUE),
          tags$label("Deleted Samples:"),
          verbatimTextOutput("brush_info", placeholder = TRUE)
        ),
        dataTableOutput("table")
      )
    )
  )),
  tabPanel("Final Plot", tabsetPanel(
    tabPanel(
      title = "Final Plot Download",
      esquisse_ui(id = "esquisse", header = FALSE),
      tags$hr(),
      actionButton("Refresh_p", "Refresh Plot", class = "btn btn-default btn-lg btn-block")
    ),
    tabPanel(
      title = "Output",
      tags$b("Code:"),
      verbatimTextOutput("code"),
      tags$b("Filters:"),
      verbatimTextOutput("filters"),
      tags$b("Data:"),
      verbatimTextOutput("data")
    )
  )),
  tabPanel("Settings", sidebarLayout(
    sidebarPanel(
      width = 2,
      numericInput(
        "Corehunter_Me",
        "Set CoreHunter Memory - Default is 0.5G",
        value = 10,
        min = 1,
        step = 1
      ),
      tags$div(
        'class' = "btn-group-vertical btn-group-sm",
        actionButton("check_memory", "Check heap memory"),
        br(),
        actionButton("set_mem", "Set Heap memory")
      )
    ),
    mainPanel(verbatimTextOutput("java_mem"))
  ))
)