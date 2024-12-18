# Load global variables
source("global.R")

server <- function(input, output, session) {
  #Session based variables
  histo.df <- NULL
  RV <<- NULL
  tab2 <<- NULL
  pca <<- NULL
  random_ids <<- NULL
  snpset <<- NULL
  geno <<- NULL
  snpset.id <<- NULL #Global variable that holds a set of SNP Ids
  snps <<- NULL
  deleted_samples <<- NULL
  lasso_data <<- NULL
  pop_data <<- NULL #Variable holding the population data
  pop_data_live <<- NULL
  gds_filepath <<- "" #Path to selected GDS file
  gds_data <<- NULL #GDS data
  data_path <- "../GDS/" #Path where all GDS files are saved
  vcf_path <- "../VCFs/" #Path where VCF files can be saved and uploaded directly to Lamington
  histogroupcol<<-NULL #variable to hold column in population data table
  histosamplecol<<-NULL 

  
  #Check if paths exists
  if (!dir.exists(data_path))
  {
    dir.create(data_path, recursive = T)
  }
  
  if (!dir.exists(vcf_path))
  {
    dir.create(vcf_path, recursive = T)
  }
  
  ######################### Convert vcf to GDS format ############################
  output$output_text <- renderPrint(invisible())
  
  observeEvent(input$convert_button, {
    GDS_filename <- paste(input$GDS_filename, ".gds", sep = "")
    output_gds_path <- file.path(data_path, GDS_filename)
    
    tryCatch({
      withProgress(message = 'Converting VCF file to GDS format', value = 0.1, {
        if (input$file_location == "Client") {
          if (is.null(input$vcf_desktop)) {
            show_toast(
              title = "No VCF file selected!",
              text = "Select and load a VCF file from the Desktop",
              type = "error",
              position = "center"
            )
            req(input$vcf_desktop != "")
          }
          file_select <- input$vcf_desktop$datapath
        }
        else
        {
          if (input$vcf_server == "") {
            show_toast(
              title = "No VCF file selected!",
              text = "Select and load a VCF file from the Server",
              type = "error",
              position = "center"
            )
            req(input$vcf_server != "")
          }
          file_select <- file.path(vcf_path, input$vcf_server)
        }
        incProgress(0.2, detail = "VCF File Selected")
        
        output_text <- capture.output({
          cat("Converting VCF to GDS...\n")
          incProgress(0.1, detail = "Processing..This may take a while...")
          snpgdsVCF2GDS_R(file_select, output_gds_path, method = "biallelic.only")
          cat("Conversion complete.\n")
          incProgress(0.3, detail = "Completing..")
        })
        output$output_text <- renderPrint({
          output_text
        })
      })
    }
    , error = function(e) {
      output$output_text <- renderPrint({
        return (e)
      })
      
    })
  })
  
  observeEvent(input$refresh_list_server, {
    updateSelectInput(session,
                      "file2",
                      choices = list.files(path = data_path, pattern = "\\.(vcf|vcf.gz)$"))
  })
  
  ######################### Refresh file and Check GDS file ############################
  
  output$file_summary <- renderPrint(invisible())
  
  observeEvent(c(input$refresh_list, input$convert_button), {
    RV <<- NULL
    updateSelectInput(session,
                      "gds_filelist",
                      choices = list.files(path = data_path, pattern = "\\.gds$"))
  })
  
  observeEvent(input$check_file, {
    if (is.null(input$gds_filelist) || input$gds_filelist == "") {
      showWarningToast("No GDS file selected")
      req(input$gds_filelist != "")
    }
    RV <<- NULL
    gds_filepath <<- file.path(data_path, input$gds_filelist)
    tryCatch({
      #Capture GDS summary details
      gds_summary <- capture.output({
        cat("Checking GDS file...\n")
        gdsinfo <- snpgdsSummary(gds_filepath)
        cat("Summary complete.\n")
      })
      
      output$file_summary <- renderPrint({
        gds_summary
      })
      sam <- length(gdsinfo$sample.id)
      gds.samples<<-gdsinfo$sample.id
      snp_len <- length(gdsinfo$snp.id)
      updateNumericInput(session, "sample_size", value = sam)
      output$GDS_sample <- renderText({
        sam
      })
      
      output$GDS_snps <- renderText({
        snp_len
      })
      
    }, error = function(e) {
      # return a safeError if a parsing error occurs
      output$file_summary <- renderPrint({
        return (paste0(
          "Conversion of VCF file to GDS file failed with the error:",
          e
        ))
      })
    })
  })
  
  ##Code to render population data upon selection of GDS data
  output$file_summary <- renderTable({
    req(input$csv_file)
    tryCatch({
      b <- colnames(pop_data)
      c <- colnames(tab2)
      updateSelectInput(session, "data_primarykey", choices = c)
      updateSelectInput(session, "sample_id", choices = b)
      updateSelectInput(session, "population_key", choices = b)
      cat("Checking population data\n")
    }, error = function(e) {
      stop(paste0("Error when rendering population details:", safeError(e)))
    })
  })
  
  
  
  #########################Genotype Matrix############################
  output$autosome_out <- renderText({
    paste("Autosome.only :", input$autosome)
  })
  
  eventsFilt <- reactive({
    list(input$update_popcolum,input$selectpopfiltOption)
  })
  
  #Add options to filter GDS matrix by samples
  observeEvent(input$filtsamples,{
    if(input$filtsamples)
    {
      output$filtsampselected <- reactive({
        return(TRUE)
      })
      outputOptions(output, "filtsampselected", suspendWhenHidden = FALSE)
      output$txtsamplefilt<-renderUI({
        textAreaInput("samplestofiltlist",label = "Enter sample IDs to retain:")
      })
      if(!is.null(pop_data))
      {
        
        output$radioSampleIdfilt<-renderUI({
          radioButtons("selectpopfiltOption","Or Select from population table", 
                      choices = c("Filter by group"="group",
                                  "Select from table"="selfromtable"),
                      selected ="group" 
                      )
        })
        
        output$selectsampleIdcol <- renderUI({
          selectInput("selectsampleidcolumn",
                      label = "Select sample ID column:",
                      choices = colnames(pop_data))
        })
      }
    }
    else{
      output$filtsampselected <- reactive({
        return(F)
      })
      outputOptions(output, "filtsampselected", suspendWhenHidden = FALSE)
      output$txtsamplefilt<-renderUI({
        return()
      })
      output$radioSampleIdfilt<-renderUI({
        return()
      })
    }
  })
  
  observeEvent(input$samplestofiltlist,{
    updateRadioButtons(session,"selectpopfiltOption",selected = character(0))
  })
  
  observeEvent(eventsFilt(),{
    req(input$gds_filelist)
    req(input$selectpopfiltOption)
    output$gds_summary <- renderPrint("Convert GDS to Genotype Matrix")
    if(input$selectpopfiltOption=="group")
    {
      options.choose<-colnames(pop_data)
      
      output$selectsamplegroups<-renderUI({
        selectInput("selSampGrpCol",
                    label = "Select group column",
                    choices=options.choose)
      })
      
      output$final_df<-renderDT({
        return()
      })
    }
    else
    {
      output$selectsamplegroups<-renderUI({
        return()
      })
      output$selectgrouptoretain <- renderUI({
        return()
      })
      output$final_df<-renderDT({
        getDataTable(pop_data,editable = F,caption = "Population data table")
      })
    }
  })
  
  observeEvent(input$selSampGrpCol, {
    if (input$selectpopfiltOption=="group")
    {
      output$selectgrouptoretain <- renderUI({
        grps <- unique(pop_data[, input$selSampGrpCol])
        selectInput(
          "retaingroup",
          label = "Select group samples to retain",
          choices = grps,
          multiple = T
        )
      })
    }
  })
  
  output$check_geno <- reactive({
    check <- is.null(snps) & input$gds_filelist != ''
    cat("Check SNP DF status;", check, "\n")
    return(check)
  })
  observeEvent(input$geno_filter, {
    output$final_df <- renderDT({
      return()
    })
  })
  
  output$gds_summary <- renderPrint("Convert GDS to Genotype Matrix")
  
  #########################Generate Genotype Matrix#############################
  observeEvent(input$get_geno, {
    tryCatch({
      withProgress(message = 'Converting GDS to Genotype Matrix', value = 0.1, {
        if (is.null(input$gds_filelist) || input$gds_filelist == "")
        {
          showWarningToast("No GDS file selected")
          req(input$gds_filelist)
        }
        
        incProgress(0.1, detail = "Starting")
        
        nosnps <- "No SNPs"
        if (file.exists(gds_filepath)) {
          f <- snpgdsOpen(gds_filepath)
          cat(gds_filepath, "\n")
          samplesids <- read.gdsn(index.gdsn(f, "sample.id"))
          sample_size <- input$sample_size
          if (!is.null(deleted_samples))
          {
            print(paste0(
              "Excluding deleted samples from sample list:",
              sum(samplesids %in% tab2$sample.id)
            ))
            samplesids <- samplesids[samplesids %in% tab2$sample.id]
          }
          random_ids <<- samplesids[1:sample_size]
          random_ids <<- random_ids[!is.na(random_ids)]
          deleted_samples <<- NULL
          incProgress(0.1, detail = "Reading File")
          
          if (input$geno_filter == 'sample') {
            cat("Checking input.......\n")
            fields <- c("sample_size",
                        "ld_cutoff",
                        "maf_rate",
                        "missing_rate",
                        "autosome")
            
            if (input$filtsamples)
            {
              if (input$samplestofiltlist != "")
              {
                samples.filt <- trimws(strsplit(input$samplestofiltlist, "\\n")[[1]])
                random_ids <<- samples.filt
              }
              else if (input$selectpopfiltOption == "group")
              {
                cat("Filtering based on samples:\n")
                random_ids <<- as.character(pop_data[pop_data[, input$selSampGrpCol] %in% as.character(input$retaingroup), input$selectsampleidcolumn])
              }
              else if (input$selectpopfiltOption == "selfromtable")
              {
                print(input$final_df_rows_selected)
                random_ids <<- as.character(pop_data[input$final_df_rows_selected, input$selectsampleidcolumn])
              }
            }
            
            if(length(random_ids)==0)
            {
              show_toast("No IDs selected",text = "No selected sample IDs")
              req(length(random_ids))
            }
            
            gds_summary <- capture.output({
              cat("Input size:", input$sample_size, "\n")
              cat("Checking GDS file...\n")
              cat("Samples included in filter:",
                  random_ids,
                  "\nFilters used:\n")
              for (n in fields)
                cat(n, "=", input[[n]], "\n")
              incProgress(0.1, detail = "Getting SNPSet")
              snpset <<- snpgdsLDpruning(
                f,
                sample.id = random_ids,
                ld.threshold = input$ld_cutoff,
                maf = as.numeric(input$maf_rate),
                missing.rate = as.numeric(input$missing_rate),
                autosome.only = input$autosome
              )
              snpset.id <<- unlist(unname(snpset))
              incProgress(0.1, detail = "Converting SNPSet ID")
              s <- snpgdsGetGeno(
                f,
                sample.id = random_ids,
                snp.id = snpset.id,
                with.id = TRUE
              )
              incProgress(0.1, detail = "Converting Geno")
              
              # # Add row and col names
              snps <<- data.frame(s$genotype)
              rownames(snps) <<- s$sample.id
              colnames(snps) <<- s$snp.id
              #Render the table
              cat("Summary complete.\n")
              head(snps)
              incProgress(0.1, detail = "Setting DataFrame")
            })
          }
          else if (input$geno_filter == 'complete')
          {
            gds_summary <- capture.output({
              cat("Checking GDS file...\n")
              cat("Samples:", samplesids, "\n")
              incProgress(0.1, detail = "Getting SNPSet")
              random_ids <<- samplesids
              snpset <<- read.gdsn(index.gdsn(f, "snp.id"))
              cat("SNPset:", snpset, "\n")
              incProgress(0.1, detail = "Converting SNPSet ID")
              snpset.id <<- unlist(unname(snpset))
              
              s <- snpgdsGetGeno(
                f,
                sample.id = random_ids ,
                snp.id = snpset.id,
                with.id = TRUE
              )
              incProgress(0.1, detail = "Converting Geno")
              
              # Add row and col names
              snps <<- data.frame(s$genotype)
              rownames(snps) <<- s$sample.id
              colnames(snps) <<- s$snp.id
              cat("Summary complete.\n")
              incProgress(0.1, detail = "Setting DataFrame")
            })
          }
          else
          {
            if (is.null(snpset.id) || is.null(random_ids))
            {
              if (input$pca_datafiles == TRUE)
                msg <- "PCA File or SNPSetID File not uploaded successfully"
              else
                msg <- "No sample ID or SNPSet ID found"
              showWarningToast(msg)
              req(random_ids, snpset.id)
            }
            
            incProgress(0.1, detail = "Sample ID & SNPSet ID ok")
            
            gds_summary <- capture.output({
              cat("Checking GDS file...\n")
              incProgress(0.2, detail = "Regenerate Genotype Matrix")
              s <- snpgdsGetGeno(
                f,
                sample.id = random_ids,
                snp.id = snpset.id,
                with.id = TRUE
              )
              
              # Add row and col names
              incProgress(0.2, detail = "Genotype Matrix Complete")
              snps <<- data.frame(s$genotype)
              rownames(snps) <<- s$sample.id
              colnames(snps) <<- s$snp.id
              cat("Summary complete.\n")
              incProgress(0.2, detail = "Saving DataFrame")
            })
          }
          snpgdsClose(f)
          incProgress(0.3, detail = "Saving DataFrame")
          output$gds_summary <- renderPrint({
            if (is.null(gds_summary)) {
              return()
            }
            gds_summary
          })
        }
        else{
          output$gds_summary <- renderPrint({
            return ("Selected file not found.")
          })
        }
        
      })
    }, error = function(e) {
      output$gds_summary <- renderPrint({
        snpgdsClose(f)
        return (paste0("Error generating genotype matrix:", safeError(e)))
      })
    })
  })
  
  
  ################################# Run PCA Process #################################
  
  output$pca_summary <- renderPrint(return())
  output$pca_dt <- renderDataTable(return())
  output$pca_status <- reactive({
    !is.null(pca)
  })
  outputOptions(output, "pca_status", suspendWhenHidden = FALSE)
  
  
  observeEvent(input$get_PCA, {
    if (is.null(input$gds_filelist) || input$gds_filelist == "") {
      showWarningToast("No GDS file selected")
    }
    req(input$gds_filelist, input$PCA_Thread)
    withProgress(message = 'Principal Component Analysis Started', value = 0.1, {
      tryCatch({
        gds_filepath <- file.path(data_path, input$gds_filelist)
        incProgress(0.1, detail = "Detect File Path")
        if (file.exists(gds_filepath)) {
          f <- snpgdsOpen(gds_filepath)
          incProgress(0.1, detail = "Starting PCA function")
          pca_summary <- capture.output({
            cat("Checking GDS file...\n")
            pca <<-
              snpgdsPCA(
                f,
                sample.id = random_ids,
                snp.id = snpset.id,
                autosome.only = as.logical(input$autosome),
                num.thread = input$PCA_Thread
              )
            cat("Summary complete.\n")
          })
          incProgress(0.5, detail = "PCA Complete")
          
          tab <- data.frame(
            sample.id = pca$sample.id,
            EV1 = pca$eigenvect[, 1],
            EV2 = pca$eigenvect[, 2],
            EV3 = pca$eigenvect[, 3],
            EV4 = pca$eigenvect[, 4],
            stringsAsFactors = FALSE
          )
          incProgress(0.2, detail = "Converting to Data Frame Tab2")
          snpgdsClose(f)
          tab2 <<- tab
          incProgress(0.1, detail = "Complete")
          output$pca_summary <- renderPrint({
            pca_summary
          })
          output$pca_dt <- renderDataTable(datatable(tab2, editable = F))
          c <- colnames(tab2)
          updateSelectInput(session, "data_primarykey", choices = c)
        }
        else {
          output$pca_summary <- renderPrint({
            showWarningToast("No GDS file found")
            return("GDS File not found")
          })
        }
        output$pca_status <- reactive({
          !is.null(pca)
        })
         outputOptions(output, "pca_status", suspendWhenHidden = FALSE)
      }, error = function(e) {
        snpgdsClose(f)
        output$pca_summary <- renderPrint({
          return (e)
        })
      })
    })
  })
  
  observeEvent(input$add_toPCA, {
    tryCatch({
      colnames(pop_data)[colnames(pop_data) == input$sample_id] <-
      input$data_primarykey
      
      tab2 <<- inner_join(tab2, pop_data[,c(input$data_primarykey,as.character(input$population_key))], by =input$data_primarykey)
      output$pca_dt <- renderDataTable(datatable(tab2, editable = F))
    }, error = function(e) {
      stop(safeError(e))
    })
    output$final_df <- renderTable({
      if (input$disp == "head") {
        return(head(tab2))
      }
      else {
        return(tab2)
      }
    })
  })
  
  ################################# Run Core Hunter #################################
  
  output$hunter_summary <- renderText("")
  
  observeEvent(input$get_core, {
    library(corehunter)
    core <- list()
    always_vector <- NULL
    never_vector <- NULL
    
    withProgress(message = 'Calculating Core Hunter Cores', value = 0.1, {
      tryCatch({
        if (!is.null(input$id_check_2))
        {
          core <-
            c(as.numeric(input$id_check_2),
              as.numeric(strsplit(input$core_values, ",")[[1]]))
        }
        else{
          if (is.null(input$core_values) || input$core_values == "") {
            showWarningToast("No core values selected")
          }
          
          req(input$core_values)
          core <- as.numeric(strsplit(input$core_values, ",")[[1]])
        }
        incProgress(0.1, detail = "Core selected")
        
        if ("Always Selected" %in% input$select_cores)
        {
          if (is.null(input$always_cores) || input$always_cores == "") {
            showWarningToast("Always Selected core values missing")
          }
          req(input$always_cores)
          always_vector <- unname(strsplit(input$always_cores, ' ')[[1]])
          incProgress(0.1, detail = "Always Selected Cores Confirmed")
        }
        
        if ("Never Selected" %in% input$select_cores)
        {
          if (is.null(input$never_cores) || input$never_cores == "") {
            showWarningToast("Never Selected core values missing")
          }
          req(input$never_cores)
          never_vector <- unname(strsplit(input$never_cores, ' ')[[1]])
          incProgress(0.1, detail = "Never Selected Cores Confirmed")
        }
        if (input$core_adv == FALSE)
        {
          objectives <- objective()
          incProgress(0.1, detail = "Objectives added")
        }
        else{
          objectives <-
            objective(input$obj_type,
                      input$obj_measure,
                      as.numeric(input$obj_weight))
          incProgress(0.1, detail = "Objectives added")
        }
        
        hunter_summary <- capture.output({
          cat("Calculating Core hunter cores...\n")
          incProgress(0.1, detail = "Reading Genotypes")
          geno <<- genotypes(snps, format = "biparental")
          incProgress(0.1, detail = "Initialize Core Hunter data")
          mydata <- coreHunterData(genotypes = geno)
          core <-
            sampleCoresAtSizes(mydata,
                               core,
                               always_vector,
                               never_vector,
                               objectives)
          incProgress(0.3, detail = "Cores Calculated")
          cat("Complete.\n")
        })
        output$hunter_summary <- renderPrint({
          hunter_summary
        })
      }, error = function(e) {
        output$hunter_summary <- renderPrint({
          return(e)
        })
      })
    })
  })
  
  #########################Add population data ##################
  output$csvtest <- reactive({
    !is.null(input$csv_file)
  })
  outputOptions(output, "csvtest", suspendWhenHidden = FALSE)
  output$pop_table <- renderTable(NULL)
  output$tab2_df <- renderTable(NULL)
  #This needs to be updated - doesn't ,m
  observeEvent(input$data_file, {
    req(input$data_file)
    tryCatch({
      deleted_samples <<- NULL
      tab2 <<- read.csv(input$data_file$datapath)
      random_ids <<- tab2[, 1]
      output$tab2_df <- renderTable({
        return(head(tab2))
      })
    }, error = function(e) {
      stop(safeError(e))
    })
  })
  
  observeEvent(input$data_file2, {
    tryCatch({
      snpset.id <<- readRDS(input$data_file2$datapath)
    }, error = function(e) {
      stop(safeError(e))
    })
  })
  
  observeEvent(input$show_df, {
    output$gds_summary <- renderPrint({
      capture.output({
        cat("Showing PCA Dataframe\n")
        incProgress(0.1, detail = "Reading dataframe")
      })
    })
    output$final_df <- renderDT({
      datatable(tab2, editable = F)
    })
  })
  #List of events to listen on to read population file
  toListen <- reactive({
    list(input$csv_file, input$header, input$sep, input$quote)
  })
  
  #Read population file
  observeEvent(toListen(), {
    req(input$csv_file)
    tryCatch({
      pop_data <<- read.csv(
        input$csv_file$datapath,
        header = input$header,
        sep = input$sep,
        quote = input$quote
      )
      b <- colnames(pop_data)
      updateSelectInput(session, "sample_id", choices = b)
      updateSelectInput(session, "population_key", choices = b)
      #Render the table
      output$pop_table <- renderDT({
        getDataTable(pop_data,editable=T)
      })
    }, error = function(e) {
      stop(safeError(e))
    })
  })
  
  output$contents <- renderTable(NULL)
  
  output$final_df <- renderDT(NULL)
  
  #Event to handle the addition of a new column to the population data table
  observeEvent(input$add_toDF, {
    # Reactive value to store the data
    pop_data_live <- reactiveValues(df = pop_data)
    # Update the data after editing the table
    if (is.null(input$column_name) || input$column_name == "") {
      show_toast(
        title = "No name provided!",
        text = "Enter the name of the column to add",
        type = "error",
        position = "center"
      )
      
      req(input$column_name)
    }
    pop_data_live$df[, input$column_name] <- ""  # Add a new column with no values
    #Render the table
    output$pop_table <- renderDT({
      getDataTable(pop_data_live$df,editable=T,caption = "Population data")
    })
    pop_data <<- pop_data_live$df
    b <- colnames(pop_data)
    updateSelectInput(session, "sample_id", choices = b)
    updateSelectInput(session, "population_key", choices = b)
  })
  
  #Option to edit individual cells
  observeEvent(input$pop_table_cell_edit, {
    info <- input$pop_table_cell_edit
    i <- info$row
    j <- info$col
    v <- info$value
    pop_data[i, j] <<- DT::coerceValue(v, pop_data[i, j])
  })
  
  #Check if rows have been selected
  observeEvent(input$pop_table_rows_selected, {
    #Check if rows have been selected
    if (length(input$pop_table_rows_selected) > 0)
    {
      updateRadioButtons(session, "updateopt", selected = "ids_selected")
      updateTextAreaInput(session, "pastedIds", value = "")
    }
    
    output$pop_colname <- reactive({
      !is.null(input$column_name) || input$column_name != ''
    })
    outputOptions(output, "pop_colname", suspendWhenHidden = FALSE)
  })
  
  #Check option used for updating table
  observeEvent(input$updateopt, {
    if (input$updateopt == "ids_pasted")
    {
      dataTableProxy("pop_table") %>% selectRows(NULL)
    }
    else{
      updateTextAreaInput(session, "pastedIds", value = "")
    }
  })
  
  # Render the "Clear Selection" button conditionally
  output$clearButton <- renderUI({
    if (length(input$pop_table_rows_selected) > 0) {
      actionButton("clear_selection", "Clear Selection")
    }
  })
  
  #Conditionally render text input and update button
  output$selsampleID <- renderUI({
    if (input$pastedIds != "")
    {
      selectInput("selsampleId",
                  label = "Select Sample ID column:",
                  choices = colnames(pop_data))
    }
  })
  
  output$popColnameText <- renderUI({
    if ((input$pastedIds != "" ||
         length(input$pop_table_rows_selected)) &
        input$column_name != "")
    {
      textInput("pop_groupname", "New value:", value = "")
    }
  })
  
  output$popColnameBtn <- renderUI({
    if ((input$pastedIds != "" ||
         length(input$pop_table_rows_selected)) &
        input$column_name != "")
    {
      actionButton("update_popcolum", "Update selected rows")
    }
  })
  
  
  observeEvent(input$clear_selection, {
    dataTableProxy("pop_table") %>% selectRows(NULL)
  })
  
  #Option to new add values to new column
  observeEvent(input$update_popcolum, {
    tryCatch({
      if (input$pop_groupname == "" || is.null(input$pop_groupname))
      {
        show_toast(
          title = "No value provided!",
          text = "Enter new value to be applied to selected IDs",
          type = "error",
          position = "center"
        )
        req(input$pop_groupname)
      }
      
      info <- input$pop_table_rows_selected
      
      if (input$pastedIds != "")
      {
        info <- trimws(strsplit(input$pastedIds, "\\n")[[1]])
        info <- which(pop_data[, input$selsampleId] %in% info)
      }
      pop_data[info, input$column_name] <<- input$pop_groupname
      #Render the table
      output$pop_table <- renderDT({
        getDataTable(pop_data,editable=T,caption = "Population data")
      })
      
    }, error = function(e) {
      stop(safeError(e))
    })
  })
  
  
  
  ######################### Download Button on PCA Graph ########################
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(sub("\\.gds$", "", input$gds_filelist), ".csv")
    },
    content = function(file) {
      write.csv(tab2, file, row.names = FALSE)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste0(sub("\\.gds$", "", input$gds_filelist), "_snpsetid.rds")
    },
    content = function(file) {
      saveRDS(snpset.id, file)
    }
  )
  
  #############################PCA graph ######################################
  
  output$scatterPlot <- renderScatterD3({
    input$data_file
    input$delete_rows
    input$refresh_selection
    req(tab2)
    
    col_var <- if (input$title == "")
      NULL
    else
      tab2[, input$title]
    symbol_var <-  NULL
    size_var <-  NULL
    auto_label <-  NULL
    zoom_on <- NULL
    
    scatterD3(
      data = NULL,
      x = tab2[, input$x],
      y = tab2[, input$y],
      xlab = input$x,
      ylab = input$y,
      hover_size = 4,
      hover_opacity = 1,
      col_var = col_var,
      col_lab = input$title,
      ellipses = input$ellipse,
      symbol_var = symbol_var,
      zoom_on = zoom_on,
      zoom_on_level = 3,
      labels_positions = auto_label,
      point_opacity = input$opacity,
      point_size = input$point_size,
      left_margin = 90,
      lasso = TRUE,
      tooltip_text = tab2[, "sample.id"],
      caption = list(
        title = "",
        subtitle = "",
        text = ""
      ),
      lasso_callback = "function(sel) { Shiny.setInputValue('lassoSelection', sel.data().map(function(d) {return d.tooltip_text})); }"
    )
  })
  
  output$table <- DT::renderDataTable(NULL)
  
  observe({
    # Watch for changes in lasso selection and perform actions in Shiny
    lasso_data <<- input$lassoSelection
    table_fun()
    
  })
  
  table_fun <- function() {
    if (!is.null(lasso_data)) {
      output$table <- DT::renderDataTable(getDataTable(tab2[tab2$sample.id %in% lasso_data, ],displength = 10))
      output$sel_samples <- renderPrint({
        cat(lasso_data)
      })
    }
    else{
      output$table <- DT::renderDataTable(getDataTable(tab2[tab2$sample.id %in% NULL, ], displength = 10))
    }
  }
  
  ######################### Getting points from Graph ########################
  
  observeEvent(
    c(
      input$refresh_selection,
      input$data_file,
      input$get_PCA,
      input$get_core,
      input$add_toDF
    ),
    {
      req(tab2)
      b <- colnames(tab2)
      columns_to_remove <- c("sample.id", "EV1", "EV2", "EV3", "EV4")
      b <- b[!b %in% columns_to_remove]
      updateSelectInput(session, "title", choices = b)
    }
  )
  
  observeEvent(input$delete_rows, {
    req(tab2, lasso_data)
    tryCatch({
      b <- tab2[tab2$sample.id %in% lasso_data, ]
      b <- b[input$table_rows_selected, 1]
      deleted_samples <<- append(deleted_samples, b)
      
      tab2 <<- tab2[!tab2$sample.id %in% b, ]
      random_ids <<- random_ids[!random_ids %in% b]
      lasso_data <<- lasso_data[!lasso_data %in% b]
    }, error = function(e) {
      stop(safeError(e))
    })
    
    output$brush_info <- renderPrint({
      cat(deleted_samples)
    })
    output$del_samples <- renderPrint({
      cat(deleted_samples)
    })
    table_fun()
  })
  
  ######################### Code for esquisse  ########################
  
  data_rv <- reactiveValues(data = tab2, name = "tab2")
  
  observe({
    data_rv$data <- tab2
    data_rv$name <- "tab2"
  })
  
  observeEvent(input$Refresh_p, {
    results <-
      esquisse_server(
        id = "esquisse",
        default_aes = c("fill", "color", "size", "shape", "group", "facet"),
        import_from = "env",
        data_rv = tab2
      )
    
    output$code <- renderPrint({
      results$code_plot
    })
    
    output$filters <- renderPrint({
      results$code_filters
    })
    
    output$data <- renderPrint({
      str(results$data)
    })
  })
  
  ################################ Code for histograms for MAF and Missing Rate ####################################
  
  observeEvent(input$csv_file,{
    output$pophistoColname <- renderUI({
      if (!is.null(pop_data))
      {
        selectInput("selgrouphisto",
                    label = "Select group by column:",
                    choices = colnames(pop_data))
      }
    })
    
    output$pophistosample <- renderUI({
      if (!is.null(pop_data))
      {
        selectInput("selsamplehisto",
                    label = "Select Sample ID column:",
                    choices = sort(colnames(pop_data)))
      }
    })
    
    output$pophistoupdateBtn<-renderUI({
      if (!is.null(pop_data))
      {
        actionButton("btnupdatehisto",label = "Update histogram")
      }
    })
    
  })
  
  #Check if column for grouping samples has been selected
  observeEvent(input$selgrouphisto,{
    req(pop_data)
    output$pophistogroups <- renderUI({
      if (!is.null(pop_data))
      {
        selectInput("selhistogroups",
                    label = "Select groups to display:",
                    choices = unique(pop_data[,input$selgrouphisto]),
                    multiple = T)
      }
    })
  })
  
  
  observe({
    tryCatch({
      selectedOption <- input$gds_filelist
      RV <<- NULL
      output$GDS_sample <- renderText({
        invisible()
      })
      output$GDS_snps <- renderText({
        invisible()
      })
    }, error = function(e) {
      stop(paste0("Error when selecting MAF and Missing rate: ", safeError(e)))
    })
  })
  
  getHistoData<-reactive(
  {
      input$btnupdatehisto
      gds_filepath <- file.path(data_path, input$gds_filelist)          
      output$filename_histo <- renderPrint({
        cat(input$gds_filelist)
      })
      
      tryCatch({
        genofile <- snpgdsOpen(gds_filepath, readonly = TRUE)  
        if (is.null(RV)) 
        {
          RV <<- snpgdsSNPRateFreq(genofile, with.snp.id = TRUE)
          histo.df <<- data.frame(RV)

        }
        else if(!is.null(histogroupcol)&!is.null(histosamplecol))
        {   histo.df<<-NULL
            RV<<-NULL
            snp.len<-NULL
            samps.len<-NULL
            groups.selected<-as.character(input$selhistogroups)
            if(length(groups.selected)==0)
            {
              groups.selected<-unique(as.character(pop_data[,histogroupcol]))
            }
            if(length(groups.selected)>4)
            {
              show_toast(
                title = "Too many groups!",
                text = "You have selected too many groups to visualise, try a maximum of 4",
                type = "error",
                position = "center"
              )
              histogroupcol<-histosamplecol<-NULL
              
              req(length(groups.selected)<=4)
            }
            for(g in groups.selected)
            {
               samples.g<-intersect(gds.samples,pop_data[pop_data[,histogroupcol] %in% g,histosamplecol])
               print(gds.samples)
               if(length(samples.g)==0)
               {
                  show_toast("No matching samples found",text = "No samples found!",position = "center")
                  req(length(samples.g)>0)
               }
                 
               RV.g<-snpgdsSNPRateFreq(genofile, with.snp.id = TRUE,sample.id = samples.g)
               
               df.g<-data.frame(RV.g)
               df.g$Group<-g
               histo.df<<-rbind(histo.df,df.g)
               snp.len<-paste0(snp.len,paste0(g,nrow(df.g),collapse = "="),collapse = ";")
               samps.len<-paste0(samps.len,paste0(g,length(samples.g),collapse = ";"),collapse = "=")
               RV<<-sapply(names(RV.g),function(x){
                 if(is.null(RV))
                   RV.g[[x]]
                 else
                   c(RV[[x]],RV.g[[x]])
               },USE.NAMES = T,simplify = F)
            }
          
      }
        snpgdsClose(genofile)  
        
      },error=function(e){
        snpgdsClose(genofile)
        stop(paste0("Error eading snpgds data file: ", safeError(e)))
      })
  })
  
  plotHisto<-function(histo.type)
  {
    tryCatch({
      xlab<-"Missing rate"
      gtitle<-"Missing Rate"
      if(histo.type!="MissingRate")
      {
        xlab="Minor allele frequency"
        gtitle<-"Minor Allele Frequency (MAF)"
      }
      
      if(!is.null(histogroupcol)&!is.null(histosamplecol))
      {
        
         ggplot() +
          geom_histogram(
            data = histo.df,
            aes(x = get(histo.type),fill=Group),
            binwidth = 0.01,
            color = "black",
            fill="black"
          ) +
          geom_histogram(
            data = histo.df[histo.df$MissingRate < input$n &
                              histo.df$MinorFreq > input$maf, ],
            aes(x = get(histo.type),fill=Group),
            binwidth = 0.01,
            color="black",
            fill="white"
          )+
          facet_wrap(~ Group, ncol = 2) +
          theme(strip.text = element_text(size = 20))+
          xlab(xlab) +
          ggtitle(gtitle) +
          theme_light()+
          theme(legend.position = "none")
      }
      else
      {
        ggplot() +
          geom_histogram(
            data = histo.df,
            aes(x = get(histo.type)),
            binwidth = 0.01,
            color = "black",
            fill = "black"
          ) +
          geom_histogram(
            data = histo.df[histo.df$MissingRate < input$n &
                              histo.df$MinorFreq > input$maf, ],
            aes(x = get(histo.type)),
            binwidth = 0.01,
            fill = "white"
          ) + xlab(xlab) +
          ggtitle(gtitle) +
          theme_light()
      }
    },error=function(e){
      stop(paste0("Error while generating histograms:",safeError(e)))
    })
  }
  
  output$histo <- renderPlot({
    tryCatch({
      req(input$gds_filelist)
      getHistoData()
      updateNumericInput(session, inputId = "missing_rate", value = input$n)
      updateNumericInput(session, inputId = "maf_rate", value = input$maf)

      plotHisto("MissingRate")
      
    }, error = function(e) {
      histogroupcol<-histosamplecol<-NULL
      stop(paste0("Error when generating histogram: ", safeError(e)))
    })
    
  })
  
  output$histo2 <- renderPlot({
    req(input$gds_filelist)
    getHistoData()
    plotHisto("MinorFreq")
  })
  
  #Check for inputs
  observeEvent(input$btnupdatehisto,{
    histogroupcol<<-input$selgrouphisto
    histosamplecol<<-input$selsamplehisto
  })
  
  ################################ Code to Get and set rJava Heap memory ####################################
  
  observeEvent(input$check_memory, {
    tryCatch({
      output$java_mem <- renderPrint({
        J("java.lang.Runtime")$getRuntime()$maxMemory() / (1024 ^ 4)
      })
    }, error = function(e) {
      output$java_mem <- renderPrint({
        return (e)
      })
    })
  })
  
  observeEvent(input$set_mem, {
    req(input$Corehunter_Me)
    mem_size <- paste("-Xmx", input$Corehunter_Me, "G", sep = "")
    
    tryCatch({
      output$java_mem <- renderPrint({
        options(java.parameters = mem_size)
        library(corehunter)
      })
    }, error = function(e) {
      output$java_mem <- renderPrint({
        return (e)
      })
    })
  })
  
  ########################################################################################################
}