# Load global variables
source("global.R")

server <- function(input, output, session) {
  #Session based variables
  histo.df <- NULL
  RV <- NULL
  tab2 <- NULL
  pca <- NULL
  random_ids <- NULL
  snpset <- NULL
  geno <- NULL
  snpset.id <- NULL #Global variable that holds a set of SNP Ids
  snps <- NULL
  deleted_samples <- NULL
  lasso_data <- NULL
  pop_data <- NULL #Variable holding the population data
  pop_data_live <- NULL
  gds_filepath <- "" #Path to selected GDS file
  gds_data <- NULL #GDS data
  data_path <- "../GDS/" #Path where all GDS files are saved
  vcf_path <- "../VCFs/" #Path where VCF files can be saved and uploaded directly to Lamington
  db_path <- "../DB/" #Path where database is stored
  pca_data <- "../PCA/" #Path where PCA results are stored
  histogroupcol <- NULL #variable to hold column in population data table
  histosamplecol <- NULL
  gds_files <- NULL
  pca_files<- NULL
  pca_start <- NULL

  # Reactive value to track user session
  user_session <- reactiveValues(user = NULL)
  
  #Options to create VCFs, GDS and database directory
  tryCatch({
    if (!dir.exists(db_path))
    {
      dir.create(db_path, recursive = T)
    }
    lamington_db <- paste0(db_path, "/lamington_db.db")
    
    #Check if paths exists
    if (!dir.exists(data_path))
    {
      dir.create(data_path, recursive = T)
    }
    
    if (!dir.exists(vcf_path))
    {
      dir.create(vcf_path, recursive = T)
    }
    if(!dir.create(pca_data))
    {
      dir.create(pca_data,recursive = T)
    }
  }, error = function(e)
  {
    show_toast(
      "Error!",
      paste0("Error while creating directories", safeError(e)),
      type = "error",
      position = "center"
    )
  })
  
  ######################### Convert vcf to GDS format ############################
  output$output_text <- renderPrint(invisible())
  
  eventstochecklogin <- reactive({
    list(
      input$refresh_gds_list,
      input$vcf_location,
      input$vcf_desktop,
      input$vcf_server,
      input$GDS_filename
    )
  })
  
  observeEvent(eventstochecklogin, {
    tryCatch({
      checkUser(user_session$user, session)
    }, error = function(e)
    {
      show_toast("Error!",
                 safeError(e),
                 position = "center",
                 type = "error")
    })
  })
  
  observeEvent(input$convert_button, {
    tryCatch({
      GDS_filename <- paste(input$GDS_filename, ".gds", sep = "")
      output_gds_path <- file.path(data_path, GDS_filename)
      checkUser(user_session$user, session)
      withProgress(message = 'Converting VCF file to GDS format', value = 0.1, {
        if (input$vcf_location == "Client") {
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
        
        incProgress(0.1, detail = "Processing..This may take a while...")
        con <- dbConnect(SQLite(), lamington_db)
        file_exists <- dbGetQuery(con,
                                  paste0("SELECT 1 FROM data WHERE dataid = '", GDS_filename, "'"))
        if (nrow(file_exists) > 0) {
          cat("GDS exists!\n")
          show_toast(
            title = "GDS file name already exists",
            text = paste0(
              "A GDS file with that name '",
              GDS_filename,
              "' already exists"
            ),
            type = "error",
            position = "center"
          )
          dbDisconnect(con)
          return()
        }
        
        output_text <- capture.output({
          cat("Converting VCF to GDS...\n")
          snpgdsVCF2GDS_R(file_select, output_gds_path, method = "biallelic.only")
          # Insert new user
          dbExecute(
            con,
            "INSERT INTO data (dataid,username,path) VALUES (?,?,?)",
            list(
              GDS_filename,
              user_session$user$username,
              output_gds_path
            )
          )
          dbDisconnect(con)
          cat("Conversion complete.\n")
          incProgress(0.3, detail = "Completing..")
        })
        output$output_text <- renderPrint({
          output_text
        })
      })
    }
    , error = function(e) {
      if (!is.null(con))
        dbDisconnect(con)
      output$output_text <- renderPrint({
        print(paste0("Error while uploading VCF file:", safeError(e)))
      })
      
    })
  })
  
  observeEvent(input$refresh_list_server, {
    tryCatch({
      updateSelectInput(session,
                        "file2",
                        choices = list.files(path = data_path, pattern = "\\.(vcf|vcf.gz)$"))
    }, error = function(e)
    {
      show_toast(
        "Error!",
        paste0("Error while refreshing VCF files: ", safeError(e)),
        type = "error",
        position = "center"
      )
    })
    
  })
  
  
  ######################### Refresh file and Check GDS file ############################
  
  output$file_summary <- renderPrint(invisible())
  
  observeEvent(c(input$refresh_gds_list, input$convert_button), {
    RV <<- NULL
    tryCatch({
      con <- dbConnect(SQLite(), lamington_db)
      gds_files <<- dbGetQuery(
        con,
        paste0(
          "SELECT * FROM data WHERE username='",
          user_session$user$username,
          "'"
        )
      )
      gds_files <<- gds_files[order(gds_files$dataid, decreasing = F), ]
      updateSelectInput(session, "gds_filelist", choices = gds_files$dataid)
      dbDisconnect(con)
    }, error = function(e)
    {
      if (!is.null(con))
      {
        dbDisconnect(con)
      }
      output$file_summary <- renderPrint({
        print("Error while loading GDS files:", safeError(e))
      })
      
    })
  })
  
  observeEvent(input$check_gdsfile, {
    tryCatch({
      if (is.null(input$gds_filelist) || input$gds_filelist == "") {
        showWarningToast("No GDS file selected")
        req(input$gds_filelist != "")
      }
      RV <<- NULL
      gds_filepath <<- file.path(data_path, input$gds_filelist)
      
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
      gds.samples <<- gdsinfo$sample.id
      snp_len <- length(gdsinfo$snp.id)
      updateNumericInput(session, "sample_size", value = sam)
      output$GDS_sample <- renderText({
        sam
      })
      
      output$GDS_snps <- renderText({
        snp_len
      })
      
      #Option to check and load pca_data
      pca_start<<-gsub("[^a-zA-Z0-9 ]", "_",input$gds_filelist)
      print(pca_start)
      
      
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
  
  observeEvent(c(input$refresh_gds_list, input$check_gdsfile), {
    tryCatch({
      output$deleteGDSButton <- renderUI({
        tagList(tags$hr(), actionButton(
          "deleteGDS",
          HTML(
            "<span class='glyphicon glyphicon-remove'></span> Remove File"
          )
        ))
      })
    }, error = function(e)
    {
      show_toast(
        "Error!",
        paste0(
          "Error while adding button to delete GDS file: ",
          safeError(e)
        ),
        type = "error",
        position = "center"
      )
    })
  })
  rv_delete <- reactiveValues(confirm_delete = NULL)
  
  observeEvent(input$deleteGDS, {
    tryCatch({
      # Check if a file has been uploaded
      if (!is.null(input$gds_filelist) && input$gds_filelist != "") {
        # Show a confirmation dialog
        confirmSweetAlert(
          session = session,
          inputId = "confirm_delete_gds",
          type = "warning",
          title = paste0(
            "Are you sure, you want to delete '",
            input$gds_filelist,
            "'"
          ),
          text = "This action cannot be undone!",
          btn_labels = c("Cancel", "Delete")
        )
        
      } else {
        # If no file is uploaded, display a message
        output$file_summary <- renderText("No file to delete.")
      }
    }, error = function(e)
    {
      show_toast("Error!",
                 safeError(e),
                 type = "error",
                 position = "center")
    })
  })
  
  observeEvent(input$confirm_delete_gds, {
    rv_delete$confirm_delete <- input$confirm_delete_gds # Store the confirmation result
    tryCatch({
      if (rv_delete$confirm_delete) {
        gds_filepath <<- file.path(data_path, input$gds_filelist)
        # User confirmed deletion, proceed with deleting the file
        if (!file.exists(gds_filepath))
        {
          show_toast(
            "File not found",
            paste0(
              "This file '",
              input$gds_filelist,
              "' was not found on server, its records will nonetheless be removed from the database"
            ),
            position = "center"
          )
        }
        else
        {
          file.remove(gds_filepath)
        }
        gds_filepath <<- NULL
        #delete from database
        con <- dbConnect(SQLite(), lamington_db)
        gds_delete <- dbExecute(con,
                                paste0(
                                  "DELETE FROM data WHERE dataid='",
                                  input$gds_filelist,
                                  "'"
                                ))
        dbDisconnect(con)
        
        if (gds_delete > 0)
        {
          output$file_summary <- renderText(paste0(input$gds_filelist, " file deleted successfully!"))
          gds_files <<- gds_files[gds_files$dataid != input$gds_filelist, ]
          print(gds_files)
          updateSelectInput(session, "gds_filelist", choices = gds_files$dataid)
          output$deleteGDSButton <- renderUI({
            NULL
          })
        }
        else
          output$file_summary <- renderText(paste0(gds_delete, " file deleted!"))
        # Optionally reset the file input
      } else {
        # User cancelled deletion
        output$file_summary <- renderText("Deletion cancelled.")
      }
      rv_delete$confirm_delete <- NULL  # Reset the confirmation after processing
    }, error = function(e)
    {
      show_toast(
        "Error!",
        paste0("Error while removing file: ", safeError(e)),
        type = "error",
        position = "center"
      )
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
      show_toast(
        "Error!",
        text = paste0("Error when rendering population details:", safeError(e)),
        position = "center"
      )
    })
  })
  
  
  
  #########################Genotype Matrix############################
  output$autosome_out <- renderText({
    paste("Autosome.only :", input$autosome)
  })
  
  eventsFilt <- reactive({
    list(input$update_popcolum, input$selectpopfiltOptionEntered)
  })
  
  #Add options to filter GDS matrix by samples
  observeEvent(input$filtsamples, {
    cat("Option for filtsamples:", input$filtsamples, "\n")
    tryCatch({
      if (input$filtsamples)
      {
        output$filtsampselected <- reactive({
          return(TRUE)
        })
        outputOptions(output, "filtsampselected", suspendWhenHidden = FALSE)
        output$txtsamplefilt <- renderUI({
          textAreaInput("samplestofiltlist", label = "Enter sample IDs to retain:")
        })
        
        if (!is.null(pop_data))
        {
          output$radioSampleIdfilt <- renderUI({
            radioButtons(
              inputId = "selectpopfiltOption",
              label = "Or Select from population table",
              choices = list(
                "Filter by group" = "group",
                "Select from table" = "selfromtable"
              ),
              selected = character(0)
            )
          })
          
          output$selectsampleIdcol <- renderUI({
            selectInput(
              "selectsampleidcolumn",
              label = "Select sample ID column:",
              choices = colnames(pop_data)
            )
          })
        }
      }
      else{
        output$filtsampselected <- reactive({
          return(F)
        })
        outputOptions(output, "filtsampselected", suspendWhenHidden = FALSE)
        output$txtsamplefilt <- renderUI({
          return()
        })
        output$radioSampleIdfilt <- renderUI({
          return()
        })
      }
    }, error = function(e)
    {
      show_toast("Error!", e, type = "error", position = "center")
    })
  })
  
  observeEvent(input$samplelistEntered, {
    tryCatch({
      cat("Sample list event activated\n")
      updateRadioButtons(session, "selectpopfiltOption", selected = character(0))
      output$selectsamplegroups <- renderUI({
        return(NULL)
      })
      output$selectgrouptoretain <- renderUI(
        {
          return(NULL)
        }
      )
    }, error = function(e)
    {
      show_toast("Error!", e, type = "error", position = "center")
    })
  })
  
  observeEvent(c(input$update_popcolum, input$selectpopfiltOptionEntered), {
    req(input$gds_filelist)
    if (is.null(input$gds_filelist) ||
        input$gds_filelist == "" ||
        length(input$selectpopfiltOption) == 0)
    {
      cat("GDS file not selected\n")
      return()
    }
    output$gds_summary <- renderPrint(cat(
      "Sample filtering option set to '",
      input$selectpopfiltOption,
      "'\n"
    ))
    tryCatch({
      if (input$selectpopfiltOption == "group")
      {
        options.choose <- colnames(pop_data)
        output$selectsamplegroups <- renderUI({
          selectInput("selSampGrpCol",
                      label = "Select group column",
                      choices = options.choose)
        })
        
        output$gds_summary <- renderPrint(cat(
          "Filter by group options '",
          paste(options.choose, collapse = ", "),
          "'\n"
        ))
        output$final_df <- renderDT({
          return(NULL)
        })
      }
      else
      {
        output$selectsamplegroups <- renderUI({
          return()
        })
        output$selectgrouptoretain <- renderUI({
          return()
        })
        output$final_df <- renderDT({
          getDataTable(pop_data, editable = F, caption = "Population data table")
        })
      }
    }, error = function(e)
    {
      show_toast("Error!",
                 safeError(e),
                 type = "error",
                 position = "center")
    })
  })
  
  observeEvent(input$selSampGrpCol, {
    tryCatch({
      if (length(input$selectpopfiltOption) == 0)
      {
        req(length(input$selectpopfiltOption) > 0)
      }
      if (input$selectpopfiltOption == "group")
      {
        grps <- unique(pop_data[, input$selSampGrpCol])
        output$selectgrouptoretain <- renderUI({
          selectizeInput(
            "retaingroup",
            label = "Select group samples to retain",
            choices = NULL,
            multiple = T
          )
        })
        updateSelectizeInput(session,
                             "retaingroup",
                             choices = grps,
                             server = T)
        cat("Ran after rendering retain group\n")
      }
    }, error = function(e)
    {
      show_toast("Error!",
                 safeError(e),
                 type = "error",
                 position = "center")
    })
  })
  
  output$check_geno <- reactive({
    tryCatch({
      check <- is.null(snps) & input$gds_filelist != ''
      cat("Check SNP DF status;", check, "\n")
      return(check)
    }, error = function(e)
    {
      show_toast(
        title = "Error!",
        text = e,
        type = "error",
        position = "center"
      )
    })
  })
  
  observeEvent(input$geno_filter, {
    output$final_df <- renderDT({
      cat("Select filtering options!\n")
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
                random_ids <<- as.character(pop_data[input$final_df_rows_selected, input$selectsampleidcolumn])
                output$final_df <- renderDT({return(NULL)})
              }
            }
            
            if (length(random_ids) == 0)
            {
              show_toast("No IDs selected", text = "No selected sample IDs")
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
    tryCatch({
      if (is.null(input$gds_filelist) || input$gds_filelist == "") {
        showWarningToast("No GDS file selected")
      }
      req(input$gds_filelist, input$PCA_Thread)
      withProgress(message = 'Principal Component Analysis Started', value = 0.1, {
        gds_filepath <- file.path(data_path, input$gds_filelist)
        incProgress(0.1, detail = "Detect File Path")
        if (file.exists(gds_filepath)) {
          f <- snpgdsOpen(gds_filepath)
          incProgress(0.1, detail = "Running PCA function")
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
      })
    }, error = function(e) {
      snpgdsClose(f)
      output$pca_summary <- renderPrint({
        return (safeError(e))
      })
    })
  })
  
  observeEvent(input$add_toPCA, {
    tryCatch({
      cat("Key selected:", input$population_key, "\n")
      if (is.null(input$population_key) ||
          length(input$population_key) == 0)
      {
        show_toast(
          "Population column missing",
          "Population key/column not selected",
          type = "warning",
          position = "center"
        )
        return()
      }
      req(input$population_key)
      
      colnames(pop_data)[colnames(pop_data) == input$sample_id] <-
        input$data_primarykey
      
      tab2 <<- inner_join(tab2, pop_data[, c(input$data_primarykey,
                                             as.character(input$population_key))], by = input$data_primarykey)
      output$pca_dt <- renderDataTable(datatable(tab2, editable = F))
      output$final_df <- renderTable({
        if (input$disp == "head") {
          return(head(tab2))
        }
        else {
          return(tab2)
        }
      })
    }, error = function(e) {
      show_toast(
        title = "Error",
        type = "error",
        text = safeError(e),
        position = "center"
      )
    })
  })
  
  ################################# Run Core Hunter #################################
  
  output$hunter_summary <- renderText("")
  
  observeEvent(input$get_core, {
    library(corehunter)
    tryCatch({
      core <- list()
      always_vector <- NULL
      never_vector <- NULL
      
      withProgress(message = 'Calculating Core Hunter Cores', value = 0.1, {
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
      })
    }, error = function(e) {
      output$hunter_summary <- renderPrint({
        return(e)
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
      show_toast(
        title = "Error",
        type = "error",
        text = safeError(e),
        position = "center"
      )
    })
  })
  
  observeEvent(input$data_file2, {
    tryCatch({
      snpset.id <<- readRDS(input$data_file2$datapath)
    }, error = function(e) {
      show_toast("Error",
                 text = safeError(e),
                 type = "error",
                 position = "center")
    })
  })
  
  observeEvent(input$show_df, {
    tryCatch({
      output$gds_summary <- renderPrint({
        capture.output({
          cat("Showing PCA Dataframe\n")
          incProgress(0.1, detail = "Reading dataframe")
        })
      })
      output$final_df <- renderDT({
        datatable(tab2, editable = F)
      })
    }, error = function(e)
    {
      show_toast("Error",
                 text = safeError(e),
                 type = "error",
                 position = "center")
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
        getDataTable(pop_data, editable = T)
      })
    }, error = function(e) {
      show_toast("Error!",
                 safeError(e),
                 type = "error",
                 position = "center")
    })
  })
  
  output$contents <- renderTable(NULL)
  
  output$final_df <- renderDT(NULL)
  
  #Event to handle the addition of a new column to the population data table
  observeEvent(input$add_toDF, {
    tryCatch({
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
        getDataTable(pop_data_live$df,
                     editable = T,
                     caption = "Population data")
      })
      pop_data <<- pop_data_live$df
      b <- colnames(pop_data)
      updateSelectInput(session, "sample_id", choices = b)
      updateSelectInput(session, "population_key", choices = b)
    }, error = function(e)
    {
      show_toast("Error!",
                 safeError(e),
                 type = "error",
                 position = "center")
    })
  })
  
  #Option to edit individual cells
  observeEvent(input$pop_table_cell_edit, {
    tryCatch({
      info <- input$pop_table_cell_edit
      i <- info$row
      j <- info$col
      v <- info$value
      pop_data[i, j] <<- DT::coerceValue(v, pop_data[i, j])
    }, error = function(e)
    {
      show_toast("Error!",
                 safeError(e),
                 type = "error",
                 position = "center")
    })
  })
  
  #Check if rows have been selected
  observeEvent(input$pop_table_rows_selected, {
    tryCatch({
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
    }, error = function(e)
    {
      show_toast("Error!",
                 safeError(e),
                 type = "error",
                 position = "center")
    })
  })
  
  #Check option used for updating table
  observeEvent(input$updateopt, {
    tryCatch({
      if (input$updateopt == "ids_pasted")
      {
        dataTableProxy("pop_table") %>% selectRows(NULL)
      }
      else{
        updateTextAreaInput(session, "pastedIds", value = "")
      }
    }, error = function(e)
    {
      show_toast("Error!",
                 safeError(e),
                 type = "error",
                 position = "center")
    })
  })
  
  # Render the "Clear Selection" button conditionally
  output$clearButton <- renderUI({
    if (length(input$pop_table_rows_selected) > 0) {
      actionButton(
        "clear_selection",
        HTML(
          "<span class='glyphicon glyphicon-repeat'></span> Clear Selection"
        )
      )
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
      actionButton(
        "update_popcolum",
        HTML(
          "<span class='glyphicon glyphicon-saved'></span> Update selected rows"
        )
      )
    }
  })
  
  
  observeEvent(input$clear_selection, {
    dataTableProxy("pop_table") %>% selectRows(NULL)
  })
  
  #Option to add values to new column
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
        getDataTable(pop_data, editable = T, caption = "Population data")
      })
      
    }, error = function(e) {
      show_toast("Error!",
                 safeError(e),
                 type = "error",
                 position = "center")
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
    input$reset_selected_rows
    req(tab2)
    tryCatch({
      colors <- NULL
      col_var <- NULL
      if (input$title != "")
      {
        col_var <- tab2[, input$title]
        colors <- getBlindColours(nlevels(factor(col_var)))
        #cat(col_var," and colours selected are ", colors)
      }
      
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
        colors = colors,
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
    }, error = function(e)
    {
      show_toast(
        "Error!",
        paste0("Rendering of PCA encoutered an error:", safeError(e)),
        type = "error",
        position = "center"
      )
    })
  })
  
  output$table <- DT::renderDataTable(NULL)
  
  observe({
    # Watch for changes in lasso selection and update table below graph
    lasso_data <<- input$lassoSelection
    table_fun(tab2, lasso_data, output)
  })
  
  
  ######################### Getting points from Graph ########################
  
  observeEvent(
    c(
      input$refresh_selection,
      input$data_file,
      input$get_PCA,
      input$get_core,
      input$add_toDF,
      input$reset_selected_rows
    ),
    {
      req(tab2)
      tryCatch({
        b <- colnames(tab2)
        columns_to_remove <- c("sample.id", "EV1", "EV2", "EV3", "EV4")
        b <- b[!b %in% columns_to_remove]
        updateSelectInput(session, "title", choices = b)
      }, error = function(e)
      {
        show_toast("Error!",
                   paste0(safeError(e)),
                   type = "error",
                   position = "center")
      })
    }
  )
  
  observeEvent(input$delete_rows, {
    if (is.null(lasso_data) || is.null(input$table_rows_selected))
    {
      show_toast(
        "No rows selected!",
        "No rows have been selected for deletion or highlighted on the graph",
        type = "warning",
        position = "center"
      )
      req(tab2, lasso_data)
    }
    tryCatch({
      b <- tab2[tab2$sample.id %in% lasso_data, ]
      b <- b[input$table_rows_selected, 1]
      deleted_samples <<- append(deleted_samples, b)
      tab2 <<- tab2[!tab2$sample.id %in% b, ]
      random_ids <<- random_ids[!random_ids %in% b]
      lasso_data <<- lasso_data[!lasso_data %in% b]
      
      output$brush_info <- renderPrint({
        cat(deleted_samples)
      })
      output$del_samples <- renderPrint({
        cat(deleted_samples)
      })
      table_fun(tab2, lasso_data, output)
    }, error = function(e) {
      show_toast("Error",
                 text = safeError(e),
                 position = "center",
                 type = "error")
    })
  })
  
  observeEvent(input$reset_selected_rows, {
    tryCatch({
      refresh_PCA()
    }, error = function(e)
    {
      show_toast("Error",
                 text = safeError(e),
                 position = "center",
                 type = "error")
    })
  })
  
  ######################### Code for esquisse  ########################
  data_rv <- reactiveValues(data = tab2, name = "tab2")
  
  observe({
    data_rv$data <- tab2
    data_rv$name <- "tab2"
  })
  
  refresh_PCA <- function()
  {
    lasso_data <<- NULL
    table_fun(tab2, lasso_data, output)
    deleted_samples <<- NULL
    
    output$brush_info <- renderPrint({
      cat(deleted_samples)
    })
    output$del_samples <- renderPrint({
      cat(deleted_samples)
    })
    output$sel_samples <- renderPrint({
      cat(deleted_samples)
    })
  }
  
  observeEvent(input$Refresh_p, {
    tryCatch({
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
    }, error = function(e)
    {
      show_toast("Error",
                 text = safeError(e),
                 position = "center",
                 type = "error")
    })
  })
  
  ################################ Code for histograms for MAF and Missing Rate ####################################
  
  observeEvent(input$csv_file, {
    tryCatch({
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
      
      output$pophistoupdateBtn <- renderUI({
        if (!is.null(pop_data))
        {
          actionButton(
            "btnupdatehisto",
            label =
              HTML(
                "<span class='glyphicon glyphicon-refresh'></span> Update histogram"
              )
          )
          
        }
      })
    }, error = function(e)
    {
      show_toast("Error",
                 text = safeError(e),
                 position = "center",
                 type = "error")
    })
  })
  
  #Check if column for grouping samples has been selected
  observeEvent(input$selgrouphisto, {
    req(pop_data)
    tryCatch({
      output$pophistogroups <- renderUI({
        selectizeInput(
          "selhistogroups",
          label = "Select groups to display:",
          choices = NULL,
          multiple = T,
          options = list(maxOptions = 5)
        )
        
      })
      updateSelectizeInput(session,
                           "selhistogroups",
                           choices = unique(pop_data[, input$selgrouphisto]),
                           server = TRUE)
    }, error = function(e)
    {
      show_toast("Error",
                 text = safeError(e),
                 position = "center",
                 type = "error")
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
      show_toast(
        "Error!",
        paste0("Error when selecting MAF and Missing rate: ", safeError(e)),
        type = "error",
        position = "center"
      )
    })
  })
  
  getHistoData <- reactive({
    tryCatch({
      input$btnupdatehisto
      gds_filepath <- file.path(data_path, input$gds_filelist)
      output$filename_histo <- renderPrint({
        cat(input$gds_filelist)
      })
      genofile <- snpgdsOpen(gds_filepath, readonly = TRUE)
      if (is.null(RV))
      {
        RV <<- snpgdsSNPRateFreq(genofile, with.snp.id = TRUE)
        histo.df <<- data.frame(RV)
        
      }
      else if (!is.null(histogroupcol) & !is.null(histosamplecol))
      {
        histo.df <<- NULL
        RV <<- NULL
        snp.len <- NULL
        samps.len <- NULL
        groups.selected <- as.character(input$selhistogroups)
        if (length(groups.selected) == 0)
        {
          groups.selected <- unique(as.character(pop_data[, histogroupcol]))
        }
        if (length(groups.selected) > 4)
        {
          show_toast(
            title = "Too many groups!",
            text = "Please note that there is a limit of 4 groups for display. Select a maximum of four groups to proceed.",
            type = "warning",
            position = "center"
          )
          histogroupcol <- histosamplecol <- NULL
          
          req(length(groups.selected) <= 4)
        }
        for (g in groups.selected)
        {
          samples.g <- intersect(gds.samples, pop_data[pop_data[, histogroupcol] %in% g, histosamplecol])
          if (length(samples.g) == 0)
          {
            show_toast(text = "No matching samples found",
                       title = "No samples found!",
                       position = "center")
            req(length(samples.g) > 0)
          }
          
          RV.g <- snpgdsSNPRateFreq(genofile,
                                    with.snp.id = TRUE,
                                    sample.id = samples.g)
          
          df.g <- data.frame(RV.g)
          df.g$Group <- g
          histo.df <<- rbind(histo.df, df.g)
          snp.len <- paste0(snp.len, paste0(g, nrow(df.g), collapse = "="), collapse = ";")
          samps.len <- paste0(samps.len,
                              paste0(g, length(samples.g), collapse = ";"),
                              collapse = "=")
          RV <<- sapply(names(RV.g), function(x) {
            if (is.null(RV))
              RV.g[[x]]
            else
              c(RV[[x]], RV.g[[x]])
          }, USE.NAMES = T, simplify = F)
        }
        
      }
      snpgdsClose(genofile)
      
    }, error = function(e) {
      snpgdsClose(genofile)
      show_toast(
        "Error!",
        paste0("Error reading snpgds data file: ", safeError(e)),
        type = "error",
        position = "center"
      )
    })
  })
  
  plotHisto <- function(histo.type)
  {
    tryCatch({
      xlab <- "Missing rate"
      gtitle <- "Missing Rate"
      if (histo.type != "MissingRate")
      {
        xlab = "Minor allele frequency"
        gtitle <- "Minor Allele Frequency (MAF)"
      }
      
      if (!is.null(histogroupcol) & !is.null(histosamplecol))
      {
        ggplot() +
          geom_histogram(
            data = histo.df,
            aes(x = get(histo.type), fill = Group),
            binwidth = 0.01,
            color = "black",
            fill = "black"
          ) +
          geom_histogram(
            data = histo.df[histo.df$MissingRate < input$n &
                              histo.df$MinorFreq > input$maf, ],
            aes(x = get(histo.type), fill = Group),
            binwidth = 0.01,
            color = "black",
            fill = "white"
          ) +
          facet_wrap( ~ Group, ncol = 2) +
          theme(strip.text = element_text(size = 20)) +
          xlab(xlab) +
          ggtitle(gtitle) +
          theme_light() +
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
    }, error = function(e) {
      show_toast(
        "Error!",
        paste0("Error while generating histograms:", safeError(e)),
        type = "error",
        position = "center"
      )
    })
  }
  
  output$histo <- renderPlot({
    tryCatch({
      if (is.null(input$gds_filelist) || input$gds_filelist == "")
      {
        show_toast(
          title = "GDS missing!",
          text = "Load GDS file to generate histogram",
          type = "error",
          position = "center"
        )
        return()
      }
      
      if (input$gds_filelist != "" && !is.null(input$gds_filelist))
      {
        getHistoData()
        updateNumericInput(session, inputId = "missing_rate", value = input$n)
        updateNumericInput(session, inputId = "maf_rate", value = input$maf)
        plotHisto("MissingRate")
      }
    }, error = function(e) {
      histogroupcol <- histosamplecol <- NULL
      show_toast(
        "Error!",
        paste0("Error when generating histogram: ", safeError(e)),
        type = "error",
        position = "center"
      )
    })
    
  })
  
  output$histo2 <- renderPlot({
    req(input$gds_filelist)
    getHistoData()
    plotHisto("MinorFreq")
  })
  
  #Check for inputs
  observeEvent(input$btnupdatehisto, {
    histogroupcol <<- input$selgrouphisto
    histosamplecol <<- input$selsamplehisto
  })
  
  ################################ Code to Get and set rJava Heap memory ####################################
  
  observeEvent(input$check_memory, {
    tryCatch({
      output$java_mem <- renderPrint({
        J("java.lang.Runtime")$getRuntime()$maxMemory() / (1024^4)
      })
    }, error = function(e) {
      output$java_mem <- renderPrint({
        return (paste0(
          "Error obtaining total Java runtime max memory:",
          safeError(e)
        ))
      })
    })
  })
  
  observeEvent(input$set_mem, {
    tryCatch({
      checkUser(user_session$user, session)
      cat("Setting heap")
      if (is.null(input$Corehunter_Me) ||
          input$Corehunter_Me == 0 || input$Corehunter_Me == "")
      {
        show_toast("",
                   text = "Heapmemory should be > 0",
                   type = "error",
                   position = "center")
        return()
      }
      mem_size <- paste0("-Xmx", input$Corehunter_Me, "G")
      
      output$java_mem <- renderPrint({
        print(options(java.parameters = mem_size))
        print(library(corehunter))
      })
    }, error = function(e) {
      output$java_mem <- renderPrint({
        return (paste0("Error while setting heap memory;", safeError(e)))
      })
    })
  })
  
  ################################Code to manage user and data profiles##########################################
  tryCatch({
    con <- dbConnect(SQLite(), lamington_db)
    # Check if users table exists
    if (!dbExistsTable(con, "users")) {
      # If the table doesn't exist, create it
      dbCreateTable(
        con,
        "users",
        fields = c(
          username = "TEXT PRIMARY KEY",
          password = "TEXT",
          name = "TEXT",
          email = "TEXT"
        )
      )
      cat("Table 'users' created.\n")
    } else {
      cat("Table 'users' already exists.\n")
    }
    
    # Check if data table exists
    if (!dbExistsTable(con, "data")) {
      # Create the data table with a foreign key referencing users(userid)
      dbExecute(
        con,
        "CREATE TABLE data (
        dataid TEXT PRIMARY KEY,
        username TEXT,
        path TEXT,
        FOREIGN KEY (username) REFERENCES users(username)
      )"
      )
      cat("Table 'data' created.\n")
    } else {
      cat("Table 'data' already exists.\n")
    }
    dbDisconnect(con)
  }, error = function(e)
  {
    if (!is.null(con))
      dbDisconnect(con)
    output$userlog <- renderPrint({
      return (paste0(
        "Error while preprocessing database connection :",
        safeError(e)
      ))
    })
  })
  output$loginoption <- reactive({
    return("login")
  })
  outputOptions(output, "loginoption", suspendWhenHidden = FALSE)
  
  
  ###Function used to login to lamington
  loginFunction <- function()
  {
    tryCatch({
      if (input$username == "" || input$password == "")
      {
        show_toast(
          title = "Username/password missing!",
          text = "Username and password can't be left blank",
          type = "error",
          position = "center"
        )
        req(input$username != "" && input$password != "")
      }
      
      con <- dbConnect(SQLite(), lamington_db)
      user <- dbGetQuery(con,
                         paste0("SELECT * FROM users WHERE username = '", input$username, "'"))
      
      if (nrow(user) == 1 &&
          bcrypt::checkpw(input$password, user$password)) {
        user_session$user <- user
        output$loginoption <- reactive({
          return("log")
        })
        outputOptions(output, "loginoption", suspendWhenHidden = FALSE)
        output$signupoptions <- renderUI({
          tagList(
            textInput("name", "Name", placeholder = user$name),
            textInput("email", "Email address", placeholder = user$email),
            tags$hr(),
            actionButton(
              "logoutBtn",
              HTML(
                "<span class='glyphicon glyphicon-log-out'></span> Logout"
              )
            ),
          )
        })
        show_toast(
          title = "Login success!",
          text = paste0("Hello ", user$username),
          type = "success",
          position = "center"
        )
        
        output$userlog <- renderPrint({
          cat("Using Lamington as ", user$username, "\n")
        })
        
        output$testname <- renderUI({
          tagList(
            tags$label(HTML(
              paste0(
                "<span class='glyphicon glyphicon-user'></span> ",
                user$username
              )
            )),
            tags$label("("),
            actionLink(
              "logout",
              HTML(
                "<span class='glyphicon glyphicon-log-out'></span> Logout"
              )
            ),
            tags$label(")")
          )
        })
        
      } else {
        showModal(modalDialog(
          title = tags$div(icon("exclamation-triangle"), "Login Failed"),
          tags$p(
            style = "color: red;",
            "Invalid username and/or password. Please try again."
          ),
          easyClose = TRUE,
          fade = T,
          footer = NULL
        ))
      }
      dbDisconnect(con)
    }, error = function(e)
    {
      if (!is.null(con))
        dbDisconnect(con)
      output$userlog <- renderPrint({
        return (paste0("Login failed with :", safeError(e)))
      })
    })
  }
  
  observeEvent(input$passwordEntered, {
    loginFunction()
  })
  
  # --- Login ---
  observeEvent(input$login, {
    loginFunction()
  })
  
  # --- Signup ---
  observeEvent(input$signup, {
    tryCatch({
      # Basic input validation (add more as needed)
      output$loginoption <- reactive({
        return("signup")
      })
      outputOptions(output, "loginoption", suspendWhenHidden = FALSE)
      output$signupoptions <- renderUI({
        tagList(
          textInput("name", "Name"),
          textInput("email", "Email address"),
          tags$hr(),
          actionButton(
            "createuser",
            HTML(
              "<span class='glyphicon glyphicon-record'></span> Sign up"
            )
          ),
          actionButton(
            "cancelcreateuser",
            HTML(
              "<span class='glyphicon glyphicon-remove-circle'></span> Cancel"
            )
          )
        )
      })
      output$userlog <- renderPrint({
        cat("Adding new user details!\n")
      })
    }, error = function(e)
    {
      show_toast("Error!",
                 safeError(e),
                 type = "error",
                 position = "center")
    })
  })
  
  observeEvent(input$createuser, {
    tryCatch({
      if (input$username == "" ||
          input$password == "" || input$email == "")
      {
        show_toast(
          title = "Username/password missing!",
          text = "Please fill in all required fields: Username, Password, and Email.",
          type = "error",
          position = "center"
        )
        req(input$username != "" && input$password != "")
      }
      
      # Check if username exists
      con <- dbConnect(SQLite(), lamington_db)
      user_exists <- dbGetQuery(
        con,
        paste0(
          "SELECT 1 FROM users WHERE (username = '",
          input$username,
          "' OR email = '",
          input$email,
          "') AND email IS NOT NULL"
        )
      )
      if (nrow(user_exists) > 0) {
        show_toast(
          title = "Username already exists",
          text = paste0(
            "The username:'",
            input$username,
            "' and/or email:'",
            input$email,
            "' already exists"
          ),
          type = "error",
          position = "center",
          timer = NULL
        )
        return()
      }
      
      # Hash the password
      hashed_password <- bcrypt::hashpw(input$password)
      # Insert new user
      dbExecute(
        con,
        "INSERT INTO users (username, password,name,email) VALUES (?,?,?,?)",
        list(
          input$username,
          hashed_password,
          input$name,
          input$email
        )
      )
      
      show_toast(
        title = "Account creation success!",
        text = "Account created succefully. Please login.",
        type = "success",
        position = "center"
      )
      
      output$userlog <- renderText({
        cat("New user details added\n")
        cat("Username:", input$username, "\n")
        cat("Name:", input$name, "\n")
        cat("Email:", input$email, "\n")
      })
      
      cat("Disconnect db connection\n")
      dbDisconnect(con)
      output$signupoptions <- renderUI({
        return()
      })
      output$loginoption <- reactive({
        return("login")
      })
      outputOptions(output, "loginoption", suspendWhenHidden = FALSE)
    }, error = function(e) {
      if (!is.null(con))
        dbDisconnect(con)
      output$userlog <- renderPrint({
        cat(paste0("Error adding new user details:", safeError(e)),
            "\n")
      })
    })
  })
  
  observeEvent(input$cancelcreateuser, {
    tryCatch({
      # Basic input validation (add more as needed)
      output$loginoption <- reactive({
        return("login")
      })
      outputOptions(output, "loginoption", suspendWhenHidden = FALSE)
      output$signupoptions <- renderUI({
        return()
      })
      output$userlog <- renderPrint({
        cat("Login to lamington!\n")
      })
    }, error = function(e) {
      show_toast(
        "Error!",
        text = safeError(e),
        type = "error",
        position = "center"
      )
    })
  })
  
  
  eventstochecklogout <- reactive({
    list(input$logout, input$logoutBtn)
  })
  #----logout------------------
  observeEvent(input$logout, {
    tryCatch({
      # Clear all session data
      session$reload()
      showModal(
        modalDialog(
          title = "Logged Out",
          "You have been logged out. Lamington will restart.",
          easyClose = TRUE
        )
      )
    }, error = function(e)
    {
      show_toast(
        "Error!",
        text = safeError(e),
        type = "error",
        position = "center"
      )
    })
  })
  
  observeEvent(input$logoutBtn, {
    # Clear all session data
    session$reload()
    showModal(
      modalDialog(
        title = "Logged Out",
        "You have been logged out. Lamington will restart.",
        easyClose = TRUE
      )
    )
  })
  
}