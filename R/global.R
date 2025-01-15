# Define global variables
library(shiny)
library(ggplot2)
library(SNPRelate)
library(shinyFiles)
library(dplyr)
library(DT)
library(esquisse)
library(scatterD3)
library(shinyWidgets)
library(shinythemes)
library(shinycssloaders)
library(rhandsontable)
library(DBI)
library(shinyjs)
library(bcrypt)
library(RSQLite)
library(rJava)
library(RColorBrewer)
# Database connection
# Set the maximum upload size to 10 GB
options(shiny.maxRequestSize = 10 * 1024^3)  # 10 GB
options(shiny.fullstacktrace = TRUE)
#' @title Generate plot
#' @description Function to generate PCA plot
#' @param x_var
#' @param y_var
#' @param color
#' @param tooltip
#' @param data_id
#' @examples
#' parse_btop("4A-40-AGC25TA5")
#' @keywords btop,blast
usernamelabel <<- "Guest"
generate_plot <- function(x_var,
                          y_var,
                          color_var = "black",
                          title,
                          width,
                          height) {
  if (title != "") {
    gg <- ggplot(
      tab2,
      aes_string(
        x = x_var,
        y = y_var,
        color = color_var,
        tooltip = "sample.id",
        data_id = "sample.id"
      )
    ) +
      labs(title = title) +
      geom_point(cex = 1) +
      labs(x = x_var, y = y_var) +
      scale_color_discrete(name = color_var) +
      theme_minimal() +
      theme(legend.position = "right") +
      geom_point_interactive(size = 1, hover_nearest = TRUE)
    
  } else {
    gg <- ggplot(tab2,
                 aes_string(
                   x = x_var,
                   y = y_var,
                   tooltip = "sample.id",
                   data_id = "sample.id"
                 )) +
      labs(title = "PCA Graph") +
      geom_point(cex = 1) +
      labs(x = x_var, y = y_var) +
      theme_minimal() +
      theme(legend.position = "right") +
      geom_point_interactive(size = 1, hover_nearest = TRUE)
    
  }
  return(gg)
}

sampleCoresAtSizes <- function(data,
                               sizes,
                               always_vector,
                               never_vector,
                               objectives) {
  for (size in sizes) {
    core <- if (!is.null(always_vector) && !is.null(never_vector)) {
      sampleCore(
        data,
        size = size,
        always.selected = always_vector,
        never.selected = never_vector,
        obj = objectives
      )
    } else if (!is.null(always_vector)) {
      sampleCore(data,
                 size = size,
                 always.selected = always_vector,
                 obj = objectives)
    } else if (!is.null(never_vector)) {
      sampleCore(data,
                 size = size,
                 never.selected = never_vector,
                 obj = objectives)
    } else {
      sampleCore(data = data,
                 obj = objectives,
                 size = size)
    }
    
    core100 <- as.data.frame(core$sel)
    tab2[[paste("Core", as.character(size), sep = "")]] <<-
      ifelse(tab2$sample.id %in% core100$`core$sel`, "c", "n")
    print(core)
    
  }
  
}

showWarningToast <- function(message) {
  show_toast(
    title = "Warning!",
    text = message,
    type = "warning",
    position = "center"
  )
}


getDataTable <- function(dataframe,
                         displength = 50,
                         editable = F,
                         caption = NULL)
{
  datatable(
    dataframe,
    options = list(
      editable = editable,
      caption = caption,
      pageLength = displength,
      lengthMenu = list(c(10, 25, 50, 100, -1), c('10', '25', '50', '100', 'All')),
      paging = T
    )
  )
}

checkUser<-function(user,session)
{
  if(is.null(user$username)||user$username=="")
  {
    show_toast(
      title = "Login!",
      text = "Login to use Lamington",
      type = "error",
      position = "center"
    )
    updateNavbarPage(session,"mainnavbar",selected = "User")
    return()
  }
}

table_fun <- function(tab2,lasso_data,output) {
  if (!is.null(lasso_data)) {
    output$table <- DT::renderDataTable(getDataTable(tab2[tab2$sample.id %in% lasso_data, ], displength = 10))
    output$sel_samples <- renderPrint({
      cat(lasso_data)
    })
  }
  else{
    output$table <- DT::renderDataTable(getDataTable(tab2[tab2$sample.id %in% NULL, ], displength = 10))
  }
}

getBlindColours<-function(num_colors)
{
  # Initialize an empty vector to store colors
  colors <- c('#DC3220','#006CD1','#FFC20A','#1AFF1A')  
  
  if(num_colors>4)
  {
    # Calculate the number of palettes needed
    num_palettes <- ceiling(num_colors / max(brewer.pal.info$maxcolors[brewer.pal.info$colorblind]))
    blindsets<-brewer.pal.info[brewer.pal.info$colorblind,]
    blindsets<-blindsets[blindsets$category == 'qual',]
    blindsets<-rownames(blindsets[order(blindsets$maxcolors,decreasing = T),])
    
    # Loop through the required number of palettes
    for (i in 1:length(blindsets)) {
      # Select a colorblind-friendly palette (e.g., "Dark2")
      palette_name <- blindsets[i] 
      # Get the maximum number of colors for the current palette
      max_colors <- brewer.pal.info[palette_name, "maxcolors"]
      # Generate colors from the current palette
      colors <- unique(c(colors, brewer.pal(max_colors, palette_name)))
      if(length(colors)>num_colors)
      {
        break
      }
    }
  }
  # Print the generated colors
  return(colors[1:num_colors])
}





