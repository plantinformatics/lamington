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

# Set the maximum upload size to 10 GB
options(shiny.maxRequestSize = 10 * 1024 ^ 3)  # 10 GB

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

generate_plot <- function(x_var, y_var, color_var ="black", title, width, height) {
  
  if (title != "") {
      gg <- ggplot(tab2, aes_string(x = x_var, y = y_var, color = color_var, tooltip = "sample.id", data_id = "sample.id")) +
      labs(title = title) +
      geom_point(cex = 1) +
      labs(x = x_var, y = y_var) +
      scale_color_discrete(name = color_var) +
      theme_minimal() +
      theme(legend.position = "right") +
       geom_point_interactive(
       size = 1, hover_nearest = TRUE)
    
  } else {
     gg <- ggplot(tab2, aes_string(x = x_var, y = y_var, tooltip = "sample.id", data_id = "sample.id")) +
      labs(title = "PCA Graph") +
      geom_point(cex = 1) +
      labs(x = x_var, y = y_var) +
      theme_minimal() +
      theme(legend.position = "right") +
        geom_point_interactive(
          size = 1, hover_nearest = TRUE)
   
  }
  return(gg)
}

sampleCoresAtSizes <- function(data, sizes,always_vector,never_vector,objectives) {
  
  for (size in sizes) {

    core <- if (!is.null(always_vector) && !is.null(never_vector)) {
    sampleCore(
      data,
      size = size,
      always.selected = always_vector,
      never.selected = never_vector,
      obj=objectives
    )
    } else if (!is.null(always_vector)) {
      sampleCore(data, size = size, always.selected = always_vector, obj=objectives)
    } else if (!is.null(never_vector)) {
      sampleCore(data, size = size, never.selected = never_vector, obj=objectives)
    } else {
      
      sampleCore(data = data,obj = objectives,size = size)
    }
    
    core100 <- as.data.frame(core$sel)
    tab2[[paste("Core", as.character(size), sep = "")]] <<-
      ifelse(tab2$sample.id %in% core100$`core$sel`, "c", "n")
    print(core)

  }

}

showWarningToast <- function(message) {
  show_toast(
    title = "Careful!",
    text = message,
    type = "warning",
    position = "center"
  )
}


getDataTable<-function(dataframe,displength=50,editable=F,caption = NULL)
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


