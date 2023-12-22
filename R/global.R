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

# Set the maximum upload size to 10 GB
options(shiny.maxRequestSize = 10 * 1024 ^ 3)  # 10 GB

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