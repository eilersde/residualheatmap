library(shiny)
library(plotly)
library(shinythemes)
library(dplyr)

ui <- fluidPage(
    # Set theme
    theme = shinytheme("spacelab"),
    
    # Some help text
    h2("Interactive Residual Heat Maps"),
    h4("Heat map for model error analysis using Plotly"),
    
    # Vertical space
    tags$hr(),
    
    # Plotly Chart Area
    fluidRow(
        column(6, plotlyOutput(outputId = "base", height = "600px")),
        column(6, plotlyOutput(outputId = "select", height = "600px"))),
    
    tags$hr(),
    tags$blockquote("Select points to calculate the residual heat map on a smaller sample.")
    )