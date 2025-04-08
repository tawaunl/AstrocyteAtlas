#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")

sigsearchdef <- c("Myoc","Fxyd6","Gfap","6330403K07Rik","Igfbp5","Serpinf1",
                  "Vim","Fxyd1","S100a6")
# Define UI for application that draws a histogram
shinyUI(fluidPage(
  tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))), 
  list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))), 
  titlePanel("Mouse Interactive Shiny"),  
  navbarPage( 
    NULL, 
    # GeneScore panel ----------------
    tabPanel( 
      HTML("Gene Signature Search"), 
      h4("Gene Signature Search"), 
      "In this tab, users can visualise the average expression patterns of ", 
      "multiple genes grouped by categorical cell information or plotted on a dimensional reduction of choice", br(), 
      "The normalised expression are averaged, log-transformed and then plotted.", 
      br(),br(), 
      fluidRow( 
        column( 
          3, style="border-right: 2px solid black", 
          textAreaInput("siginp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                        height = "300px", 
                        value = paste0(sigsearchdef, collapse = ", ")) %>% 
            helper(type = "inline", size = "m", fade = TRUE, 
                   title = "List of genes to create signature", 
                   content = c("Input genes", 
                               "- Maximum 50 genes", 
                               "- Genes should be separated by comma, semicolon or newline")), 
          selectInput("siggrp", "Group by (Bar Plot only):", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1conf[grp == TRUE]$UI[11]) %>% 
            helper(type = "inline", size = "m", fade = TRUE, 
                   title = "Cell information to group cells by", 
                   content = c("Select categorical cell information to group cells by", 
                               "- Single cells are grouped by this categorical covariate", 
                               "- Plotted as the X-axis of the Bar Plot")), 
          radioButtons("sigplt", "Plot type:", 
                       choices = c("UMAP", "BarPlot"), 
                       selected = "UMAP", inline = TRUE),  
          br(), 
          actionButton("sigtogL", "Toggle to subset cells"), 
          conditionalPanel( 
            condition = "input.sigtogL % 2 == 1", 
            selectInput("sigsub1", "Cell information to subset:", 
                        choices = sc1conf[grp == TRUE]$UI, 
                        selected = sc1def$grp1), 
            uiOutput("sigsub1.ui"), 
            actionButton("sigsub1all", "Select all groups", class = "btn btn-primary"), 
            actionButton("sigsub1non", "Deselect all groups", class = "btn btn-primary") 
          ), br(), br(), 
          actionButton("sigtog", "Toggle graphics controls"), 
          conditionalPanel( 
            condition = "input.sigtog % 2 == 1", 
            radioButtons("sigcols", "Colour scheme:", 
                         choices = c("White-Red", "Blue-Yellow-Red", 
                                     "Yellow-Green-Purple"), 
                         selected = "Blue-Yellow-Red"), 
            radioButtons("sigord1", "Plot order:", 
                         choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                         selected = "Original", inline = TRUE),
            radioButtons("sigpsz", "Plot size:", 
                         choices = c("Small", "Medium", "Large"), 
                         selected = "Medium", inline = TRUE), 
            radioButtons("sigfsz", "Font size:", 
                         choices = c("Small", "Medium", "Large"), 
                         selected = "Medium", inline = TRUE)) 
        ), # End of column (6 space) 
        column(9, h4(htmlOutput("sigoupTxt")), 
               uiOutput("sigoup.ui"), 
               downloadButton("sigoup.pdf", "Download PDF"), 
               downloadButton("sigoup.png", "Download PNG"), br(), 
               div(style="display:inline-block", 
                   numericInput("sigoup.h", "PDF / PNG height:", width = "138px", 
                                min = 4, max = 20, value = 10, step = 0.5)), 
               div(style="display:inline-block", 
                   numericInput("sigoup.w", "PDF / PNG width:", width = "138px", 
                                min = 4, max = 20, value = 10, step = 0.5)) 
        )  # End of column (6 space) 
      )    # End of fluidRow (4 space) 
    ),    # End of tab (2 space) 
    # Bubble plot ----------------
    tabPanel( 
      HTML("Bubbleplot / Heatmap"), 
      h4("Gene expression bubbleplot / heatmap"), 
      "In this tab, users can visualise the gene expression patterns of ", 
      "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
      "The normalised expression are averaged, log-transformed and then plotted.", 
      br(),br(), 
      fluidRow( 
        column( 
          3, style="border-right: 2px solid black", 
          textAreaInput("sc1d1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                        height = "400px", 
                        value = paste0(sc1def$genes, collapse = ", ")) %>% 
            helper(type = "inline", size = "m", fade = TRUE, 
                   title = "List of genes to plot on bubbleplot / heatmap", 
                   content = c("Input genes to plot", 
                               "- Maximum 50 genes (due to ploting space limitations)", 
                               "- Genes should be separated by comma, semicolon or newline")), 
          selectInput("sc1d1grp", "Group by:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1conf[grp == TRUE]$UI[11]) %>% 
            helper(type = "inline", size = "m", fade = TRUE, 
                   title = "Cell information to group cells by", 
                   content = c("Select categorical cell information to group cells by", 
                               "- Single cells are grouped by this categorical covariate", 
                               "- Plotted as the X-axis of the bubbleplot / heatmap")), 
          radioButtons("sc1d1plt", "Plot type:", 
                       choices = c("Bubbleplot", "Heatmap"), 
                       selected = "Bubbleplot", inline = TRUE), 
          checkboxInput("sc1d1scl", "Scale gene expression", value = TRUE), 
          checkboxInput("sc1d1row", "Cluster rows (genes)", value = FALSE), 
          checkboxInput("sc1d1col", "Cluster columns (samples)", value = FALSE), 
          br(), 
          actionButton("sc1d1togL", "Toggle to subset cells"), 
          conditionalPanel( 
            condition = "input.sc1d1togL % 2 == 1", 
            selectInput("sc1d1sub1", "Cell information to subset:", 
                        choices = sc1conf[grp == TRUE]$UI, 
                        selected = sc1def$grp1), 
            uiOutput("sc1d1sub1.ui"), 
            actionButton("sc1d1sub1all", "Select all groups", class = "btn btn-primary"), 
            actionButton("sc1d1sub1non", "Deselect all groups", class = "btn btn-primary") 
          ), br(), br(), 
          actionButton("sc1d1tog", "Toggle graphics controls"), 
          conditionalPanel( 
            condition = "input.sc1d1tog % 2 == 1", 
            radioButtons("sc1d1cols", "Colour scheme:", 
                         choices = c("White-Red", "Blue-Yellow-Red", 
                                     "Yellow-Green-Purple"), 
                         selected = "Blue-Yellow-Red"), 
            radioButtons("sc1d1psz", "Plot size:", 
                         choices = c("Small", "Medium", "Large"), 
                         selected = "Medium", inline = TRUE), 
            radioButtons("sc1d1fsz", "Font size:", 
                         choices = c("Small", "Medium", "Large"), 
                         selected = "Medium", inline = TRUE)) 
        ), # End of column (6 space) 
        column(9, h4(htmlOutput("sc1d1oupTxt")), 
               uiOutput("sc1d1oup.ui"), 
               downloadButton("sc1d1oup.pdf", "Download PDF"), 
               downloadButton("sc1d1oup.png", "Download PNG"), br(), 
               div(style="display:inline-block", 
                   numericInput("sc1d1oup.h", "PDF / PNG height:", width = "138px", 
                                min = 4, max = 20, value = 10, step = 0.5)), 
               div(style="display:inline-block", 
                   numericInput("sc1d1oup.w", "PDF / PNG width:", width = "138px", 
                                min = 4, max = 20, value = 10, step = 0.5)) 
        )  # End of column (6 space) 
      )    # End of fluidRow (4 space) 
    ),    # End of tab (2 space) 
    
    # End of Tabs ------------------
    br(), 
    p("", style = "font-size: 125%;"), 
    p(em("This webpage was adapted from  "), a("ShinyCell", 
                                            href = "https://github.com/SGDDNB/ShinyCell",target="_blank")), 
    br(),br(),br(),br(),br() 
  ))) 


