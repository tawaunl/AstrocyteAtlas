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

sigsearchdef <- c("TNC","C3","CP","GFAP","NEAT1","CCL2","GAP43","SLC8A3")
# Define UI for application that draws a histogram
shinyUI(fluidPage(
    tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))),
    list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))),
    titlePanel("Human Downsampled Interactive Shiny"),
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
            tags$h4(
                tags$span(
                    "PLEASE NOTE: DATA IS DOWNSAMPLED. To allow for fluid plotting, data has been downsampled",
                    "maintaining proper proportions of cells within each study, donor, cluster, and diagonsis",  # Your desired text
                    style = "color: red; font-weight: bold;"  # Inline CSS for red, bold text
                ),
                align = "center"  # Optional alignment
            ),
            br(),br(),
            fixedRow(
                column(
                    7,
                    h4("Dimension Reduction"),
                    selectInput(
                        "sc1a2drX",
                        "X-axis:",
                        choices = sc1conf[dimred == TRUE]$UI,
                        selected = sc1def$dimred[1],
                        width = "200px" # Fixed width
                    ),
                    selectInput(
                        "sc1a2drY",
                        "Y-axis:",
                        choices = sc1conf[dimred == TRUE]$UI,
                        selected = sc1def$dimred[2],
                        width = "200px" # Fixed width
                    )
                ),
                column(5,
                       h4("Cell information"),
                       selectInput("sc1a2inp2", "Cell information:",
                                   choices = sc1conf$UI,
                                   selected = sc1def$meta2) %>%
                           helper(type = "inline", size = "m", fade = TRUE,
                                  title = "Cell information to colour cells by",
                                  content = c("Select cell information to colour cells",
                                              "Categorical covariates have a fixed colour palette",
                                              paste0("- Continuous covariates are coloured in a ",
                                                     "Blue-Yellow-Red colour scheme, which can be ",
                                                     "changed in the plot controls"))),
                       actionButton("sc1a2tog2", "Toggle plot controls"),
                       conditionalPanel(
                           condition = "input.sc1a2tog2 % 2 == 1",
                           radioButtons("sc1a2col2", "Colour (Continuous data):",
                                        choices = c("White-Red", "Blue-Yellow-Red",
                                                    "Yellow-Green-Purple"),
                                        selected = "Blue-Yellow-Red"),
                           radioButtons("sc1a2ord2", "Plot order:",
                                        choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                        selected = "Original", inline = TRUE),
                           checkboxInput("sc1a2lab2", "Show cell info labels", value = FALSE)
                       )
                )),
            fixedRow(
                column(
                    3, style="border-right: 2px solid black; padding-left: 10px;",
                    textAreaInput("siginp", HTML("List of gene names <br />
                                          (Max 50 genes, separated <br />
                                           by , or ; or newline):"),
                                  height = "200px", width = "200px",
                                  value = paste0(sigsearchdef, collapse = ", ")) %>%
                        helper(type = "inline", size = "m", fade = TRUE,
                               title = "List of genes to create signature",
                               content = c("Input genes",
                                           "- Maximum 50 genes",
                                           "- Genes should be separated by comma, semicolon or newline")),
                    selectInput("siggrp", "Group by (Bar Plot only):",
                                choices = sc1conf[grp == TRUE]$UI,
                                selected = sc1conf[grp == TRUE]$UI[10],width = "200px") %>%
                        helper(type = "inline", size = "m", fade = TRUE,
                               title = "Cell information to group cells by",
                               content = c("Select categorical cell information to group cells by",
                                           "- Single cells are grouped by this categorical covariate",
                                           "- Plotted as the X-axis of the Bar Plot")),
                    radioButtons("sigplt", "Plot type:",
                                 choices = c("DimRed", "BarPlot"),
                                 selected = "DimRed", inline = TRUE),
                    br(),
                    actionButton("sigtogL", "Toggle to subset cells"),
                    conditionalPanel(
                        condition = "input.sigtogL % 2 == 1",
                        selectInput("sigsub1", "Cell information to subset:",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp2),
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
                column(4, h4(htmlOutput("sigoupTxt")),
                       plotOutput("sigoup",height = "800px"),
                       downloadButton("sigoup.pdf", "Download PDF",),
                       downloadButton("sigoup.png", "Download PNG"), br(),
                       div(style="display:inline-block",
                           numericInput("sigoup.h", "PDF / PNG height:", width = "100px",
                                        min = 4, max = 20, value = 10, step = 0.5)),
                       div(style="display:inline-block",
                           numericInput("sigoup.w", "PDF / PNG width:", width = "100px",
                                        min = 4, max = 20, value = 10, step = 0.5))
                ),# End of column (6 space)
                column(
                    4, style="border-left: 2px solid black; padding-right: 10px;",
                    column(12,
                           plotOutput("sc1a2oup2",height = "800px"),
                           downloadButton("sc1a2oup2.pdf", "Download PDF"),
                           downloadButton("sc1a2oup2.png", "Download PNG"), br(),
                           div(style="display:inline-block",
                               numericInput("sc1a2oup2.h", "PDF / PNG height:", width = "100px",
                                            min = 4, max = 20, value = 6, step = 0.5)),
                           div(style="display:inline-block",
                               numericInput("sc1a2oup2.w", "PDF / PNG width:", width = "100px",
                                            min = 4, max = 20, value = 8, step = 0.5))
                    )  # End of column (6 space)
                ) # End of fluidRow (4 space)
            )),    # End of tab (2 space)
        # Bubble plot ----------------
        tabPanel(
            HTML("Bubbleplot / Heatmap"),
            h4("Gene expression bubbleplot / heatmap"),
            "In this tab, users can visualise the gene expression patterns of ",
            "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(),
            "The normalised expression are averaged, log-transformed and then plotted.",
            br(),br(),
            tags$h4(
                tags$span(
                    "PLEASE NOTE: DATA IS DOWNSAMPLED. To allow for fluid plotting, data has been downsampled",
                    "maintaining proper proportions of cells within each study, donor, cluster, and diagonsis",  # Your desired text
                    style = "color: red; font-weight: bold;"  # Inline CSS for red, bold text
                ),
                align = "center"  # Optional alignment
            ),
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
                                selected = sc1conf[grp == TRUE]$UI[10]) %>%
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
                                    selected = sc1def$grp2),
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


