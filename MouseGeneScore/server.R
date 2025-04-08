#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
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
library(ggplot2) 
library(ggrepel) 
library(hdf5r) 
library(ggdendro) 
library(gridExtra) 
sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")
sc1gene = readRDS("sc1gene.rds")
sc1meta = readRDS("sc1meta.rds")



### Useful stuff 
# Colour palette 
cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"), 
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)], 
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C", 
               "#2C728E","#3B528B","#472D7B","#440154")) 
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple") 

# Panel sizes 
pList = c("400px", "700px", "800px") 
names(pList) = c("Small", "Medium", "Large") 
pList2 = c("500px", "700px", "900px") 
names(pList2) = c("Small", "Medium", "Large") 
pList3 = c("800px", "1000px", "1200px") 
names(pList3) = c("Small", "Medium", "Large") 
sList = c(18,24,30) 
names(sList) = c("Small", "Medium", "Large") 
lList = c(5,6,7) 
names(lList) = c("Small", "Medium", "Large") 

# Function to extract legend 
g_legend <- function(a.gplot){  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")  
  legend <- tmp$grobs[[leg]]  
  legend 
}  

# Plot theme 
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){ 
  oupTheme = theme( 
    text =             element_text(size = base_size, family = "Helvetica"), 
    panel.background = element_rect(fill = "white", colour = NA), 
    axis.line =   element_line(colour = "black"), 
    axis.ticks =  element_line(colour = "black", linewidth = base_size / 20), 
    axis.title =  element_text(face = "bold"), 
    axis.text =   element_text(size = base_size), 
    axis.text.x = element_text(angle = Xang, hjust = XjusH), 
    legend.position = "bottom", 
    legend.key =      element_rect(colour = NA, fill = NA) 
  ) 
  if(!XYval){ 
    oupTheme = oupTheme + theme( 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  } 
  return(oupTheme) 
} 


# Define server logic required to draw a histogram

searchGeneSignature <- reactive({
  input$searchButton # Make reactive response on button click
  isolate({
    genes <- strsplit(input$gene_sig_input, "[,;\\n]+")[[1]] %>% trimws()
    genes <- genes[genes %in% names(sc1gene)]
    if (length(genes) == 0) {
      shiny::showNotification("No valid gene signatures found!", type = "error")
    }
    return(genes)
  })
})

scGeneList <- function(inp, inpGene){ 
  geneList = data.table(gene = unique(trimws(strsplit(inp, ",|;|
")[[1]])), 
                        present = TRUE) 
  geneList[!gene %in% names(inpGene)]$present = FALSE 
  return(geneList) 
} 



scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                       inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, 
                       inpcols, inpfsz, save = FALSE){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Identify genes that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!")) 
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
  
  # Prepare ggData 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData = data.table() 
  for(iGene in geneList$gene){ 
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE] 
    colnames(tmp) = c("sampleID", "sub") 
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]] 
    tmp$geneName = iGene 
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=))) 
    ggData = rbindlist(list(ggData, tmp)) 
  } 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!")) 
  
  # Aggregate 
  ggData$val = expm1(ggData$val) 
  ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(sampleID)), 
                  by = c("geneName", "grpBy")] 
  ggData$val = log1p(ggData$val) 
  
  # Scale if required 
  colRange = range(ggData$val) 
  if(inpScl){ 
    ggData[, val:= scale(val), keyby = "geneName"] 
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
  } 
  
  # hclust row/col if necessary 
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val") 
  tmp = ggMat$geneName 
  ggMat = as.matrix(ggMat[, -1]) 
  rownames(ggMat) = tmp 
  if(inpRow){ 
    hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat)))) 
    ggRow = ggplot() + coord_flip() + 
      geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)), 
                         labels = unique(ggData$grpBy), expand = c(0, 0)) + 
      scale_x_continuous(breaks = seq_along(hcRow$labels$label), 
                         labels = hcRow$labels$label, expand = c(0, 0.5)) + 
      sctheme(base_size = sList[inpfsz]) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.y = element_blank(), 
            axis.text.x = element_text(color="white", angle = 45, hjust = 1)) 
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label) 
  } else { 
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene)) 
  } 
  if(inpCol){ 
    hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat))))) 
    ggCol = ggplot() + 
      geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_x_continuous(breaks = seq_along(hcCol$labels$label), 
                         labels = hcCol$labels$label, expand = c(0.05, 0)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)), 
                         labels = unique(ggData$geneName), expand=c(0,0)) + 
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(), 
            axis.text.y = element_text(color = "white")) 
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label) 
  } 
  
  # Actual plot according to plottype 
  if(inpPlt == "Bubbleplot"){ 
    # Bubbleplot 
    ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) + 
      geom_point() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_size_continuous("proportion", range = c(0, 8), 
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) + 
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(color = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank(), legend.box = "vertical") 
  } else { 
    # Heatmap 
    ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) + 
      geom_tile() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(fill = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank()) 
  } 
  
  # Final tidy 
  ggLeg = g_legend(ggOut) 
  ggOut = ggOut + theme(legend.position = "none") 
  if(!save){ 
    if(inpRow & inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                   layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                   layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      grid.arrange(ggOut, ggLeg, heights = c(7,2),  
                   layout_matrix = rbind(c(1),c(2)))  
    }  
  } else { 
    if(inpRow & inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                  layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                  layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      arrangeGrob(ggOut, ggLeg, heights = c(7,2),  
                  layout_matrix = rbind(c(1),c(2)))  
    }  
  } 
  return(ggOut) 
} 

scSigSearch <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                       inpsub1, inpsub2, inpH5, inpGene,inpord,
                       inpcols, inpfsz,inpsiz, save = FALSE){ 
  
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Identify genes that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!")) 
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
  # New prepare
  ggData = inpMeta[, c(inpConf[UI == sc1def$dimred[1]]$ID, inpConf[UI == sc1def$dimred[2]]$ID, 
                       inpConf[UI == inpsub1]$ID,inpConf[UI == inpGrp]$ID),  
                   with = FALSE]
  colnames(ggData) = c("X", "Y", "sub","group")
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y))
  # Read gene signature expression levels
  h5file <- H5File$new(inpH5, mode = "r")
  h5data <- h5file[["grp"]][["data"]]
  
  signature_expr <- sapply(geneList$gene, function(gene_name) {
    gene_id <- inpGene[gene_name]
    if (is.na(gene_id)) return(rep(NA, nrow(ggData)))
    expr_values <- h5data$read(args = list(gene_id, quote(expr=)))
    pmax(0, expr_values) # Clamp to non-negative values
  })
  
  h5file$close_all()
  valid_exprs <- signature_expr[, colSums(is.na(signature_expr)) == 0]
  if (ncol(valid_exprs) > 0) ggData$val <- rowMeans(valid_exprs, na.rm = TRUE)
  # Handle subsetting
  bgCells = FALSE
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)) {
    bgCells = TRUE
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2]
  }
  
  # Plotting
   # Actual plot according to plottype 
  if(inpPlt == "UMAP"){ 
    # Ordering depending on chosen method
    if (inpord == "Max-1st"){ 
      ggData = ggData[order(val, decreasing = FALSE)]
    } else if (inpord == "Min-1st") {
        ggData = ggData[order(val,decreasing = TRUE)] 
    } else if (inpord == "Random") {
          ggData = ggData[sample(nrow(ggData))] }
    #UMAP 
    ggOut = ggplot(ggData, aes(X, Y, color = val))
    if (bgCells){
      ggOut = ggOut + geom_point(data = ggData2, color = "snow2", size = 4,
                                 shape = 16)
    }
    ggOut = ggOut + geom_point(size = 4, shape = 16) +
      xlab("UMAP1") + ylab("UMAP2") +
      sctheme(base_size = sList[4], XYval = 12) +
      scale_color_gradientn("Signature Expression", colours = cList[[inpcols]]) +
      guides(color = guide_colorbar(barwidth = 15)) +coord_fixed(ratio = 1)
    
  } else { 
    # Bar plot 
    if(!is.na(inpConf[UI == inpGrp]$fCL)){ 
      ggCol = strsplit(inpConf[UI == inpGrp]$fCL, "\\|")[[1]] 
      names(ggCol) = levels(ggData$group) 
      ggLvl = levels(ggData$group)[levels(ggData$group) %in% unique(ggData$group)] 
      ggData$group = factor(ggData$group, levels = ggLvl) 
      ggCol = ggCol[ggLvl] 
    } 
    ggOut <- ggplot(ggData, aes(x = group, y = val, fill = group)) +
      scale_fill_manual("", values = ggCol)+
      geom_boxplot()+
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1)+
      theme(legend.position = "none") + ylab("Signature Expression")
  }
  
  return(ggOut) 
} 

# Start Server-------
shinyServer(function(input, output, session) { 
  
  ## Bubble plot UI ----------
  output$sc1d1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1d1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1d1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1d1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1d1oupTxt <- renderUI({ 
    geneList = scGeneList(input$sc1d1inp, sc1gene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$sc1d1oup <- renderPlot({ 
    scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
               input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
               input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
               input$sc1d1cols, input$sc1d1fsz) 
  }) 
  output$sc1d1oup.ui <- renderUI({ 
    plotOutput("sc1d1oup", height = pList3[input$sc1d1psz]) 
  }) 
  output$sc1d1oup.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1d1plt,"_",input$sc1d1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1d1oup.h, width = input$sc1d1oup.w, 
      plot = scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                        input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                        input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                        input$sc1d1cols, input$sc1d1fsz, save = TRUE) ) 
    }) 
  output$sc1d1oup.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1d1plt,"_",input$sc1d1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1d1oup.h, width = input$sc1d1oup.w, 
      plot = scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                        input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                        input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                        input$sc1d1cols, input$sc1d1fsz, save = TRUE) ) 
    }) 
  ## End Bubble plot --------------
  
  ## Signature Search UI----------------
  output$sigsub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sigsub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sigsub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sigsub1non, { 
    sub = strsplit(sc1conf[UI == input$sigsub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sigsub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sigsub1all, { 
    sub = strsplit(sc1conf[UI == input$sigsub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sigsub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sigoupTxt <- renderUI({ 
    geneList = scGeneList(input$siginp, sc1gene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oupsig = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oupsig = paste0(oupsig, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oupsig) 
    } 
  }) 
  
  # Run search-----
  output$sigoup <- renderPlot({ 
    scSigSearch(sc1conf, sc1meta, input$siginp, input$siggrp, input$sigplt, 
               input$sigsub1, input$sigsub2, "sc1gexpr.h5", sc1gene,input$sigord1,
               input$sigcols, input$sigfsz,input$sigpsz) 
  }) 
  output$sigoup.ui <- renderUI({ 
    plotOutput("sigoup", height = pList2[input$sigpsz]) 
  }) 
  output$sigoup.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sigplt,"_",input$siggrp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sigoup.h, width = input$sc1d1oup.w, 
      plot = scSigSearch(sc1conf, sc1meta, input$siginp, input$siggrp, input$sigplt, 
                         input$sigsub1, input$sigsub2, "sc1gexpr.h5", sc1gene,input$sigord1,
                         input$sigcols, input$sigfsz,input$sigpsz, save = TRUE) ) 
    }) 
  output$sigoup.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sigplt,"_",input$siggrp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sigoup.h, width = input$sigoup.w, 
      plot = scSigSearch(sc1conf, sc1meta, input$siginp, input$siggrp, input$sigplt, 
                         input$sigsub1, input$sigsub2, "sc1gexpr.h5", sc1gene,input$sigord1,
                         input$sigcols, input$sigfsz,input$sigpsz, save = TRUE) ) 
    }) 
# End Sig Search ----------
  
})