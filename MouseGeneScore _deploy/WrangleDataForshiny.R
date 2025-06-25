library(Seurat)
library(ShinyCell)
library(scCustomize)

dir <- "~/Documents/AstrocytePaper"

seu = readRDS(file.path(dir,"Astrocyteintegration_AmbientRemoved_filtered_noneuron.RDS"))


scConf = createConfig(seu)

scConf = modColours(scConf, meta.to.mod = "ClusterNames", 
                    new.colours= c("dodgerblue", "goldenrod1", "red", "darkgreen"))


scConf = modColours(scConf, meta.to.mod = "finalClusters", 
                    new.colours= c("dodgerblue", "goldenrod1", "red", "darkgreen"))


scConf = modColours(scConf, meta.to.mod = "Disease", 
                    new.colours= c("darkgreen","lightgreen", "purple","plum1"))



makeShinyApp(seu, scConf,
             shiny.title = "Mouse Interactive Shiny",
             default.gene1 = "Gfap", default.gene2 = "Aldoc",
             default.multigene =  c("Gpc5","Lsamp" ,"Ntm","Wdr17","Nrxn1",
                                    "Gfap","Id3","Mt2","Id4","Apoe","Cst3",
                                    "Aldoc","Mt1","Ckb","Meg3","Dlg2","Nrg3",
                                    "Nkain2", "Nrxn3")) 

