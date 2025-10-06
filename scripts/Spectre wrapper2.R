# modified 9.23.2025
#renv::init()
#renv::snapshot() 

directory <- dirname(rstudioapi::getActiveDocumentContext()$path)       # Finds the directory where this script is located
setwd(directory)  
print(paste("Current working directory:", getwd()))

source("read.cytofFiles.R")
source("run.spectre2.R")
# User Customization ----------------------------------------------------

#define your parameters and file input location
run.spectre(phenok=50,
            metaFile="sample.details.csv",
            markerFile = "ORIGINAL MARKERS.csv",
            do.plot=TRUE,
            do.summary=TRUE,
            do.batchAlign=TRUE,
            do.fcsExport=TRUE, 
            flowType = "aurora",
            coFactor =2000,
            plot.against= "CD45RO_asinh",
            ref.ctrls=c("5-0002 ELUTE_CD3, CD4 _batch-2","5-0002 ELUTE_CD3, CD4 _batch-3"))
###    The defaults for cofactor: cytof = 5,    aurora = 2000,   flow = 200
#plot.against= "marker_asinh"
#phenok <- #put klevel 
#metaFile <- "sample.details.csv" #put file name
#markerFile <- "ORIGINAL MARKERS.csv" #put file n
#all the defaults are TRUE except for batch and sumamry table
#ref.ctrls need to be sample name as entered in sample.details.csv only when do.batchAlign is true
#ref.ctrls=c("5-0035_CD3, CD4 _citru batch1","5-0035_CD3, CD4 subset_batch2","5-0035_CD3, CD4 _citru batch2", "5-0035_CD3, CD4 subset_batch1")
#when there's batch alignment, CD4 grp name should better be CD4 BC 
#do.plot is for plots in addition to basic heatmap and clusters
#renv::snapshot() 
