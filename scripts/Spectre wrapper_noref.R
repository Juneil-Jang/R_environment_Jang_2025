# modified 9.23.2025
#renv::init()
#renv::snapshot() 

## This version doesn't need the reference cells for batch. If you have, don't use this code

# Finds the directory where this script is located
directory <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(directory)  
print(paste("Current working directory:", getwd()))

## source
source("run.spectre_noref.R")

# User Customization ----------------------------------------------------

#define your parameters and file input location
run.spectre_noref(phenok=50,
             metaFile="sample.details.csv",
             markerFile = "ORIGINAL MARKERS.csv",
             meta_col = c('Sample', 'Group', 'Batch', "Donor"),
             do.plot=TRUE,
             do.summary=TRUE,
             do.batchAlign=TRUE,
             do.fcsExport=TRUE, 
             flowType = "aurora",
             coFactor =2000,
             plot.against= "CD45RO_asinh"
             )
###    The defaults for cofactor: cytof = 5,    aurora = 2000,   flow = 200
#plot.against= "marker_asinh" ex) CD45RO
#phenok <- #put klevel 
#metaFile <- "sample.details.csv" #put file name
#meta_col <- column names of metadata Type names of Sample, Group, Batch, Donor column
#if not typed, Sample, Group, Batch, Donor are used by default
#markerFile <- "ORIGINAL MARKERS.csv" #put file n
#all the defaults are TRUE
#do.plot is for plots in addition to basic heatmap and clusters

#renv::snapshot() 