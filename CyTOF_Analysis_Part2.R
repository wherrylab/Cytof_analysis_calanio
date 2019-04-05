## This script details how we process Phenograph analysis on CyTOF data using the cytofkit package developed by others previously

library(stringr)
library(cluster)
#Install cytofkit if not previously installed
#source("https://bioconductor.org/biocLite.R")
#biocLite("cytofkit")
library(cytofkit)
library(flowCore)
##### Run Phenograph analysis on multiple samples ------------------------------------------------------------------------------------------------------------------------------------------------

# Select a folder containing your fcs files to analyse (the subsetted one exported from FlowJo)
dir <- "."
# Identify files
files <- list.files(dir ,pattern='.fcs$', full=TRUE)
isFCSfile(files)
# Extract flow parameters list
fcs_file=read.FCS(filename=files[1], transformation=NULL, which.lines=NULL,
                  alter.names=FALSE, column.pattern=NULL, invert.pattern = FALSE,
                  decades=0, ncdf = FALSE, min.limit=NULL, truncate_max_range = TRUE, dataet=NULL, emptyValue=TRUE)

#Create a parameter text file after cleaning up names of parameters in one of the fcs file
p=toString(parameters(fcs_file)$name)
p
p=strsplit(p, ",")[[1]]
pceil
q=toString(parameters(fcs_file)$desc)
q
q=strsplit(q, ",")[[1]]
q
r=paste(p, q, sep='<')
r
r=paste(r, '>')
r
r=str_replace_all(string=r, pattern=" ", repl="")
r
write(r,"parameters_T.txt")

# Open the generated txt file with a text reader, and keep only markers of interest for analysis
#then
paraFile <- list.files(dir, pattern='.txt$', full=TRUE)
# Identify parameters to be included in the analysis
parameters <- as.character(read.table(paraFile, header = T)[,1])
parameters

# Run phenograph function
cytofkit(fcsFiles = files, markers = parameters[1:10],
         projectName = "cytofkit_Project1",
         transformMethod = c("cytofAsinh"),
         mergeMethod = c("min"),
         dimReductionMethod = c("tsne"),
         clusterMethods = c("Rphenograph"),
         visualizationMethods = c("tsne"),
         progressionMethod = c("NULL"),
         resultDir = getwd(), 
         saveResults = TRUE,
         saveObject = TRUE, openShinyAPP = TRUE)

# This will generate a folder that you can store in a general CyTOF_cytofkit folder
# The processing of generated data is explained in CyTOF_Analysis_Part3

