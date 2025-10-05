####Assignment for Module_II Class 3#####
#On RStudio, I created a new project "Module_II"
#First I will create the folders to store my data, script, and results:

dir.create("raw_data")
dir.create("scripts")                        
dir.create("results")

#In this assignment, I will be doing Microarray data analysis
##Installing the required packages to do so:
 #1) check if BiocManager is installed:

if(!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
#2) Use the BiocManager to install the required packages (GEOquery, Affy, Arrayqualitymetrics):
BiocManager::install(c("GEOquery", "affy", "arrayQualityMetrics"))

#3) Install dplyr from CRAN, for data manupulation:
install.packages("dplyr")

#4) Load those packages to be able to use their functions:

library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(dplyr)

#### Data retrieval ####
##I will be using the data set of this study (GSE118916)
#Getting the series matrix data (Processed data), "NOT MANDATORY FOR THIS WORK"
series_data <- getGEO("GSE118916", GSEMatrix = TRUE)

##Let's extract and explore the series matrix data:

# Extract expression data matrix (genes/probes × samples)
# Rows corresponds to probes and columns corresponds to samples
expression_data <- exprs(series_data$GSE118916_series_matrix.txt.gz)


# Extract feature (probe annotation) data
# Rows corresponds to probes and columns corresponds to samples
feature_data <-  fData(series_data$GSE118916_series_matrix.txt.gz)


# Extract phenotype (sample metadata) data
# Rows corresponds to samples and columns corresponds to probes
phenotype_data <-  pData(series_data$GSE118916_series_matrix.txt.gz)

##I will use this phenotype_data to rename the column names of my processed/filtered data (Later)

# Check missing values in sample annotation
sum(is.na(phenotype_data$source_name_ch1)) 

##Raw Data:
#I've already downloaded the raw data as a tar file (compressed), I moved it to raw_data directory

#Untar the raw data into CEL files

untar("raw_data/GSE118916_RAW.tar", exdir = "raw_data/CEL_files")
#Read the CEL files as AffyBatch object:
raw_data <- ReadAffy(celfile.path = "raw_data/CEL_files")
#By displaying raw_data, note that the annotation type is 'primeview'
#we'll need this later so that we use the specific annotation package for mapping 
#probe ID to genes


####Quality check before pre-processing ####

#We can do QC using multiple ways including:
#1. Box plot
#2. MA Plots
#3. Heatmaps and distance matrices
#4. PCA
##We can alternatively use one method that combine them all (arrayQualityMetrics):
#This will create an interactive HTML report that summarizes the QC results
 
arrayQualityMetrics(expressionset = raw_data, 
                    outdir = "results/QC_raw_data",
                    force = TRUE,
                    do.logtransform = TRUE)


#From the HTML report, shown that there are 3 outliers as follows:
#2	GSM3351221_A3607_PrimeView_.CEL.gz	(x	x) [2 outliers in array/sample "2"]
#17	GSM3351236_A3608_PrimeView_.CEL.gz		(x	) [1 outlier in array/smple "17]

####Data Normalization####

##Robust Multi-array Average "RMA':
normalized_data <- rma(raw_data)

#QC after normalization (excluding the log transformation):
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "results/QC_Normalized_data",
                    force = TRUE)

#No outliers detected after normalization

###Extracting the expression values from normalized data:

processed_data <- as.data.frame(exprs(normalized_data))
dim(processed_data)  # This command to check the dimensions of data frame (nrows:49495 and ncolmns:30)

##In this expression data frame, each row is the probe and each column is the sample 
#This is our main data set that we'll use for further downstream analyses

##BUT before proceeding with it, we should filter it as some probe values may show a little variation
#This may cause problems in the downstream analysis or the ML task
#By removing them, we focus on genes that are biologically relevant

###Filter low-variant transcripts ###

#1. Calculating median intensity per probe across samples:
#Gives us an overall sense of expression level of each probe

Row_median <- rowMedians(as.matrix(processed_data))
Row_median

#Now plot the distribution of the median intensities of probes:
hist(Row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")
threshold <- 4
abline(v = threshold, col = "red", lwd = 2) 

#Keep only probes above the threshold:
good_probes <-  Row_median > threshold 
filtered_data <- processed_data[good_probes,]

#Renaming the filtered expression data with the sample metadata:
colnames(filtered_data) <- rownames(phenotype_data)
#Now the sample names (filtered_data column names) are clearer and less bulky
## Overwrite the processed data with the filtered data (that are now ready for further analysis)
processed_data <- filtered_data
##Now our processed data has only the good probes with transcripts above the threshold


###Phenotype data preparation###
#This contains the sample-level metadata
#We will focus only on the type of tissues of the samples, whether the tissue is normal or has cancer
#The column source_name_ch1 has this information
#we need to define experimental groups for statistical analysis"

#If we check the data type of the values in this column:
class(phenotype_data$source_name_ch1)
#We will find it "Character"

#Change it into factor and define the experimental groups (Normal, Cancer):
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("human normal stomach tissue","human gastric tumor"),
                 labels = c("Normal", "Cancer")) #For more convenience
class(groups)
levels(groups)

save.image(file = "Shahenda_Dawoud_class3B_assignment.RData")






