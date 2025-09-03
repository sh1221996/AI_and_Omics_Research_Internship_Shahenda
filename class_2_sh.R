# Assignment 2
# --------------------------
# In this assignment you will work with the results of differential gene expression (DGE) analysis. 
#The analysis produces two key measures for each gene:

# log2FoldChange (log2FC): 
# Indicates the magnitude and direction of change in gene expression. 
# Positive values suggest higher expression(upregulated gene) in the experimental condition compared to control. 
# Negative values suggest lower expression (downregulated gene). 
# The absolute value reflects the strength of the change.

# Adjusted p-value (padj): 
# Represents the statistical significance of the observed difference, corrected for multiple testing. 
# A smaller value indicates stronger evidence that the observed difference is not due to chance.


#=================
# Write a function classify_gene() 

# that takes:
#   - logFC (log2FoldChange)
#   - padj  (adjusted p-value)

# and returns:
#   - "Upregulated" if log2FC > 1 and padj < 0.05
#   - "Downregulated" if log2FC < -1 and padj < 0.05
#   - "Not_Significant" otherwise

classify_gene <- function(logFC, padj) {
  ifelse(logFC > 1 & padj < 0.5, "Upregulated", 
  ifelse(logFC < -1 & padj < 0.05, "Downregulated", "Not_Significant"))
}

#===============
input_dir <- "row.data"
output_dir <- "results"
 if(!dir.exists(output_dir)) {
   dir.create(output_dir)
 }
#provide the files to be processed

files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")
 
#prepare a list for later results

 result_list <- list()
 
 
# Then:
 #   - Apply it in a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)
 #   - Replace missing padj values with 1

 #=== create the loop function
 
 for (file_names in files_to_process) {
  cat("\nProcessing", file_names, "\n")
   
  #create file input path 
  file_input_path <- file.path(input_dir, file_names)
  
  #importing the datasets 
  
  Gn_exp_data <- read.csv(file_input_path, header = TRUE) 
   cat("files imported. checking for missing values...\n")
  
  #handling missing values in padj  
  if("padj" %in% names(Gn_exp_data)) {
    missing_count <- sum(is.na(Gn_exp_data$padj))
    cat("Missing values in 'padj':", missing_count, "\n")
    # replace missing values
    Gn_exp_data$padj[is.na(Gn_exp_data$padj)] <- 1 
  }
  } 
   
  
#   - Add a new column 'status'
 
 Gn_exp_data$status <- classify_gene(Gn_exp_data$logFC, Gn_exp_data$padj)
 cat("Classification is completed")
#   - Save processed files into Results folder
 result_list[[file_names]] <- Gn_exp_data
#   - Print summary counts of significant, upregulated, and downregulated genes
#   - Use table() for summaries
 print(table(Gn_exp_data$status))

 output_file_path <- file.path(output_dir, paste0("classified_gn_exps", file_names))
 
 #save the result as csv file
 write.csv(Gn_exp_data, output_file_path, row.names = FALSE)

#======================
 
# Data Availability
# The input files are available in the GitHub repository:
#      DEGs_Data_1.csv
#      DEGs_Data_2.csv
# Each file contains three columns: 
 # Gene_Id	
 # padj	
 # logFC
 
 
save.image(file = "Shahenda_Dawoud_class2_assignment.RData")

