#  Module I: Getting Started with R
# Class Ic: Syntax in R
# Includes:
# 1. Variables
#2. Comments
# 3. Keywords

#****
  #### 1. VARIABLES ####
  # Example: store a gene name
Gene <- "TP53"
# Retrieve the value in the console
Gene
print(Gene)

# Store numeric values in one variable (vector)
expression_levels <- c(3.4, 2.3, 5.7, 8.6)

# Import a CSV file as a variable

raw_data <- read.csv(file.choose())
raw_data


#### Rules for Variable Naming ####

# Must start with a letter
1protein <- 34  
#Error: unexpected symbol in "1protein"
# Invalid, as it starts with a number

protein1 <- 34    # Valid, number at the end

# No spaces in names

patient id <- "p01"  
#Error: unexpected symbol in "patient id"
# Invalid: contains space

# Use underscore (_) or dot (.) instead
patient_id <- "p01"
patient.id <- "p02"

# R is case-sensitive
Glucose_level <- 98
glucose_level <- 98  # Treated as a different variable

# Variables can be overwritten without warning
glucose_level <- c(110, 90, 120)


# Example of variable overwriting in a data

# Removing Columns

# First create a copy of raw_data to preserve the original dataset
copy_data <- raw_data

# Remove the 'smoker' column from the original raw_data
# This overwrites raw_data permanently, and the column cannot be restored
# unless a copy (like 'copy_data') was created beforehand
raw_data$smoker <- NULL

# Alternatively, create a new dataset without the last column
# This keeps the original 'copy_data' intact and stores the modified version in clean_data
clean_data <- copy_data[, -6]

####2. Comments####
# A hashtag before any writings, create a comment and it's very usefulto clarify what're you gonna code.
#R doesn't consider  the comment as a code, so it won't be run or stored in the environment

#####New Heading####
## Four hashtags before and after any comment create a heading

#To be continued#
save.image("Shahenda_Dawoud_class_Ic_Practice.RData")






