#### Tasks ####

# 1. Set Working Directory
# Create a new folder on your computer "AI_Omics_Internship_2025".

## ** DONE **

# 2. Create Project Folder
# In RStudio, create a new project named "Module_I" in your "AI_Omics_Internship_2025" folder.

getwd()

## ** in the console: "D:/AI&biotech/AI_Omics_Internship_2025/Module_I" **

# Inside the project directory, create the following subfolders using R code:
# raw_data, clean_data, scripts, results or Tasks, plots etc

dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")


# ---------------------------------------------------------------------------
# 3. Download "patient_info.csv" dataset from GitHub repository

## ** DONE **

# load the dataset into your R environment

data <- read.csv(file.choose())

# Inspect the structure of the dataset using appropriate R functions

View(data)
str(data)

# Identify variables with incorrect or inconsistent data types.

## **I found that gender, diagnosis, and smoker variables 
# have the wrong data type. **

# Convert variables to appropriate data types where needed
 data$gender <- as.factor(data$gender)
data$gender <- factor(data$gender, levels = c("Male", "Female"))

data$diagnosis <- as.factor(data$diagnosis)
data$diagnosis <- factor(data$diagnosis, levels = c("Normal", "Cancer"))

data$smoker <- as.factor(data$smoker)
data$smoker <- factor(data$smoker, levels = c("Yes", "No"))

# Create a new variable for smoking status as a binary factor:
# 1 for "Yes", 0 for "No"

data$smoker_status <- ifelse(data$smoker == "Yes", 1, 0)
data
str(data)

# Save the cleaned dataset in your clean_data folder with the name patient_info_clean.csv
write.csv(data, file = "D:/AI&biotech/AI_Omics_Internship_2025/Module_I/clean_data/patient_info_clean.csv")
# Save your R script in your script folder with name "class_Ib"
## ** DONE **
# Upload "class_Ib" R script into your GitHub repository

save.image(file = "Shahenda_Dawoud_Class_Ib_Assignment.RData")
