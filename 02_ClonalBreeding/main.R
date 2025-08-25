#Execute all
rm(list = ls())

reps=50
main_dir <- "C:/Users/Pheno/OneDrive - UNIVERSIDAD DE MURCIA/Escritorio/Almond_sim/almondsim/02_ClonalBreeding"
# Cambiar a la carpeta del esquema
setwd(main_dir)
source(file = "compatible_crosses.R")
pipeline<-TRUE


# ###### Burn-in #####

if (dir.exists("burn_in_folder")) {
  unlink("burn_in_folder", recursive = TRUE)  # Eliminar la carpeta y su contenido
}
dir.create("burn_in_folder")

setwd(file.path(main_dir, "00_Burn_in"))
source(file = "00RUNME.R")




### Run scenarios #####
setwd(file.path(main_dir))

# schemes <- c("01_phenotypicselection","02_pedigreeselection",
# "03_genomicselection",
# "10_phenotypicselection_sf_introgress_100", "10_phenotypicselection_sf_introgress_25",
# "11_genomicselection_sf_introgress_100","11_genomicselection_sf_introgress_25",
# "12_pedigreeselection_sf_introgress_25", "12_pedigreeselection_sf_introgress_100")

schemes <- c("01_PhenotypicSelection",
"03_GenomicSelection", "03_GenomicSelection_OCS_5","03_GenomicSelection_OCS_10","03_GenomicSelection_OCS_15",
"03_GenomicSelection_OCS_20","03_GenomicSelection_OCS_25","03_GenomicSelection_OCS_30")

# schemes <- c("03_genomicselection_ocs")

# schemes <- c("10_PhenotypicSelection_sf_introgress_100", "10_PhenotypicSelection_sf_introgress_25",
#              "11_GenomicSelection_sf_introgress_100","11_GenomicSelection_sf_introgress_25",
#              "12_PedigreeSelection_sf_introgress_25", "12_PedigreeSelection_sf_introgress_100")



#Removes previous results
if (file.exists("results_all_schemes.txt")) {
  file.remove("results_all_schemes.txt")  # Elimina el archivo
} else {
}


# Abre (o crea) el archivo de errores para escribir en él
error_file <- "errores.txt"
if (file.exists(error_file)) {
  file.remove(error_file)  # Elimina el archivo anterior si existe
}
file.create(error_file)



for (REP in 1:reps){
  
  for (scheme in schemes) {
    print(scheme)
    # Intenta ejecutar el bloque de código dentro del tryCatch
    tryCatch({
      # print(scheme)
      print(REP)
      main_dir <- "C:/Users/Pheno/OneDrive - UNIVERSIDAD DE MURCIA/Escritorio/Almond_sim/almondsim/02_ClonalBreeding"
      # Cambiar a la carpeta del esquema
      setwd(file.path(main_dir, scheme))
      
      keep_objects <- c("scheme", "main_dir","schemes", "REP","reps","stages","nBurnin","pipeline","randCrossGamSI",
                        "self_compatible_allele", "selectind_selfcompatible", "randCrossGamSI_selfcomp_compulsory")
      
      # Eliminar todo menos los objetos especificados
      rm(list = setdiff(ls(), keep_objects))
      

      load(paste0("../burn_in_folder/Burnin_", REP, ".RData"))
      # Ejecutar el script
      source(file = "00RUNME.R")
      
      # Combinar los resultados
      combined_df <- do.call(rbind, results)
      combined_df$scenario[combined_df$year >= nBurnin+1] <- scenarioName
      # Volver al directorio principal
      setwd("../")
      
      # Comprobar si el archivo ya existe y escribir los resultados
      if (file.exists("results_all_schemes.txt")) {
        write.table(combined_df, "results_all_schemes.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
      } else {
        write.table(combined_df, "results_all_schemes.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
      }
      

      
    }, error = function(e) {
      # Captura el error y escribe el mensaje en la consola y en el archivo de errores
      error_message <- paste("Error en el esquema:", scheme, "\n", 
                             "Mensaje de error:", conditionMessage(e), "\n\n")
      
      # Imprimir el error en consola
      cat(error_message)
      
      # Escribir el error en el archivo 'errores.txt'
      cat(error_message, file = "errores.txt", append = TRUE)
    })
  }

}




# Set the main directory
setwd(main_dir)
library(ggplot2)
library(dplyr)

# Read the data
data <- read.table("results_all_schemes.txt", header = TRUE, sep = "\t")

# Create a directory for saving the plots if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Calculate the mean of 'meanG' by 'scenario' and 'year'
data_summary <- data %>%
  group_by(scenario, year) %>%
  summarize(mean_meanG = mean(meanG, na.rm = TRUE))

# Create the line plot for 'meanG'
p_meanG <- ggplot(data_summary, aes(x = year, y = mean_meanG, color = scenario, group = scenario)) +
  geom_line(size = 1) +   # Thin lines
  labs(x = "Year", y = "Mean of MeanG", title = "MeanG by Scenario") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),  # 60% of 18
    axis.title = element_text(size = 8.4),             # 60% of 14
    axis.text = element_text(size = 7.2),              # 60% of 12
    legend.text = element_text(size = 7.2),            # 60% of 12
    legend.title = element_text(size = 8.4)            # 60% of 14
  )

# Export the plot as a JPG with higher resolution
jpeg("plots/meanG_plot.jpg", width = 2000, height = 1000, res = 300)
print(p_meanG)
dev.off()

# Calculate the mean of 'varG' by 'scenario' and 'year'
data_summary_varG <- data %>%
  group_by(scenario, year) %>%
  summarize(mean_varG = mean(varG, na.rm = TRUE))

# Create the line plot for 'varG'
p_varG <- ggplot(data_summary_varG, aes(x = year, y = mean_varG, color = scenario, group = scenario)) +
  geom_line(size = 1) +   # Thin lines
  labs(x = "Year", y = "Mean of VarG", title = "VarG by Scenario") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 8.4),
    axis.text = element_text(size = 7.2),
    legend.text = element_text(size = 7.2),
    legend.title = element_text(size = 8.4)
  )

# Export the plot as a JPG with higher resolution
jpeg("plots/varG_plot.jpg", width = 2000, height = 1000, res = 300)
print(p_varG)
dev.off()

# Calculate the mean of 'accSel' by 'scenario' and 'year'
data_summary_accSel <- data %>%
  group_by(scenario, year) %>%
  summarize(mean_accSel = mean(accSel, na.rm = TRUE))

# Create the line plot for 'accSel'
p_accSel <- ggplot(data_summary_accSel, aes(x = year, y = mean_accSel, color = scenario, group = scenario)) +
  geom_line(size = 1) +   # Thin lines
  labs(x = "Year", y = "Mean of AccSel", title = "Accuracy of Selection by Scenario") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 8.4),
    axis.text = element_text(size = 7.2),
    legend.text = element_text(size = 7.2),
    legend.title = element_text(size = 8.4)
  )

# Export the plot as a JPG with higher resolution
jpeg("plots/accSel_plot.jpg", width = 2000, height = 1000, res = 300)
print(p_accSel)
dev.off()

# Calculate the mean of 'HHI' ACT1 by 'scenario' and 'year'
data_summary_HHIact1 <- data %>%
  group_by(scenario, year) %>%
  summarize(mean_HHIact1 = mean(HHIACT1, na.rm = TRUE))

# Create the line plot for 'varG'
p_HHIact1 <- ggplot(data_summary_HHIact1, aes(x = year, y = mean_HHIact1, color = scenario, group = scenario)) +
  geom_line(size = 1) +   # Thin lines
  labs(x = "Year", y = "Mean of HHIact1", title = "HHIact1 by Scenario") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 8.4),
    axis.text = element_text(size = 7.2),
    legend.text = element_text(size = 7.2),
    legend.title = element_text(size = 8.4)
  )

# Export the plot as a JPG with higher resolution
jpeg("plots/HHIact1_plot.jpg", width = 2000, height = 1000, res = 300)
print(p_HHIact1)
dev.off()


# Calculate the mean of 'HHI' ECT1 by 'scenario' and 'year'
data_summary_HHIect1 <- data %>%
  group_by(scenario, year) %>%
  summarize(mean_HHIect1 = mean(HHIECT1, na.rm = TRUE))

# Create the line plot for ''
p_HHIect1 <- ggplot(data_summary_HHIect1, aes(x = year, y = mean_HHIect1, color = scenario, group = scenario)) +
  geom_line(size = 1) +   # Thin lines
  labs(x = "Year", y = "Mean of HHIect1", title = "HHIect1 by Scenario") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 8.4),
    axis.text = element_text(size = 7.2),
    legend.text = element_text(size = 7.2),
    legend.title = element_text(size = 8.4)
  )

# Export the plot as a JPG with higher resolution
jpeg("plots/HHIect1_plot.jpg", width = 2000, height = 1000, res = 300)
print(p_HHIect1)
dev.off()

# Calculate the mean of allelesSI by 'scenario' and 'year'
data_summary_allelesSI <- data %>%
  group_by(scenario, year) %>%
  summarize(mean_allelesSI = mean(allelesSI, na.rm = TRUE))

# Create the line plot for ''
p_allelesSI <- ggplot(data_summary_allelesSI, aes(x = year, y = mean_allelesSI, color = scenario, group = scenario)) +
  geom_line(size = 1) +   # Thin lines
  labs(x = "Year", y = "Mean of allelesSI", title = "allelesSI by Scenario") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 8.4),
    axis.text = element_text(size = 7.2),
    legend.text = element_text(size = 7.2),
    legend.title = element_text(size = 8.4)
  )

# Export the plot as a JPG with higher resolution
jpeg("plots/allelesSI_plot.jpg", width = 2000, height = 1000, res = 300)
print(p_allelesSI)
dev.off()

# Calculate the mean of allelesSI by 'scenario' and 'year'
data_summary_He_chr6 <- data %>%
  group_by(scenario, year) %>%
  summarize(mean_He_chr6 = mean(He_chr6, na.rm = TRUE))

# Create the line plot for 'He_chr6 '
p_He_chr6  <- ggplot(data_summary_He_chr6 , aes(x = year, y = mean_He_chr6 , color = scenario, group = scenario)) +
  geom_line(size = 1) +   # Thin lines
  labs(x = "Year", y = "Mean of He_chr6 ", title = "He_chr6  by Scenario") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 8.4),
    axis.text = element_text(size = 7.2),
    legend.text = element_text(size = 7.2),
    legend.title = element_text(size = 8.4)
  )

# Export the plot as a JPG with higher resolution
jpeg("plots/He_chr6 _plot.jpg", width = 2000, height = 1000, res = 300)
print(p_He_chr6)
dev.off()


#BOXPLOTS LAST YEAR

# Filter the data for the last year
last_year <- max(data$year)
data_last_year <- data %>% filter(year == last_year)

# Boxplot for 'meanG'
p_meanG_box <- ggplot(data_last_year, aes(x = scenario, y = meanG, fill = scenario)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot without outlier shapes for clarity
  geom_jitter(color = "black", size = 1, width = 0.2) + # Points for individual values
  labs(x = "Scenario", y = "MeanG", title = "Boxplot of MeanG (Last Year)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.text.x = element_text(size = 7.2, angle = 30, hjust = 1),
    axis.title = element_text(size = 8.4),
    legend.position = "none" # Remove legend since the x-axis already shows scenarios
  )

# Export the plot as a JPG
jpeg("plots/meanG_boxplot.jpg", width = 2000, height = 1000, res = 300)
print(p_meanG_box)
dev.off()

# Boxplot for 'varG'
p_varG_box <- ggplot(data_last_year, aes(x = scenario, y = varG, fill = scenario)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(color = "black", size = 1, width = 0.2) +
  labs(x = "Scenario", y = "VarG", title = "Boxplot of VarG (Last Year)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.text.x = element_text(size = 7.2, angle = 30, hjust = 1),
    axis.title = element_text(size = 8.4),
    legend.position = "none"
  )

# Export the plot as a JPG
jpeg("plots/varG_boxplot.jpg", width = 2000, height = 1000, res = 300)
print(p_varG_box)
dev.off()

# Boxplot for 'accSel'
p_accSel_box <- ggplot(data_last_year, aes(x = scenario, y = accSel, fill = scenario)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(color = "black", size = 1, width = 0.2) +
  labs(x = "Scenario", y = "Accuracy of Selection (accSel)", title = "Boxplot of accSel (Last Year)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.text.x = element_text(size = 7.2, angle = 30, hjust = 1),
    axis.title = element_text(size = 8.4),
    legend.position = "none"
  )

# Export the plot as a JPG
jpeg("plots/accSel_boxplot.jpg", width = 2000, height = 1000, res = 300)
print(p_accSel_box)
dev.off()



# Load the data
data <- read.table("results_all_schemes.txt", header = TRUE, sep = "\t")

# Define the breeding stages in order
stages <- c("Parents","ECT3", "ECT2", "ECT1", "ACT3", "ACT2", "ACT1", 
            "HPT4", "HPT3", "HPT2", "HPT1", "Seedlings", "F1")

# Define end year (customizable)
end_year <- 43 # Change this to whatever maximum year you want

# Helper: find the HHI column name for a given stage
get_hhi_col <- function(stage, df) {
  grep(paste0("^HHI.*", stage, "$"), names(df), value = TRUE)
}

# Function to summarize one year for a subset of data
summarize_year <- function(df_subset, target_year) {
  df_year <- subset(df_subset, year == target_year)
  
  meanG_vals <- sapply(stages, function(stage) {
    mean(df_year[[paste0("meanG_", stage)]], na.rm = TRUE)
  })
  varG_vals <- sapply(stages, function(stage) {
    mean(df_year[[paste0("varG_", stage)]], na.rm = TRUE)
  })
  hhi_vals <- sapply(stages, function(stage) {
    hhi_col <- get_hhi_col(stage, df_year)
    if (length(hhi_col) == 1) {
      mean(df_year[[hhi_col]], na.rm = TRUE)
    } else {
      NA
    }
  })
  
  data.frame(
    Stage = stages,
    meanG = meanG_vals,
    varG = varG_vals,
    HHI   = hhi_vals,
    row.names = NULL
  )
}

# Load Burn-in year 40 as baseline
burnin_40 <- summarize_year(subset(data, scenario == "Burn_in"), 40)
names(burnin_40)[2:4] <- c("meanG_Burnin40", "varG_Burnin40", "HHI_Burnin40")

# Define scenarios to compare
scenarios <- c("GS", "Pheno", "PedigreSelection")

# Loop over each scenario
for (sc in scenarios) {
  # Start with Burn-in year 40 data
  combined <- burnin_40
  
  # Loop over years 41 to end_year
  for (yr in 41:end_year) {
    df_tmp <- summarize_year(subset(data, scenario == sc), yr)
    
    # Rename columns to reflect scenario and year
    names(df_tmp)[2:4] <- paste0(c("meanG_", "varG_", "HHI_"), sc, yr)
    
    # Merge with accumulated table
    combined <- merge(combined, df_tmp, by = "Stage", all = TRUE)
  }
  
  # Write output for current scenario
  out_file <- paste0("summary_Burnin40_vs_", sc, "41_", end_year, ".txt")
  write.table(
    combined, file = out_file,
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}
