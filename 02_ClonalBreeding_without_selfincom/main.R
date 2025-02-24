#Execute all

reps=20
main_dir <- "C:/Users/Pheno/OneDrive - UNIVERSIDAD DE MURCIA/Escritorio/Almond_sim/almondsim/02_ClonalBreeding_without_selfincom"
# Cambiar a la carpeta del esquema
setwd(main_dir)
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
# schemes <- c("01_PhenotypicSelection", "02_PedigreeSelection", "03_GenomicSelection",
#              "04_GenomicSelection_costeffective", "05_GenomicSelection_early_selection",
#              "06_GenomicSelection_early_selection_extreme")

# schemes <- c("01_PhenotypicSelection",  "03_GenomicSelection",
#              "04_GenomicSelection_costeffective", "05_GenomicSelection_early_selection",
#              "06_GenomicSelection_early_selection_extreme")



schemes <- c("01_PhenotypicSelection","03_GenomicSelection")

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
  
  # Iniciar la escritura en el archivo de errores

  
  for (scheme in schemes) {
    
    # Intenta ejecutar el bloque de código dentro del tryCatch
    tryCatch({
      
      main_dir <- "C:/Users/Pheno/OneDrive - UNIVERSIDAD DE MURCIA/Escritorio/Almond_sim/almondsim/02_ClonalBreeding_without_selfincom"
      # Cambiar a la carpeta del esquema
      setwd(file.path(main_dir, scheme))
      
      keep_objects <- c("scheme", "main_dir","schemes", "REP","reps","nBurnin","pipeline")
      
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


