##Primero hay que instalar las librerias necesarias y depues cargarlos 
install.packages("readxl")
install.packages("dplyr")
install.packages("data.table")
install.packages("fgsea")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")
a
BiocManager::install("fgsea", force = TRUE)
a

##Cargar las librearias
library(readxl)
library(dplyr)

# Define la ruta del archivo Excel
ruta_excel <- "D:\\Carpetas\\Experimentos\\Celulas Endoteliales\\RNA\\RNAseq\\Resultados\\Nuevo analisis Jose\\All compare\\Tablaall _compare_LOG2.xlsx"

# Define la ruta donde guardar los archivos .rnk
ruta_guardado <- "D:\\Carpetas\\Experimentos\\Celulas Endoteliales\\RNA\\RNAseq\\Resultados\\Nuevo analisis Jose\\All compare\\ALL_COMPARE_filtrado.xlsx"

# Lee el archivo Excel
datos <- read_excel(ruta_excel)

# Definir los nombres de las columnas con los log2FC
comparaciones <- c(
  "DKOvsWTP_log2FoldChange",
  "KOMvsWTP_log2FoldChange",
  "KOSvsWTP_log2FoldChange",
  "WTPvsWTG_log2FoldChange",
  "KOPvsWTG_log2FoldChange",
  "KOMvsKOP_log2FoldChange",
  "KOSvsKOP_log2FoldChange",
  "DKOvsKOP_log2FoldChange",
  "KOSvsKOM_log2FoldChange",
  "DKOvsKOM_log2FoldChange",
  "DKOvsKOS_log2FoldChange",
  "KOMvsWTG_log2FoldChange",
  "KOSvsWTG_log2FoldChange",
  "DKOvsWTG_log2FoldChange",
  "KOPvsWTP_log2FoldChange"
)

# Genera los archivos .rnk
for (comp in comparaciones) {
  rnk <- datos %>%
    select(gene_name, !!sym(comp)) %>%
    filter(!is.na(!!sym(comp))) %>%
    arrange(desc(!!sym(comp)))
  
  write.table(rnk,
              file = paste0(ruta_guardado, comp, ".rnk"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}
