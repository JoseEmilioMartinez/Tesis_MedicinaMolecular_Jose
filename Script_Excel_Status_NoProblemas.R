
##Calcular promedio FPKM por grupo para cada gen.
##Definir un estado simple: ejemplo, si el promedio en un grupo es 2x mayor que en otro, marcar "aumentado" o "disminuido".


# Cargar librerías necesarias
library(readxl)
library(dplyr)
library(openxlsx)

# Ruta al archivo
ruta_archivo <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis nuevo/GCKIII/Heatmap GCKIII genes.xlsx"

# Cargar datos
genes_heatmap <- read_excel(ruta_archivo)

# Asegurar que gene_name es carácter
genes_heatmap$gene_name <- as.character(genes_heatmap$gene_name)

# Identificar columnas por grupo
DKO_cols <- grep("^DKO", colnames(genes_heatmap), value = TRUE)
WTG_cols <- grep("^WTG", colnames(genes_heatmap), value = TRUE)
KOS_cols <- grep("^KOS", colnames(genes_heatmap), value = TRUE)
KOM_cols <- grep("^KOM", colnames(genes_heatmap), value = TRUE)

# Limpiar y convertir columnas a numérico
clean_numeric <- function(x){
  x <- gsub(",", ".", as.character(x))   # reemplazar coma decimal por punto
  x <- gsub("#¡NUM!", NA, x)             # reemplazar errores de Excel por NA
  as.numeric(x)
}

genes_heatmap[DKO_cols] <- lapply(genes_heatmap[DKO_cols], clean_numeric)
genes_heatmap[WTG_cols] <- lapply(genes_heatmap[WTG_cols], clean_numeric)
genes_heatmap[KOS_cols] <- lapply(genes_heatmap[KOS_cols], clean_numeric)
genes_heatmap[KOM_cols] <- lapply(genes_heatmap[KOM_cols], clean_numeric)

# Convertir columnas FPKM a numéricas (en caso de que no lo sean)
genes_heatmap[WTG_cols] <- lapply(genes_heatmap[WTG_cols], as.numeric)
genes_heatmap[WTP_cols] <- lapply(genes_heatmap[WTP_cols], as.numeric)
genes_heatmap[KOP_cols] <- lapply(genes_heatmap[KOP_cols], as.numeric)
genes_heatmap[DKO_cols] <- lapply(genes_heatmap[DKO_cols], as.numeric)

# Calcular promedio FPKM por grupo
genes_heatmap <- genes_heatmap %>%
  mutate(
    avg_DKO = rowMeans(select(., all_of(DKO_cols)), na.rm = TRUE),
    avg_WTG = rowMeans(select(., all_of(WTG_cols)), na.rm = TRUE),
    avg_KOS = rowMeans(select(., all_of(KOS_cols)), na.rm = TRUE),
    avg_KOM = rowMeans(select(., all_of(KOM_cols)), na.rm = TRUE)
  )

# Definir estados para cada comparación usando umbral 2x fold change
genes_heatmap <- genes_heatmap %>%
  mutate(
    status_DKO_vs_WTG = case_when(
      avg_DKO >= 2 * avg_WTG ~ "aumentado",
      avg_DKO <= 0.5 * avg_WTG ~ "disminuido",
      TRUE ~ "no diferencial"
    ),
    status_KOS_vs_WTG = case_when(
      avg_KOS >= 2 * avg_WTG ~ "aumentado",
      avg_KOS <= 0.5 * avg_WTG ~ "disminuido",
      TRUE ~ "no diferencial"
    ),
    status_KOM_vs_WTG = case_when(
      avg_KOM >= 2 * avg_WTG ~ "aumentado",
      avg_KOM <= 0.5 * avg_WTG ~ "disminuido",
      TRUE ~ "no diferencial"
    ),
    status_DKO_vs_KOM = case_when(
      avg_DKO >= 2 * avg_KOM ~ "aumentado",
      avg_DKO <= 0.5 * avg_KOM ~ "disminuido",
      TRUE ~ "no diferencial"
    ),
    status_DKO_vs_KOS = case_when(
      avg_DKO >= 2 * avg_KOS ~ "aumentado",
      avg_DKO <= 0.5 * avg_KOS ~ "disminuido",
      TRUE ~ "no diferencial"
    )
  )

#exportar a excel 
resultado <- genes_heatmap %>%
  select(gene_name, avg_WTG, avg_WTP, avg_DKO, avg_KOP,
         status_DKO_vs_WTG, status_KOP_vs_WTP, status_DKO_vs_KOP, status_WTG_vs_WTP)

# Ruta donde quieres guardar el archivo
ruta_guardado <- "D:/Carpetas/Experimentos/Celulas Endoteliales/RNA/RNAseq/Resultados/Analisis nuevo/GCKIII/Genes_status_GCKIII.xlsx"

# Guardar resultado en la ruta especificada
write.xlsx(resultado, ruta_guardado, overwrite = TRUE)

cat("Archivo guardado en:", ruta_guardado, "\n")

